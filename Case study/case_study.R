library(limma)
library(edgeR)
library(Homo.sapiens)
library(ggplot2)
library(stringr)
library(reshape2)
library(eulerr)
library(gridExtra)

if	(!requireNamespace("BiocManager",	quietly	=	TRUE))
  install.packages("BiocManager")
packages <- c("limma", "Homo.sapiens", "stringr")
BiocManager::install(packages)

#Set working directory (change to fit your own thing)
#If you remove the paste function and just put down getwd() it should return the folder your R file is in
WORK_DIR <- paste(getwd(),'/Case study', sep='') ##Change this to match your working directory
COUNTS_PATH <- paste(WORK_DIR,'/counts/', sep='')

#Create count folder in your working directory
dir.create(COUNTS_PATH)

#Get url from GEO page, under download (http) 
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123496&format=file&file=GSE123496%5FHuman%5FMSNL%5Fcounts%2Ecsv%2Egz"

#Download counts files - check "destfile" - adjust the file name as specified under: "Supplementary file"
destination_file <- paste(COUNTS_PATH, "GSE123496_Human_MSNL_counts.csv.gz", sep='')
utils::download.file(url, destfile= destination_file, mode="wb") 

#Unzip the 'gz' files
gzipped_files <- list.files(COUNTS_PATH, pattern='.gz', full.names = TRUE)
for(i in paste(gzipped_files)) {
  R.utils::gunzip(i, overwrite=TRUE)
}
  
#Get URL from - Series Matrix File(s)
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123496/matrix/GSE123496_series_matrix.txt.gz"
destination_file <- paste(WORK_DIR, "GSE123496_series_matrix.txt.gz", sep='/')
utils::download.file(url, destfile=destination_file, mode="wb") 
R.utils::gunzip(destination_file, overwrite=TRUE)

#Open the META-data file with excel and determine number of rows to skip
metadata <- as.data.frame(read.delim('Case study/GSE123496_series_matrix.txt', skip=30, sep='\t', header=FALSE))
##Make sure you check that the file location is correct, due to how I work it is somewhat different than standard

#Select useful information (rows) from the metadata
row_useful <- c(1,2,10,11)
metadata <- metadata[row_useful,2:ncol(metadata),]

#Remove uneccesary indicators from the metadata
rownames_metadata <- c("Sample", "GEO_accession", "Region", "Condition")
row.names(metadata) <- rownames_metadata

metadata[3,] <- str_replace(metadata[3,], 'tissue: ', '')
metadata[3,] <- str_replace(metadata[3,], 'corpus callosum', 'corpus_callosum')
metadata[3,] <- str_replace(metadata[3,], 'frontal cortex', 'frontal_cortex')
metadata[3,] <- str_replace(metadata[3,], 'internal capsule', 'internal_capsule')
metadata[3,] <- str_replace(metadata[3,], 'parietal cortex', 'parietal_cortex')

metadata[4,] <- str_replace(metadata[4,], 'disease state: ', '')
metadata[4,] <- str_replace(metadata[4,], 'healthy control', 'control')

#Making a transposed version of the metadata for easier use later
metadata_trans <- as.data.frame(t(metadata))

#Load counts file
counts <- read.csv(paste(COUNTS_PATH, 'GSE123496_Human_MSNL_counts.csv', sep=''), header=TRUE)

#Update sample names of counts file
colnames(counts) <- metadata[1,]

#Filtering out lowly expressed genes
cpms	<- 1000000*counts/colSums(counts)

conditions <- levels(factor(metadata_trans$Condition))
regions <- levels(factor(metadata_trans$Region))
highlyExpressed <- NULL

#Retrieve highly expressed genes per group as described in the lecture
for(i in conditions) {
  sub1 <- as.vector(metadata_trans$Sample[metadata_trans$Condition == i])
  subset1 <- cpms[, colnames(cpms) == sub1]
  
  for(o in regions) {
    sub2 <- as.vector(metadata_trans$Sample[metadata_trans$Region == o])
    subset2 <- subset1[, colnames(subset1) %in% sub2]
    highLowExpression <- rowSums(subset2>=1)>=round(ncol(subset2)/2)
    highExpression <- names(highLowExpression[!highLowExpression == FALSE])
    
    highlyExpressed <- c(highlyExpressed, highExpression)
  }
}

#Remove duplicates from generated list
highlyExpressed <- highlyExpressed[!duplicated(highlyExpressed)]

#Assign filtered counts and cpms
filteredCounts	<- counts[row.names(counts) %in% highlyExpressed,]
filteredCpms	<- cpms[row.names(cpms) %in% highlyExpressed,]

#Load filtered data into DGEList object
dgeList <- DGEList(count=filteredCounts)

#Assign metadata
brainregion = as.character(metadata[3,])
patientcondition = as.character(metadata[4,])

dgeList$samples$brainregion <- as.factor(brainregion)
dgeList$samples$patientcondition <- as.factor(patientcondition)

#Assign gene symbol
geneid <- rownames(dgeList)
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL"), keytype="ENSEMBL")
genes <- genes[!duplicated(genes$ENSEMBL),]
dgeList$genes <- genes

#Calculate normalization factors using the TMM method
dgeList <- calcNormFactors(dgeList, method = "TMM")



###DGE between conditions
#Define the model used for differential gene expression analysis total between different conditions
condition <- dgeList$samples$patientcondition
design_t <- model.matrix(~0 + condition)

#Make  contrast.matrix. e.g. (determine comparisons of interest)
contr.matrix_t <- makeContrasts(
  condition_comp = conditioncontrol-conditionMS, 
  levels = colnames(design_t))

#Apply the statistical models
v_t <- voom(dgeList, design_t, plot=FALSE)
vfit_t <- lmFit(v_t, design_t)
vfit_t <- contrasts.fit(vfit_t, contrasts=contr.matrix_t)

lfc_cutoff = 1 
tfit_t <- treat(vfit_t, lfc=lfc_cutoff)

#Identify differentially expressed genes - p-value cutoff and logFoldChange cutoff
Pvalue_cutoff = 0.05
dt_t <- decideTests(tfit_t, p.value =Pvalue_cutoff)

Amean_t <- rowMeans(filteredCpms)
symbols_t <- dgeList$genes$SYMBOL[dgeList$genes$ENSEMBL %in% names(tfit_t$coefficients[,1])]
plot_t <- data.frame(tfit_t$coefficients[,1], Amean_t, tfit_t$p.value[,1], p.adjust(tfit_t$p.value[,1], method="BH"), dt_t[,1], symbols_t)
colnames(plot_t) <- c("logFC", "Amean", "p_value", "adj_p_value", "up_down", "symbol")

plot_t$symbol[is.na(plot_t$symbol)] <- row.names(plot_t[is.na(plot_t$symbol),])

volcplot_t <- ggplot(plot_t, aes(x=logFC, y=-log(p_value,2), col=factor(up_down)))
volcplot_t <- volcplot_t + geom_point()
print(volcplot_t)



###DGE between brain regions
#Define the model used for differential gene expression analysis of a brain regions between different conditions
regions <- dgeList$samples$brainregion
condition <- dgeList$samples$patientcondition
design <- model.matrix(~0 + condition:regions)

#Change design colnames
colnames_design <- c('CTRLcc', 'MScc',
                     'CTRLfc', 'MSfc',
                     'CTRLhc', 'MShc',
                     'CTRLic', 'MSic',
                     'CTRLpc', 'MSpc')
colnames(design) <- colnames_design

#Make  contrast.matrix. e.g. (determine comparisons of interest)
contr.matrix <- makeContrasts(
  cc = MScc-CTRLcc, 
  fc = MSfc-CTRLfc, 
  hc = MShc-CTRLhc, 
  ic = MSic-CTRLic,
  pc = MSpc-CTRLpc,
  levels = colnames(design))

#Apply the statistical models
v <- voom(dgeList, design, plot=FALSE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)

lfc_cutoff = 1 
tfit <- treat(vfit, lfc=lfc_cutoff)

#Identify differentially expressed genes - p-value cutoff and logFoldChange cutoff
Pvalue_cutoff = 0.05
dt <- decideTests(tfit, p.value = Pvalue_cutoff)

#Make volcano plot
contr_names <- colnames(contr.matrix)
region_names <- levels(factor(metadata_trans$Region))

for(i in 1:5){
  Amean <- rowMeans(filteredCpms[,colnames(filteredCpms) %in% metadata_trans$Sample[metadata_trans$Region == region_names[i]]])
  symbols <- dgeList$genes$SYMBOL[dgeList$genes$ENSEMBL %in% names(tfit$coefficients[,i])]
  plot <- data.frame(tfit$coefficients[,i], Amean, tfit$p.value[,i], p.adjust(tfit$p.value[,i], method="BH"), dt[,i], symbols)
  colnames(plot) <- c("logFC", "Amean", "p_value", "adj_p_value", "up_down", "symbol")
  
  plot$symbol[is.na(plot$symbol)] <- row.names(plot[is.na(plot$symbol),])
  
  assign(paste(contr_names[i], "_plot", sep=""), plot)
}

amplot <- ggplot(fc_plot, aes(x=log(Amean,2), y=logFC, col=factor(up_down)))
amplot <- amplot + geom_point(alpha=0.4)
print(amplot)

volcplot <- ggplot(fc_plot, aes(x=logFC, y=-log(p_value,2), col=factor(up_down)))
volcplot <- volcplot + geom_point()
print(volcplot)



##Generate dataframe containing information for the heatmap
dif_gen_list <- NULL
regions <- levels(factor(metadata_trans$Region))

for(i in 1:5) {
  dif_gen <- row.names(subset(dt, dt[,i] == -1 | dt[,i] == 1))
  dif_gen_counts <- dgeList$counts[row.names(dgeList$counts) %in% dif_gen,]
  dif_gen_counts <- dif_gen_counts[,colnames(dif_gen_counts) %in% metadata_trans$Sample[metadata_trans$Region == regions[i]]]
  assign(paste("dif_gen_", contr_names[i],  sep=""), dif_gen_counts)
}
dif_gen_list <- list(dif_gen_cc, dif_gen_fc, dif_gen_hc, dif_gen_ic, dif_gen_pc)

for(i in 1:5){
  heatmap_symbols <- dgeList$genes$SYMBOL[dgeList$genes$ENSEMBL %in% row.names(dif_gen_list[[i]])]
  heatmap_symbols[is.na(heatmap_symbols)] <- row.names(dif_gen_list[[i]][is.na(heatmap_symbols),])
  heatmap_symbols[duplicated(heatmap_symbols)] <- row.names(dif_gen_list[[i]])[duplicated(heatmap_symbols)]
  dif_gen_list[[i]] <- cbind(dif_gen_list[[i]], heatmap_symbols)
}

test <- melt(dif_gen_list[[1]], variable = "heatmap_symbols")
colnames(test) <- c("Y", "X", "cpm")
ggplot(test, aes(x = X, y=Y, fill=cpm)) + geom_tile()


##Generate objects and .txt files containing the up and downregulated genes per region for later analysis and upload to metascape
for(i in 1:5){
  dif_gen_up_ensembl <- row.names(tfit[dt[,i]==1,])
  dif_gen_up_logfc <- tfit$coefficients[dt[,i]==1,i]
  
  dif_gen_down_ensembl <- row.names(tfit[dt[,i]==-1,])
  dif_gen_down_logfc <- tfit$coefficients[dt[,i]==-1,i]
  
  write.table(data.frame(Gene = dif_gen_up_ensembl, logFC = dif_gen_up_logfc), 
              file = paste(WORK_DIR, "/genelist_up_", contr_names[i], ".txt", sep=""), sep = "\t", row.names = FALSE)
  
  write.table(data.frame(Gene = dif_gen_down_ensembl, logFC = dif_gen_down_logfc), 
              file = paste(WORK_DIR, "/genelist_down_", contr_names[i], ".txt", sep=""), sep = "\t", row.names = FALSE)
}



##Generate venn diagram for the different groups
png(file = "venn_diagram_up.png", width = 10, height = 10, units = "in", res=300, pointsize = 12, bg = "white")	
vennDiagram(dt, circle.col=c("#3498DB", "#E74C3C", "#2ECC71", "#FFA500", "#1ABC9C"),include = c("up"), 
            main="Number of significantly upregulated genes different groups")
dev.off()

png(file = "venn_diagram_down.png", width = 10, height = 10, units = "in", res=300, pointsize = 12, bg = "white")	
vennDiagram(dt[,c(1,2,3,4,5)], circle.col=c("#3498DB", "#E74C3C", "#2ECC71", "#FFA500", "#1ABC9C"),include = c("down"), 
                main="Number of significantly downregulated genes different groups")
dev.off()
