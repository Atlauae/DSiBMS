library(limma)
library(edgeR)
library(Homo.sapiens)
library(ggplot2)
library(stringr)

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
metadata[3,] <- str_replace(metadata[3,], 'tissue: ', '')
metadata[3,] <- str_replace(metadata[3,], 'corpus callosum', 'corpus_callosum')
metadata[3,] <- str_replace(metadata[3,], 'frontal cortex', 'frontal_cortex')
metadata[3,] <- str_replace(metadata[3,], 'internal capsule', 'internal_capsule')
metadata[3,] <- str_replace(metadata[3,], 'parietal cortex', 'parietal_cortex')

metadata[4,] <- str_replace(metadata[4,], 'disease state: ', '')
metadata[4,] <- str_replace(metadata[4,], 'healthy control', 'control')

#Load counts file
counts <- read.csv(paste(COUNTS_PATH, 'GSE123496_Human_MSNL_counts.csv', sep=''), header=TRUE)

#Update sample names of counts file
colnames(counts) <- metadata[1,]

#Filtering out lowly expressed genes
cpms	<- 1000000*counts/colSums(counts)

conditions <- levels(factor(metadata[4,]))
regions <- levels(factor(metadata[3,]))
highlyExpressed <- NULL

for(i in conditions) {
  sub1 <- as.vector(metadata[1, metadata[4,] == i])
  subset1 <- cpms[, colnames(cpms) == sub1]
  
  for(o in regions) {
    sub2 <- as.vector(metadata[1, metadata[3,] == o])
    subset2 <- subset1[, colnames(subset1) %in% sub2]
    highLowExpression <- rowSums(subset2>=1)>=round(ncol(subset2)/2)
    highExpression <- names(highLowExpression[!highLowExpression == FALSE])
    
    highlyExpressed <- c(highlyExpressed, highExpression)
  }
}

highlyExpressed <- highlyExpressed[!duplicated(highlyExpressed)]
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
  cc = CTRLcc-MScc, 
  fc = CTRLfc-MSfc, 
  hc = CTRLhc-MShc, 
  ic = CTRLic-MSic,
  pc = CTRLpc-MSpc,
  levels = colnames(design))

#Apply the statistical models
v <- voom(dgeList, design, plot=FALSE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)

lfc_cutoff = 1 
tfit <- treat(vfit, lfc=lfc_cutoff)

#Identify differentially expressed genes - p-value cutoff and logFoldChange cutoff
Pvalue_cutoff = 0.05
dt <- decideTests(tfit, p.value =Pvalue_cutoff)

#Make volcano plot
cc_Amean <- rowMeans(filteredCpms[,colnames(filteredCpms) == metadata[1, metadata[3,] == "corpus_callosum"]])
cc_plot <- data.frame(tfit$coefficients[,1], cc_Amean, tfit$p.value[,1], dt[,1])
colnames(cc_plot) <- c("logFC", "Amean", "p_value", "up_down")

amplot <- ggplot(cc_plot, aes(x=log(Amean,2), y=logFC, col=factor(up_down)))
amplot <- amplot + geom_point(alpha=0.4)
print(amplot)

volcplot <- ggplot(cc_plot, aes(x=logFC, y=-log(p_value,2), col=factor(up_down)))
volcplot <- volcplot + geom_point()
print(volcplot)
