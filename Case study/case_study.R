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
metadata[4,] <- str_replace(metadata[4,], 'disease state: ', '')

#Load counts file
counts <- read.csv(paste(COUNTS_PATH, 'GSE123496_Human_MSNL_counts.csv', sep=''), header=TRUE)

#Filtering out lowly expressed genes
cpms	<- 1000000*counts/colSums(counts)
highlyExpressed	<- rowSums(cpms>=1)>=25
filteredCounts	<- counts[highlyExpressed,]
filteredCpms	<- cpms[highlyExpressed,]

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