###############################################################
##                        Part 1                             ##
###############################################################

##1. Loading the data
#Set working dir
setwd("C:\\Users\\Lenovo\\Documents\\github repositories\\DSiBMS\\DGE tut")

#Load Data
annotation	<- read.delim('white_samples.txt',	header=TRUE)
readCounts	<- read.delim('white_readcounts.txt',	header=TRUE, 
                         row.names='Gene')

#How many genes and samples do we have?
dim(readCounts)


##2. Check read depth per sample
librarySizes	<- colSums(readCounts)
barplot(librarySizes,	
        names=colnames(readCounts),	
        las=2,
        main="Library	sizes")


##3. Filter out lowly expressed genes using option 2
#	Filtering	option	2:	Remove	lowly	expressed
#	those	that	do	not	show	100 counts-per-million	(CPM)
#	expression	level	in	3	or	more	samples
cpms	<-1000000*readCounts/colSums(readCounts)
highlyExpressed	<- rowSums(cpms>=100)>=3
readCounts	<- readCounts[highlyExpressed,]
cpms	<- cpms[highlyExpressed,]
dim(readCounts)


##4. Build a principal component map
#Transpose matrix and scale the CPM values to log space
pcaDat	<- prcomp(log(t(cpms)))

#check	the	proportion	of	variance	explained	by	principal	components	PC1,	PC2,	PC3	…
varExplained	<- pcaDat$sdev^2/sum(pcaDat$sdev^2)
varExplained[1:5]
barplot(varExplained)

#Plot PCA
plot(pcaDat$x, pch=16,	cex=0.5)
text(pcaDat$x,	labels=annotation$Name,	cex=0.7)




###############################################################
##                        Part 2                             ##
###############################################################

##1. Install and load edgeR library
#Install	edgeR	library	if	you	did	not	install	it	yet
if	(!requireNamespace("BiocManager",	quietly	=	TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

#Load edgeR
library(edgeR)


##2. Loading the data
#Code is given in brightspace file, this has however already been done within this file and will therefore not be repeated


##3. Filtering	out	lowly	expressed	genes
#Code is given in brightspace file, this has however already been done within this file but is changed here
cpms	<- 1000000*readCounts/colSums(readCounts)
highlyExpressed	<- rowSums(cpms>=1)>=3
filteredCounts	<- readCounts[highlyExpressed,]
filteredCpms	<- cpms[highlyExpressed,]
dim(filteredCounts)


##4. Perform differential expression analysis
#	Define	sample	groups	according	to	annotation	table	(young	or	old)
group	<- annotation$Type

#	Load	egdeR	library	and	create	an	object	containing	read	count	data
dgeList	<- DGEList(count=filteredCounts,	group=group)

#	Calculate	normalization	factors	using	the	TMM	method
dgeList	<- calcNormFactors(dgeList)

#Check calculated normalization	factors	by	querying	the	data	inside	this object
dgeList$samples$norm.factors

#	Define	the	model	used	for	differential	gene	expression	analysis	
design	<- model.matrix(~group)

#	Estimate	dispersions
dgeList	<- estimateDisp(dgeList,	design)

#	Fitting	the	model
fit	<- glmQLFit(dgeList,	design)

#	Perform	statistical	testing	using	quasi-likelihood	F-test
qlfTest	<- glmQLFTest(fit,	coef=2)

#The same command as before has been given, but modified to be more explicit IT DOES NOT DO ANYTHING DIFFERENT AT ALL. See below:
qlfTest	<- glmQLFTest(fit,	coef="groupYoung")

#Display	the	10	most	significant	genes
topTags(qlfTest)

#	Save	all	results	into	a	tab-delimited	file
dgeResult	<- topTags(qlfTest,	n='all')
write.csv(dgeResult,	file='white_old_vs_young.txt')


#Read the saved file and load it into the variable dgeTable
dgeTable	<- read.csv('white_old_vs_young.txt')

#See how many differentially expressed genes we can observe
significantGenes	<- dgeTable[dgeTable$FDR<=0.01,]
upGenes	<- significantGenes[significantGenes$logFC>0,]
downGenes	<- significantGenes[significantGenes$logFC<0,]


##5. Visualizing	results	of	differential	gene	expression	analysis
#	Draw	volcano plot
plot(dgeTable$logFC,	-log10(dgeTable$FDR),	pch=16,	cex=0.4,	col='black')
points(upGenes$logFC,	-log10(upGenes$FDR),	pch=16,	cex=0.5,	col='red')
points(downGenes$logFC,	-log10(downGenes$FDR),	pch=16,	cex=0.5,	
       col='blue')

# MA	plot
plot(dgeTable$logCPM,	dgeTable$logFC,	pch=16,	cex=0.4,	col='black')
points(upGenes$logCPM,	upGenes$logFC,	pch=16,	cex=0.5,	col='red')
points(downGenes$logCPM,	downGenes$logFC,	pch=16,	cex=0.5,	col='blue')

#Extract data for heatmap
significantCpms<-filteredCpms[rownames(filteredCpms)	%in% significantGenes$X,]

#Define a	blue-white-red	palette	with	300	different	colors
palette<-colorRampPalette(c('blue','white','red'))(300)

#Draw the heatmap
heatmap(as.matrix(significantCpms), scale	=	'row',	col=palette,	cexRow=0.2,	
        cexCol=0.7)
