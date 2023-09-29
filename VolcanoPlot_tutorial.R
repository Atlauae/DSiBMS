"
Tutorial: Visualization of Differential Expression of Genes
 using the ggplot2 library

Expression of many proteins and almost all expression levels of genes can be measured in any living species. 
Typical expriments to study expression levels in different condition, cell lines, mutations, etc.  
is to compared it (the Target) against a control.
This will result in a ratio value for each gene or protein: (level in Target) / (level in Control)

When differential protein or gene expression analysis is done on proteome or transcriptome data, 
typical results will contain at least the following values:

1) Gene or protein identifiers: ID (or sometimes called 'Key')
2) The ratio of expression values of the ID between Control (C) and Target (T): Fold value = T/C
3) The probability that the ID is not different in the Control and Target is described as a p-value ( 1 means that C and T are equal)
4) The average of the absolute expression value of the ID  (C + T)/2

Note: *Factors are used to describe replicates


Data Science in Biomedicine, 2023

Anne de Jong

"








##############################################################################################
##                        1. set working directory                                          ##
##############################################################################################

" 
The working directory is the folder on your PC or laptop where all the data is read and stored
NOTE: In R the syntax for folder structures is Linux compatible, meaning that subfolders 
are separated by forward slashes, compared to backslashes in windows. In R (R-studio) backslash 
is an escape character and to use it you need to write it as double back slash.
exmample Linux /data/user/tutorial
example windows C:\\data\\user\\tutorial

==> Escape Characters are common in programming languages, for more info see: 
https://www.w3schools.com/python/gloss_python_escape_characters.asp


==> use Ctrl-Enter to execute a line in R-Studio at the cursor position



==> FIRST STEPs: 
  1) goto: http://ngs.molgenrug.nl/courses/DataScienceinBiomedicine/VolcanoPCA_tutorial/
  2) Download data files you need to your laptop folder
  3) Adapt the line below to your meet your working folder"

setwd("C:\\Users\\Lenovo\\Documents\\github repositories\\DSiBMS")



"Show the content of the current Working Directory"
dir()



##############################################################################################
##                        2. Load data                                                      ##
##############################################################################################

"
The data table that we are going to load is a tab delimited file with column names (the experiments) 
in the first row and row names (the genes) in the first column.
The R function read.table can be used to load the data in the variable count_data.

Select below one of the two data sets
Be aware to load the proper associated Factor file. later in the script

"

url="http://ngs.molgenrug.nl/courses/DataScienceinBiomedicine/VolcanoPCA_tutorial/"


"The data"
DGE_Filename <- "trex2.Significant_Changed_Genes.1.A_F71Y-WT.txt"
DGE_url <- paste(url, DGE_Filename, sep="")
DGE_data <- read.table(DGE_url, header = TRUE, row.names = 1)

View(DGE_data)

"The Columns;
logFC       Log2 of the Fold change
logCPM      Log2 of the expression value
LR          -log2(likelyhood ratio): The higher the value the more significant
pvalue      pvalue (of a t-test)
adj-pvalue  pvalue corrected for multiple testing
Fold        Fold Change
minFDR      -log2(adj-pvalue)
"





##############################################################################################
##                         3. Loading libraries                                             ##
##############################################################################################

"
As you learned the library ggplot2 is a powerful method to create nice plots
Let start making a nice PCA plot on the basis of this library
"


"If the ggplot2 is not loaded, means that it's NOT available yet on your computer.
In this case you need to install in once and can be done using the BiocManager
as is shown below.
Alternatively you can use R-Studio menu: Tools->Install packages
"

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("ggplot2"))     { BiocManager::install("ggplot2",     ask=FALSE, dependencies = TRUE)}
if(!require("statmod"))     { BiocManager::install("statmod",     ask=FALSE, dependencies = TRUE)}

"First we need to load the ggplot2 library"
library(ggplot2)



##############################################################################################
##                         4. Visualizing data                                              ##
##############################################################################################


"
Add a theme for your Volcano plots.
A theme is optional, but very powerful to build your own style of plots
"
theme_Volcano_plot = theme(
  axis.text.x =	element_text(size = 20, hjust = 1), 
  axis.text.y = 	element_text(size = 20, hjust = 1),
  title = 		element_text(size = 36, face = "bold", hjust = 0.5),
  axis.title.x = 	element_text(size = 26, hjust = 0.5),
  axis.title.y = 	element_text(size = 26, hjust = 0.5)
)


" The basic Plot"
VolcanoPlot <- ggplot(DGE_data,  aes(x=logFC,y=minFDR)) + geom_point(pch = 18, cex = 4) 
print(VolcanoPlot)




##############################################################################################
##                         5. Adding information to the plot                                ##
##############################################################################################

"
Here we will indicate some areas in the plot

First we need to know the borders of the Volcano Plot, omitting missing values
na.rm=T  means Not Available (na) Remove is true (rm=T)
Using na.rm=T is very common in R and enables to deal with unmeasured sample and incomplete tables
"
xmin <- min(DGE_data$logFC, na.rm=T)
xmax <- max(DGE_data$logFC, na.rm=T)
ymin <- min(DGE_data$minFDR, na.rm=T)
ymax <- max(DGE_data$minFDR, na.rm=T)

"
Decide which threshold values we will use
FDR meains False Discovery Rate
"
FDR_threshold <- -log2(0.05)
Fold_threshold <- log2(2)

"Here we add grey rectangles to show non-significant areas"
VolcanoPlot <- VolcanoPlot + annotate("rect", xmin = xmin, xmax = xmax, ymin = 0, ymax = FDR_threshold, alpha = .2)
VolcanoPlot <- VolcanoPlot + annotate("rect", xmin = -Fold_threshold, xmax = Fold_threshold, ymin = FDR_threshold, ymax = ymax, alpha = .2)
VolcanoPlot

"Add a title"
VolcanoPlot <- VolcanoPlot + ggtitle("My Nice Volcano Plot")
VolcanoPlot

"To prettify the plot you can apply your Theme"
VolcanoPlot <- VolcanoPlot + theme_Volcano_plot
VolcanoPlot





##############################################################################################
##                         6. Adding expression values to the plot                          ##
##############################################################################################

"
Volcano Plots are very popular but expression values are not included in the plot
A simple solution is to draw a transparent (alpha= 1/4) circle reflecting the 
expression value: a sphere around each dot
"

DGE_data$sphere <- ifelse(DGE_data$logCPM<1, 1, DGE_data$logCPM)
VolcanoPlot <- VolcanoPlot + geom_point(size = DGE_data$sphere, alpha = 1/4)
VolcanoPlot




##############################################################################################
##                         7. Adding Classification                                         ##
##############################################################################################


"
Let's say you want to highlight certain genes/proteins in a plot.
"

"Read the genes of interest from a file
This file is a tab delimited table that contain the ID's of proteins or genes, 
a color, and a group name.
Because the ID's are divided into groups (Classes) we call this Classification of genes/proteins
"
Class_Filename <- "trex2.Class.txt"
Class_url <- paste(url, Class_Filename, sep = "")
my_class_genes<- read.table(Class_url, header = TRUE, row.names = 1)
head(my_class_genes)


"
Merge the Classification data with the DGE table
"
DGE_Table <- merge(DGE_data, my_class_genes, by="row.names",all.x=TRUE)
colnames(DGE_Table)[colnames(DGE_Table)=="GeneID.x"] <- "GeneID"
head(DGE_Table)

"
To plot this we can add a color factor on the basis of the Groups (Classes) to which genes belong
"
VolcanoPlot <- ggplot(DGE_Table,  aes(x=logFC,y=minFDR,color=factor(Group))) + geom_point(pch = 18, cex = 4) 
print(VolcanoPlot)



"==> At this point you again can add the gray bars, if you want"

"Here we add grey rectangles to show the non-significant areas"
VolcanoPlot <- VolcanoPlot + annotate("rect", xmin = xmin, xmax = xmax, ymin = 0, ymax = FDR_threshold, alpha = .2)
VolcanoPlot <- VolcanoPlot + annotate("rect", xmin = -Fold_threshold, xmax = Fold_threshold, ymin = FDR_threshold, ymax = ymax, alpha = .2)
print(VolcanoPlot)

"Add a title"
VolcanoPlot <- VolcanoPlot + ggtitle("My Nice Volcano Plot")
print(VolcanoPlot)

"To prettify the plot you can apply your Theme"
VolcanoPlot <- VolcanoPlot + theme_Volcano_plot
print(VolcanoPlot)

"Add labels"
# 1. make a subset of IDs you want to annotate. E.g., all genes with logFC > 2 and minFDR (aka -log2(pvalue)) > 50
MyLabels=subset(DGE_Table, abs(logFC) > 2 & minFDR > 50)
# 2. Check the labels you selected
head(MyLabels)
# 2. Add this subset to your plot
VolcanoPlot <- VolcanoPlot + geom_text(data=MyLabels, aes(label=Row.names),hjust=0, vjust=0)
print(VolcanoPlot)



##############################################################################################
##                         8. Saving Graphics                                               ##
##############################################################################################

"Once the plot is stored in a variable, in our case 'VolcanoPlot', you want to save it
as jpg, png or pdf"

"
Lets save it as a png of 10x10 inch with a resolution of 300dpi (dots per inch)
"
png(file = "My_super_Nice_VolcanoPlot.png", width = 10, height = 10, units = "in", res=300, pointsize = 12, bg = "white")	
print(VolcanoPlot)
dev.off()



##############################################################################################
##                         9. Assignment                                                    ##
##############################################################################################

"
Available contrast the files: T2-T1.txt, T3-T1.txt, T4-T1.txt, T5-T1.txt, T6-T1.txt

1) Choose one and make a Volcano Plot
2) Also Make an MA plot: Expression x Ratio ( see lectures)
3) Annotate (add names) some dots
  Hint: 
4) Optional: Make a multiplot
  hint: https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
  or more sophisticated: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

Export the result as one pdf and upload it in Brightspace
One per Couple!!

"

#Retreive data for exercise, using T4-T1.txt
url="http://ngs.molgenrug.nl/courses/DataScienceinBiomedicine/VolcanoPCA_tutorial/"
DGE_Filename <- "T4-T1.txt"
DGE_url <- paste(url, DGE_Filename, sep="")
DGE_data <- read.table(DGE_url, header = TRUE, row.names = 1)

#Retrieve genes of interest for classification of genes (using the class file of the example volcano plot)
Class_Filename <- "trex2.Class.txt"
Class_url <- paste(url, Class_Filename, sep = "")
my_class_genes<- read.table(Class_url, header = TRUE, row.names = 1)
head(my_class_genes)

#Merge the Classification data with the DGE table
DGE_Table <- merge(DGE_data, my_class_genes, by="row.names",all.x=TRUE)
colnames(DGE_Table)[colnames(DGE_Table)=="GeneID.x"] <- "GeneID"

View(DGE_Table)


#Create custom theme for the plot



#Create basic VolcanaPlot for initial data investigation
VolcanoPlot <- ggplot(DGE_data,  aes(x=logFC,y=minFDR)) + geom_point(pch = 16, cex = 2)
print(VolcanoPlot)

#Create Volcano plot with added color factor based on gene class
VolcanoPlot <- ggplot(DGE_Table,  aes(x=logFC,y=minFDR,color=factor(Group))) + geom_point(pch = 16, cex = 2) 


##Set color palatte to match colors given in dataset
#Retrieve colors in correct order from the dataset
dataset_color <- unique(DGE_Table$Color[order(factor(DGE_Table$Group))])

#Define color for non-specified genes
dataset_color[is.na(dataset_color)] <- "black"

#Set dataset_color as the color palatte to be used in the plot
VolcanoPlot <- VolcanoPlot + scale_color_manual(values=dataset_color)


##Adding informative visual information the plot
#Define axis values
xmin <- min(DGE_Table$logFC, na.rm=T)
xmax <- max(DGE_Table$logFC, na.rm=T)
ymin <- min(DGE_Table$minFDR, na.rm=T)
ymax <- max(DGE_Table$minFDR, na.rm=T)

#Decide which threshold values to use
FDR_threshold <- -log2(0.05)
Fold_threshold <- log2(2)

#Add in gray areas indicating insignificant results and dashed lines bordering these areas
VolcanoPlot <- VolcanoPlot + annotate("rect", xmin = xmin, xmax = xmax, ymin = 0, ymax = FDR_threshold, alpha = .2)
VolcanoPlot <- VolcanoPlot + annotate("rect", xmin = -Fold_threshold, xmax = Fold_threshold, ymin = FDR_threshold, ymax = ymax, alpha = .2)
VolcanoPlot <- VolcanoPlot + geom_segment(x = xmin, xend = xmax, y = FDR_threshold, yend = FDR_threshold, linetype="dashed")
VolcanoPlot <- VolcanoPlot + geom_segment(x = -Fold_threshold, xend = -Fold_threshold, y = 0, yend = ymax, linetype="dashed")
VolcanoPlot <- VolcanoPlot + geom_segment(x = Fold_threshold, xend = Fold_threshold, y = 0, yend = ymax, linetype="dashed")

#Add title to the plot
VolcanoPlot <- VolcanoPlot + ggtitle("My Nice Volcano Plot")

#Add labels to the axes
VolcanoPlot <- VolcanoPlot + labs(x = "-log2(adj. p-value)", y = "log2(Fold Change)")

#Change legend title
VolcanoPlot <- VolcanoPlot + guides(color = guide_legend(title = "Gene classes of interest"))

#Add actual expression indication
DGE_Table$sphere <- ifelse(DGE_Table$logCPM<1, 1, DGE_Table$logCPM)
VolcanoPlot <- VolcanoPlot + geom_point(size = DGE_Table$sphere, alpha = 1/8)

#Remove legend from plot and apply classis theme
VolcanoPlot <- VolcanoPlot + theme_classic() + theme(legend.position="none")

print(VolcanoPlot)



### Create an MA plot for the dataset ###

#Create Volcano plot with added color factor based on gene class
MAPlot <- ggplot(DGE_Table,  aes(x=logCPM,y=LR,color=factor(Group))) + geom_point(shape = 16, cex = 2) 


##Set color palatte to match colors given in dataset

#Set dataset_color as the color palatte to be used in the plot
MAPlot <- MAPlot + scale_color_manual(values=dataset_color)


##Adding informative visual information the plot

#Define axis values
xmin_MA <- min(DGE_Table$logCPM, na.rm=T)
xmax_MA <- max(DGE_Table$logCPM, na.rm=T)
ymin_MA <- min(DGE_Table$LR, na.rm=T)
ymax_MA <- max(DGE_Table$LR, na.rm=T)

#Decide which threshold values to use
LR_threshold <- -log2(2)

#Add in gray areas indicating insignificant results and dashed lines bordering these areas
MAPlot <- MAPlot + annotate("rect", xmin = xmin_MA, xmax = xmax_MA, ymin = LR_threshold, ymax = -LR_threshold, alpha = .2)
MAPlot <- MAPlot + geom_segment(x = xmin_MA, xend = xmax_MA, y = LR_threshold, yend = LR_threshold, linetype="dashed")
MAPlot <- MAPlot + geom_segment(x = xmin_MA, xend = xmax_MA, y = -LR_threshold, yend = -LR_threshold, linetype="dashed")

#Add title to the plot
MAPlot <- MAPlot + ggtitle("My Nice MA Plot")

#Add labels to the axes
MAPlot <- MAPlot + labs(x = "log2(CPM)", y = "-log2(LR)")

#Apply classis theme
MAPlot <- MAPlot + theme_classic()

print(MAPlot)

### Plot MA and Volcano plot together ###

#Load required package
library(grid)
library(gtable)

#Combine the plots into one
g_MAPlot <- ggplotGrob(MAPlot)
g_VolcanoPlot <- ggplotGrob(VolcanoPlot)
g_MAPlot$widths <- unit.pmax(g_MAPlot$widths, g_VolcanoPlot$widths)
g_VolcanoPlot$widths <- unit.pmax(g_MAPlot$widths, g_VolcanoPlot$widths)
combined_plot <- rbind(g_MAPlot, g_VolcanoPlot, size = "first")


#Generate combined plot
grid.newpage()
grid.draw(combined_plot)

#Save combined plot as image
png(file = "Combined_MA_Plot_Volcano_Plot.png", width = 10, height = 10, units = "in", res=300, pointsize = 12, bg = "white")	
grid.newpage()
grid.draw(combined_plot)
dev.off()