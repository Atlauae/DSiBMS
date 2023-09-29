"
Tutorial for Principle Component Analysis (PCA) and visualization using the ggplot2 library
For background on PCA, see: Lecture - 04 - Principle Component Analysis

Data Science in Biomedicine 2023

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

setwd("G:\\My Drive\\WERK\\Cursus_DataScienceInBiomedicine\\VolcanoPCA_tutorial")



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


"==> DATASET 1: Load data from your working folder"
count_data <- read.csv("TimeSeries_168_Counts.txt", sep="\t",header = TRUE, row.names = 1)

"or direct from internet 
NOTE: merging two string in R is done via the function paste()  "
filename=paste(url,"TimeSeries_168_Counts.txt",sep="")
count_data <- read.csv(filename, sep="\t", header = TRUE, row.names = 1)

# Check if the data is properly loaded:
  # Method 1:
    head(count_data)
  # Method 2:
    # In R-studio you have a window with the name "Environment" in the Upper-Right panel
    # search for the variable count_data and click it to view the content
    # NOTE: after viewing the content go back to this screen "PCA_tutorial
  # Method 3:
    View(count_data)


    
"=================================================================================="    
    
"==> DATASET 2: CodY_Counts.txt"
  # Adapt the code above to load DATASET 2 when needed for the second exercise 

"=================================================================================="  




##############################################################################################
##                        3. handling zero values                                           ##
##############################################################################################

"
Zero values always need your attention. 
Zeros can represent missing values and need to be discarded
But zeros can be real values, such as non-expressed genes.
In our case we want to keep zero values. To avoid problems in statistical R routines 
zero values can be increased to a value just above zero using a fixed arbitrary value. This is an easy method, 
but alternatively you could calculate the level of the background noise and use this value as lower 
limit.
"

"Here we choose the easy method by adding +1 to all values
In R this is very simple using the line below
"
count_data <- count_data + 1





##############################################################################################
##                        4. Formatting data for the PCA                                    ##
##############################################################################################

"
PCA is a standard function in R (with many options)
We will perform the PCA on log2 transformed data.
But note that the log function in R, by default, will take natural log of the value
In proteomics and transcriptomics data is always log2 transformed

Important is to know that the PCA function in R needs the genes in columns and experiments in rows
We can simply transform the table using the t function in R
"

"Log2 transfor the data"
log_count_data <- log(count_data,2)
head(log_count_data)

"Transform the data table, NOTE: have a look at the data using R-studio Environment"
pca_data <- t(log_count_data)





##############################################################################################
##                        5. Perform the PCA                                                ##
##############################################################################################

"
Here we use the prcomp function in R to perform the PCA
https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp
"

experiments.pca <- prcomp(pca_data, center = TRUE, scale = FALSE) 


" The function 'summary' will show the contribution of the Principle Components to the variance "
summary(experiments.pca)

"The content of the data structure experiments.pca"
head(experiments.pca$rotation)
head(experiments.pca$center)
head(experiments.pca$x)


"If we plot the contribution of the Princple Components to the variance 
we can observe that the first 2 PC's contribute the most to the Variance"
plot(experiments.pca, type = "l")

"We can plot the first 2 Principle Components for all experiments using the plot function of R"
plot(experiments.pca$x, type = "p")



##############################################################################################
##                         6. Loading libraries                                             ##
##############################################################################################

"
As you learned, library ggplot2 is a powerful method to create nice plots
Let start making a nice PCA plot on the basis of this library
"

"If the ggplot2 is not present in your R installation, you need to install in once.
Most simple method is to use BiocManager as is shown below.
See also: https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html
Alternatively you can use R-Studio menu: Tools->Install packages
"

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("ggplot2"))     { BiocManager::install("ggplot2",     ask=FALSE, dependencies = TRUE)}
if(!require("statmod"))     { BiocManager::install("statmod",     ask=FALSE, dependencies = TRUE)}



"To be able to use ggplot2 in R you need to load the library"
library(ggplot2)


##############################################################################################
##                         7. Visualizing data                                              ##
##############################################################################################

"
Add a theme for your PCA plots.
A theme is optional, but very powerful to build your own style of graphs
"
theme_PCA_plot = theme(
  axis.text.x =	element_text(size = 20, hjust = 1), 
  axis.text.y = 	element_text(size = 20, hjust = 1),
  title = 		element_text(size = 36, face = "bold", hjust = 0.5),
  axis.title.x = 	element_text(size = 26, hjust = 0.5),
  axis.title.y = 	element_text(size = 26, hjust = 0.5)
)



# 7a. Get the data from the PCAnalysis
"
get the experiments x PC data from experiments.pca
"
my_pca_data <- as.data.frame(experiments.pca$x)





# 7b. We can plot all experiments but we can also decide to combine experiments that are replicates

" Add a factors to describe the replicates and merge this with the PCA data
In our example we use one Factor (Time Points) but often Multiple Factors are used (e.g. Time and Treatments) 
MF : Multiple Factors
" 

"For DATASET 1 the MF is:"
MF_Filename="TimeSeries_168_Factors.txt"


# OPTION 1 load from file or url
MF_url=paste(url, MF_Filename, sep="")
MF <- read.table(MF_url, header=T)

# OPTION2 if you downloaded the file: 
MF <- read.table(MF_Filename, header=T)

my_pca_data$experiment <- row.names(my_pca_data)
my_pca_data <-  merge(my_pca_data, MF, by = "experiment")
head(my_pca_data)


"For DATASET 2 the MF is:"
# MF_Filename="CodY_Factors.txt"


# 7c Let's make nice plots 

# Plot PC1 x PC2 of Time Series: DATASET 1
PCAplot <- ggplot(my_pca_data,  aes(x=PC2,y=PC3, label=factor(Time), color = factor(Time))) + geom_point(pch = 18, cex = 4) 
PCAplot <- PCAplot + geom_text(hjust = -0.4, size = 3, check_overlap = TRUE, aes(colour = factor(Time)))
PCAplot <- PCAplot + ggtitle("PCA plot")
print(PCAplot)

# Plot PC1 x PC2 of CodY: DATASET 2
PCAplot <- ggplot(my_pca_data,  aes(x=PC2,y=PC3, label=factor(strain), color = factor(strain))) + geom_point(pch = 18, cex = 4) 
PCAplot <- PCAplot + geom_text(hjust = -0.4, size = 3, check_overlap = TRUE, aes(colour = factor(strain)))
PCAplot <- PCAplot + ggtitle("PCA plot")
print(PCAplot)




##############################################################################################
##                         8. Saving Graphics                                               ##
##############################################################################################

"Once the graph is stored in a variable, in our case 'PCAplot', you can save the images
as jpg, png or pdf"

"let's save it as a png of 10x10 inch with a resolution of 300dpi (dots per inch)"
png(file = "My_super_Nice_PCA_plot.png", width = 10, height = 10, units = "in", res=300, pointsize = 12, bg = "white")	
print(PCAplot)
dev.off()


# Assignment:
# Make 2 nice PCA plots; one of data set 1 and one on the basis of data set 2
# Upload these 2 plots as couple assignment




##############################################################################################
##                         9. PCA plots of genes                                            ##
##############################################################################################



"
1) Perform the PCA on genes using the prcomp function of R
To do this we need to rotate the table using the t function of R
MORE about transpose: https://www.statology.org/transpose-data-frame-in-r/
"
genes.pca <- prcomp(t(pca_data), center = TRUE, scale = TRUE) 
summary(genes.pca)


"
2) Get the PC (Principle Components) values from the prcomp object. The name of this object is genes.pca
"
my_pca_data <- as.data.frame(genes.pca$x)



"
3) As seen before we classify genes
"

CLASS_filename = "trex2.Class.txt"
CLASS_url = paste(url,CLASS_filename, sep="")

my_class_genes <- read.csv(CLASS_url,sep="\t")
colnames(my_class_genes) =c("GeneID","color","Group")
row.names(my_class_genes) <- my_class_genes$GeneID
head(my_class_genes)

"
4) Add the Classification data to the PCA result
"
my_pca_data <- merge(my_pca_data, my_class_genes, by="row.names",all.x=TRUE)
head(my_pca_data)

"
6) Add the PCA rotation axes
What are rotation axes and what is the purpose of drawing these axis?
"
my_pca_rotation <- NULL
my_pca_rotation <- as.data.frame(genes.pca$rotation)


"
7) The PCA plot
"
scale=10  # this helps to scale the rotation axes 
PCAplot <- ggplot(my_pca_data,  aes(x=PC1,y=PC2, colour = factor(Group) )) + geom_point() 
PCAplot <- PCAplot + geom_segment(data=my_pca_rotation, mapping=aes(x=0, y=0, xend=PC1*scale, yend=PC2*scale), size=1.2, colour = "darkgreen")
PCAplot <- PCAplot + geom_text(data=my_pca_rotation, aes(label=row.names(my_pca_rotation), x=PC1*scale, y=PC2*scale), size=5, colour = "darkred") 
PCAplot


##############################################################################################
##                         9. Assignment                                                    ##
##############################################################################################

"
For the logbook the only 2 thing you need to report is:
no intro
no methods

1) A good legend for each plot made by this script which explain the content of the plot properly
It is allowed to modify plots, e.g, colors, sizes of dots etc

2) What can you conclude from the PCA plots


"







