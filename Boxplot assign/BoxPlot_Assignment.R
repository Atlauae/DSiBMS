#Load libraries
library(ggplot2)
library(ggsignif)
library(reshape2)
library(stats)

#Set working directory
setwd("C:\\Users\\Lenovo\\Documents\\github repositories\\DSiBMS\\Boxplot assign")
WORK_DIR <- getwd()

#Import dataset into R as dataframe
PCMarkers <- read.csv("Debernardi et al 2020 data.csv", header = TRUE, sep = ",")
PCMarkers_doc <- read.csv("Debernardi et al 2020 documentation.csv", header = TRUE, sep = ",")
View(PCMarkers)
View(PCMarkers_doc)

##Most important markers mentioned in original research paper are: Creatinine, YVLE1, REG1B, TFF1
##Creatinine is used as a reference to normalize the data between patients



#Select only the data relevent for our boxplot
PCMarkers_sel <- data.frame(PCMarkers$diagnosis)
colnames(PCMarkers_sel) <- c("diagnosis")
PCMarkers_sel$LYVE1 <- PCMarkers$LYVE1
PCMarkers_sel$REG1B <- PCMarkers$REG1B
PCMarkers_sel$TFF1 <- PCMarkers$TFF1
View(PCMarkers_sel)

#Normalize data based on creatinine levels
PCMarkers_sel$LYVE1 <- PCMarkers_sel$LYVE1 / PCMarkers$creatinine
PCMarkers_sel$REG1B <- PCMarkers_sel$REG1B / PCMarkers$creatinine
PCMarkers_sel$TFF1 <- PCMarkers_sel$TFF1 / PCMarkers$creatinine
View(PCMarkers_sel)

#Reshape dataframe format from wide to long
PCMarkers_sel_long <- melt(PCMarkers_sel, id="diagnosis")
View(PCMarkers_sel_long)

#Generate box plot
boxplot <- ggplot(PCMarkers_sel_long, aes(x=variable,y=log(value, 2), fill=factor(diagnosis, labels=c("Control","Benign","PDAC"))))
boxplot <- boxplot + geom_boxplot()

#Add error bar to box plot
boxplot <- boxplot + stat_boxplot(geom="errorbar")

#Add significance indicators to the plot
boxplot <- boxplot + scale_y_continuous(breaks = seq(-15,15,by=1))
boxplot <- boxplot + stat_signif(position = "identity",
                                 data=data.frame(x=c(0.875, 1.875), xend=c(1.125, 2.125),
                                                 y=c(5.8, 8.5), annotation=c("**", "NS")),
                                 aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))

#Change fill color
boxplot <- boxplot + scale_fill_brewer(palette = "PiYG")

#Add theme to box plot
boxplot <- boxplot + theme_minimal()

#Update legend and axis names
boxplot <- boxplot + guides(fill = guide_legend(title = "Diagnosis"))
boxplot <- boxplot + ylab("log2(ng/mg) creatinine")
boxplot <- boxplot + xlab("Biomarker")

print(boxplot)

##Statistical analysis based on plot, using normalized data

t.test(PCMarkers_sel[(PCMarkers_sel$diagnosis == 1),2], PCMarkers_sel[(PCMarkers_sel$diagnosis == 3),2])$p.value
t.test(PCMarkers_sel[(PCMarkers_sel$diagnosis == 1),2], PCMarkers_sel[(PCMarkers_sel$diagnosis == 3),2])
t.test(PCMarkers_sel[(PCMarkers_sel$diagnosis == 2),2], PCMarkers_sel[(PCMarkers_sel$diagnosis == 3),2])
