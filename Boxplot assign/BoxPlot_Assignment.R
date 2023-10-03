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



##Making the box plot using the formatted data

#Generate box plot
boxplot <- ggplot(PCMarkers_sel_long, aes(x=variable,y=log(value, 2), fill=factor(diagnosis, labels=c("Control","Benign","PDAC"))))
boxplot <- boxplot + geom_boxplot()

#Add error bar to box plot
boxplot <- boxplot + stat_boxplot(geom="errorbar")

#Add title to box plot
boxplot <- boxplot + ggtitle("Urine Biomarkers in Pancreatic Ductal Adenocarcinoma")

#Add significance indicators to the plot (* p < 0.05, ** p < 0.005, *** p < 0.0005)
significance_indicators <- data.frame(xmin = c(0.75, 1, 0.75, 1.75, 2, 1.75, 2.75, 3, 2.75), 
                                      xmax = c(1, 1.25, 1.25, 2, 2.25, 2.25, 3, 3.25, 3.25), 
                                      y_pos = c(7, 8, 9.5, 13, 14, 15.5, 15.5, 16.5, 18), 
                                      annotations = c("***", "***", "***", "*", "***", "***", "***", "***", "***"))

for(i in 1:nrow(significance_indicators)){
  boxplot <- boxplot + stat_signif(position = "identity",
                                   xmin = significance_indicators[i,1],
                                   xmax = significance_indicators[i,2],
                                   y_position = significance_indicators[i,3],
                                   annotations = significance_indicators[i,4],
                                   tip_length = 0)
}

#Change fill color
boxplot <- boxplot + scale_fill_brewer(palette = "PiYG")

#Add theme to box plot
boxplot <- boxplot + theme_minimal()
boxplot <- boxplot + theme(title = element_text(face = "bold", size = 14))

#Update legend and axis names
boxplot <- boxplot + guides(fill = guide_legend(title = "Diagnosis"))
boxplot <- boxplot + ylab("log2(ng/mg) creatinine")
boxplot <- boxplot + xlab("Biomarker")

print(boxplot)

#Save generated plot as image
png(file = "PDAC_boxplot.png", width = 10, height = 8, units = "in", res=300, pointsize = 12, bg = "white")
print(boxplot)
dev.off()



##Statistical analysis based on plot, using normalized data

for(i in 2:4){
  print(paste("Control vs benign:", 
              t.test(PCMarkers_sel[(PCMarkers_sel$diagnosis == 1),i], PCMarkers_sel[(PCMarkers_sel$diagnosis == 2),i])$p.value, 
              sep=" "))
  print(paste("Control vs fucked:", 
        t.test(PCMarkers_sel[(PCMarkers_sel$diagnosis == 1),i], PCMarkers_sel[(PCMarkers_sel$diagnosis == 3),i])$p.value, 
        sep=" "))
  print(paste("Benign vs fucked:", 
              t.test(PCMarkers_sel[(PCMarkers_sel$diagnosis == 2),i], PCMarkers_sel[(PCMarkers_sel$diagnosis == 3),i])$p.value,
              sep=" "))
}