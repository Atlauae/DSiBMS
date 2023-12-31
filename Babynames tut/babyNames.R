"
Tutorial R-Studio, package management and data reduction

Goal of this exercise  

1) How to install and use R-Studio
2) How to load and use libraries, in this tutorial the R library dplyr
3) How to filter large data sets
4) How to group data 
5) How to visualize data
6) How to reduce millions of data points to something that can be interpretated by our human brain


As example we have a data set of more than a million different baby names



Data Science in Biomedicine, 2023

Anne 

"









##############################################################################################
##                        1. set working directory                                          ##
##############################################################################################

" 
The working directory is the folder on your PC or laptop where all input and output is stored.
R offers the possiblility to read data directly from the internet

IMPORTANT TO KNOW: In R the syntax for folder structures is Linux based, meaning that subfolders 
are separated by forward slashes, compared to backslashes in windows. In R (R-studio) backslash 
is an escape character (for more info see link below). 
To be able to use a backslash you need to write it as double back slash.

  example Linux /data/user/tutorial
  example windows C:\\data\\user\\tutorial

==> Escape Characters are common in programming languages, for more info see: 
https://www.w3schools.com/python/gloss_python_escape_characters.asp



==> FIRST STEPs: 
  1)  goto: http://ngs.molgenrug.nl/courses/DataScienceinBiomedicine/VolcanoPCA_tutorial/
  2)  Download data files you need to your laptop folder
  3)  Set to working folder of R-studio to your laptop folder using the 'setwd' command
      Adapt the line below (after the setwd command) to your own working folder on your laptop

==> Put your mouse cursor on the line below and press Ctrl-Enter key combination 
    to execute a line in R-Studio at the cursor position
"


"NOTE: CHANGE THE LINE BELOW TO YOUR LAPTOP FOLDER!!!!!"

setwd("C:\\Users\\Lenovo\\Documents\\github repositories\\DSiBMS")


"To show the content of the current Working Directory execute the line below; again use Ctrl-Enter" 

dir()


"
==> As you just learned; you can type commands in this window and the result will be 
    shown in the window below
"





##############################################################################################
##                        2. Loading data                                                      ##
##############################################################################################

"
The data table we load is a file with Comma Separated Values (aka csv file) with column names
We will use the R function read.table to load the data in the variable babyNames
"

#here the internet link in stored in the variable 'url' to enble easy use later in the script"
url="http://ngs.molgenrug.nl/courses/DataScienceinBiomedicine/babyNames_tutorial/"



# loading the data set

babyNames_Filename <- "babyNames.db"
babyNames_url <- paste(url,babyNames_Filename, sep="")

# Option 1: from file ====> BUT BETTER IS TO USE OPTION 2 <=====
#babyNames<- read.table(babyNames_Filename, sep="," ,header=TRUE)

# Option 2: directly from the internet
babyNames<- read.table(babyNames_url, sep="," ,header=TRUE)


# To check if the data is loaded

# option 1: quickly show the first 6 lines
head(babyNames)

# option 2: use R-studio "Environment" see upper right panel
"NOTE: In the Environment windows you will see babyNames 1352203 obs. of 4 variables"
" Click on babyNames in this window to see the content of babyNames"

# option 3: Open the R-studio environment using the View function
View(babyNames)

# also try: As you see this gives an error because R is CASE SENSITIVE!
view(babyNames)




##############################################################################################
##                        3. Filter data                                                      ##
##############################################################################################


"How many names are in this dataset?"
summary(babyNames)
#Returned length for name is: 1352203, representing the amount of names in the dataset (not filtering out duplicates)

"How many names with Anne?, Try also your own name
NOTE babyNames$Name means; the column Names of data set babyNames "
babyNames[(babyNames$Name == "Anne"),]
anne_count <- 0
for(i in babyNames[(babyNames$Name == "Anne"),3]){
  anne_count <- anne_count + i
}
anne_count
#Final anne_count resulted in a number of 122359 children having a name with Anne

babyNames[(babyNames$Name == "Stefan"),]
stefan_count <- 0
for(i in babyNames[(babyNames$Name == "Stefan"),3]){
  stefan_count <- stefan_count + i
}
stefan_count
#Final stefan_count resulted in a number of 18073 children having a name with Stefan

babyNames[(babyNames$Name == "Femke"),]
femke_count <- 0
for(i in babyNames[(babyNames$Name == "Femke"),3]){
  femke_count <- femke_count + i
}
femke_count
#Final femke_count resulted in a number of 5 children having a name with Femke

"How many names with Anne in 1961?, Try also your own name"
babyNames[(babyNames$Name == "Anne" & babyNames$Year=="1961"),]
name_count <- 0
for(i in babyNames[(babyNames$Name == "Anne" & babyNames$Year=="1961"),3]){
  name_count <- name_count + i
}
name_count
#Final name_count resulted in a number of 5533 children having a name with Anne in the year 1961

"
In the previous example we used == to filter rows. 
Other common relational and logical operators are listed below.
For more see: https://www.statmethods.net/management/operators.html
In most programming languages these operators are similar
 
Operator	Meaning
==        equal to
!=	      not equal to
>	        greater than
>=	      greater than or equal to
<	        less than
<=	      less than or equal to
%in%	    contained in

"


##############################################################################################
##                         4. Loading and/or installing libraries                               ##
##############################################################################################

"
As you learned in DataCamp the library ggplot2 is a powerful method to create nice plots.
More graphical libraries are available in R but in this course we will only use ggplot2
"

# Step 1: loading ggplot2
library(ggplot2)

"
If the ggplot2 is not loading and the result is an error, means that the ggplot2 library
it's NOT available yet on your computer. 
In this case you need to install ggplot once 
"

# Option 1: Via R-Studio menu: Menu => Tools => Install packages


# Option 2: Using the BiocManager as is shown below
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("ggplot2"))   { BiocManager::install("ggplot2", ask=FALSE, dependencies = TRUE)}
if(!require("dplyr"))     { BiocManager::install("dplyr",   ask=FALSE, dependencies = TRUE)}


# If the library is installed you can load ggplot, if not done yet
library(ggplot2)



##############################################################################################
##                   5. Plotting baby name trends over time                                 ##
##############################################################################################


"
Before we make a plot we make a subset
In the example below the data set is 'babyNames' and the column is 'Name' this is written as babyNames$Name
"
Anne <- babyNames[(babyNames$Name == "Anne"),]

# NOTE: you can use the 'View' command to see the data, BUT you need will be directed to another TAB
View(Anne)


"Make a Plot of Year vs Count and using Sex as the Factor for coloring points"
AnnePlot <- ggplot(Anne,  aes(x=Year,y=Count, color=factor(Sex))) + geom_point(pch = 18, cex = 4) 
print(AnnePlot)

"You can change the plot type easily"
AnnePlot <- ggplot(Anne,  aes(x=Year,y=Count, color=factor(Sex))) + geom_line() 
print(AnnePlot)







##############################################################################################
##                   6.The powerfull 'group_by' function of the dplyr library                 ##
##############################################################################################

# Data Transformation with dplyr


library(dplyr)

# If the library is not loaded see above (chapter 4) to install it

"Examples:"
babyNames %>% group_by(Sex) %>% summarise(Count = n())
babyNames %>% group_by(Year) %>% summarise(Count = n())

babyNames %>% group_by(Year, Sex) %>% summarise(Count = n())

"Plotting the number of Boys and Girls per year"
SexPerYear <- babyNames %>% group_by(Year, Sex) %>% summarise(Count = n())
View(SexPerYear)

SexPerYearPlot <- ggplot(SexPerYear,  aes(x=Year,y=Count, color=factor(Sex))) + geom_line() 
print(SexPerYearPlot)


"Get the most popular name of every Year"
MostPopular <- babyNames %>% group_by(Year, Sex) %>%  slice(which.max(Count)) 
MostPopular[(MostPopular$Sex == "Girls"),]

PopularNames <- ggplot(MostPopular,  aes(x=Name)) + geom_bar(width = 0.25) 
PopularNames <- PopularNames + ggtitle("Number of times most popular annual name") 
print(PopularNames)


"
Above are just a few examples, the number of possibilities is endless
In programming people like to make 'Cheat Sheets'


Examples of dplyr in R Cheat Sheets:

https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf
https://nyu-cdsc.github.io/learningr/assets/data-transformation.pdf
or from DataCamp
https://www.datacamp.com/cheat-sheet/data-manipulation-with-dplyr-in-r-cheat-sheet
"

##############################################################################################
##                             ASSIGNMENT babyNames                                         ##
##############################################################################################
library(tidyverse)

"Define a research question for this data set and show the result in a nice colorful and annotated plot"

#Select all baby names from the year 1978 and removing the year column
babyNames_selection <- babyNames[(babyNames$Year == 1987),1:3]

#Seperate the baby names from 1978 in a variable containing boy and girl names
babyNames_selection_boys <- babyNames_selection[(babyNames_selection$Sex == "Boys"),]
babyNames_selection_girls <- babyNames_selection[(babyNames_selection$Sex == "Girls"),]

order(babyNames_selection$Count[(babyNames_selection$Sex == "Boys")], decreasing=TRUE)[1:10]
order(babyNames_selection$Count[(babyNames_selection$Sex == "Girls")], decreasing=TRUE)[1:10]

#Ordering baby names according to count and selecting the top 10 names with the highest count per sex and combining into 1 element
topten_boys <- babyNames_selection_boys[order(babyNames_selection_boys$Count, decreasing=TRUE),][1:10,]
topten_girls <- babyNames_selection_girls[order(babyNames_selection_girls$Count, decreasing=TRUE),][1:10,]
topten_bg <- rbind(topten_boys, topten_girls)
topten_bg$Name <- factor(topten_bg$Name, levels=topten_bg$Name)

#Create a bar plot with the count on the X axis and the names on the Y axis, coloured by sex
ggplot(topten_bg, aes(x = Count, y = Name, color = Sex, fill = Sex)) +
  facet_grid(rows = vars(Sex), scales = "free_y") +
  geom_point(size=6.5) +
  geom_segment(aes(x = 0, xend = Count, yend = Name), linewidth=3.2, position="dodge") +
  geom_text(aes(label = Count), color = "black", size = 3, hjust=-0.5) +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(breaks = seq(0, 70000, by = 10000), expand=expansion(mult = c(0,.12)), position="top") +
  scale_color_manual(values=c("skyblue","hotpink")) +
  labs(title = "Top ten most popular baby names as counted in 1987, for boys and girls") +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color="black"),
        axis.title = element_blank(),
        strip.text = element_blank(),
        legend.position = "right")
