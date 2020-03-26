#*******************************************************************************
# SVN Version ID: $Id$
#*******************************************************************************

#Syntax for plotting figures in R that combine a "static" and a stochastic
# model run with all variables against length. 

#!!!IMPORTANT!!! The runs have to be performed with the same amount of 
# E_categories, preferably with an odd number of environments

#Clears workspace to get a fresh start (avoids errors if you are plotting several results in different
# folders after each other)
rm(list=ls())

#Install libraries needed (you only need to do this once)
# might take a few minutes to install, so you have time to get a cup of coffee or tea while waiting :)
#install.packages(c("lattice","corrplot","reshape2","tidyverse","gridExtra", "GGally"), dependencies=TRUE)

#Importing libraries
lapply(c("lattice","corrplot","reshape2","tidyverse", "gridExtra", "GGally"), require, character.only = TRUE)

#Setting working directory (Add your own here)
setwd("./")
getwd() #Getting working directory

#Header of forward figures, useful for identification
# Just write "" if you don't want one 
title.forward<-"3, men E_categories = 4 & AutoCorr = 0.90/1 (4 & 5 U3-18)"


#Information from the model where c(min,max)
length.limits   <- c(10,30)
R_cat.limits    <- c(1,11)
OXF.limits      <- c(0,2000)
GHF.limits      <- c(0,200)
THF.limits      <- c(0,5)

#Output folder:
main.path <- "./Figures-4&5-Combined/"
dir.create(main.path, showWarnings = TRUE) #Creating folder

#Function of importing binary files containing reals
# following http://www.ats.ucla.edu/stat/r/faq/read_binary.htm
import.binary <- function(file_path) {
  bin <- file(file_path, "rb")
  ranking <- readBin(bin,n=1,integer())
  dim_length <- readBin(bin, n = ranking, integer())
  bin.data <- readBin(bin, n = prod(dim_length), numeric())
  bin.data <- array(bin.data,dim=dim_length)
}

#################################################
# Importing data, fixing and combining datasets #
#################################################

#Importing stochastic ind.bin to R as 3D array
ind.stoc.mx <- import.binary("./4 3, men E_categories = 4/Results/ind.bin")

#Importing statind ind.bin file to R as 3D array
ind.stat.mx <- import.binary("./5 3, men E_categories = 4 & AutoCorr = 1/Results/ind.bin")

#Melting 3D arrays into 2D data frames
ind.stoc.df <- melt(ind.stoc.mx, varnames=c("type.of.data", "timestep","individual")) #Stochastic
ind.stat.df <- melt(ind.stat.mx, varnames=c("type.of.data", "timestep","individual")) #Static
rm(ind.stoc.mx, ind.stat.mx) #Removing arrays from memory 

#Adding a new coloumn for model run type to the datasets
ind.stoc.df$run.type <- "Stochastic"
ind.stat.df$run.type <- "Static"

#Combining the two datasets by putting them on top of each other
ind.df <- rbind(ind.stoc.df,ind.stat.df)
rm(ind.stoc.df, ind.stat.df) #Removing old datasets from memory

#Turning the ID of individuals and type.of.data from real into a factor
ind.df$type.of.data <- as.factor(ind.df$type.of.data) 

#List of all names used in ind array
ind.names <- c("length", "Length (cm)",                             #1
              "R_cat", "Reserve category",                          #2      
              "reserves","Reserves (J)",                            #3
              "weight_somatic", "Structural weight (kg)",           #4
              "weight", "Total weight (kg)",                        #5
              "OXF", "Orexin Function (pg ml-1)",                               #6
              "GHF", "Growth Hormone Function (ng ml-1)",                               #7
              "THF", "Thyroid Hormone Function (ng ml-1)",                               #8
              "SMR_std", "SMR standard",                            #9
              "SMR_THF", "SMR adjusted for thyroid hormones",       #10
              "growth", "Growth per week (kg)",                              #11
              "target_intake", "Target intake given OXF in multiples of standard SMR", #12 
              "foraging_cost", "Cost of foraging for food (J timestep-1)",   #13
              "intake", "Intake (J timestep-1)",                     #14
              "SDA", "SDA - Cost of digestion (J timestep -1)",      #15
              "conversion_costs_via_reserves", "Conversion costs from reserves to metabolism (J)", #16
              "O2_used", "Oxygen used",                              #17
              "O2max_THF", "Maximum O2 potential given thyroid hormones",      #18
              "chance_of_max_length", "Chance of reaching max length",         #19
              "M_size", "Size dependent mortality per year",                   #20
              "M_size_O2", "Size dependent free-scope mortality per year",     #21
              "M_size_foraging", "Size dependent foraging mortality per year", #22    
              "survival_next", "Survival",                                     #23
              "survival_step", "Probability of surviving the current timestep",#24
              "reserves_max", "Max reserves",                                  #25
              "M_sizeindependent", "Size independent mortality per year",      #26
              "M_size_02_foraging", "Size dependent free-scope and foraging mortality per year", #27
              "environment", "Environment",                                    #28
              "end_status", "Status at the end",                               #29
              "last_timestep", "Last timestep",                                #30
              "foraging_required", "Foraging intensity required in multiples of standard SMR", #31
              "surplus_before_growth", "Surplus before growth (J timestep-1)", #32
              "reserve_diff","Change in reserves (J timestep-1)", #33
              "conversion_cost_to_growth","Conversion cost to growth rescaled to requirement (J)") #34
              
#Making list with variable names
ind.var.names  <- ind.names[c(TRUE, FALSE)]

#Making list with header names used for plotting
ind.head.names <- ind.names[c(FALSE, TRUE)]

#Making new dataset where each data type has it own column with names based on ind.var.names
ind2.df <- ind.df  %>% mutate(type.of.data = ind.var.names[type.of.data]) %>% spread(key = type.of.data, value = value) 
rm(ind.df) #Removing ind.df from memory 


#Making a new variable based on the environment they started in 
################################################################

#Add a new column to the ind2.df dataset with the environment the
# environment the individual started in called start_environment 
ind2.df <- ind2.df %>% group_by(individual) %>% mutate(start_environment = first(environment))

#Turns this into an integer for better plots
ind2.df$start_environment <- as.factor(ind2.df$start_environment)


#Making 6 groups based on when they reach max length(growth.speed)
##################################################################

#Making a new variable with where all individuals are put into
# growth speed group = NA
ind2.df$growth_speed <- NA

#Making new dataset for timestep 1 
# without the individuals that died or did not finish
ind3.df <- filter(ind2.df, end_status!=1 & timestep==1, run.type=="Stochastic")

#Dividing the number of individuals in the new dataset by 10 and rounding
# this to nearest whole number
nr.x <- round(nrow(ind3.df) / 10, digits=0)

#Calculating upper and lower limit for the individuals of average growth speed
lower.a <- (nrow(ind3.df)/2) - (nr.x/2)               #lower limit
upper.a <- ((nrow(ind3.df)/2) - (nr.x/2) + nr.x) - 1  #upper limit

#Calculating lower limit for the slow group
lower.s <- nrow(ind3.df)- nr.x + 1

#Sort the individuals after when they finished
ind3.df <- ind3.df[order(ind3.df$last_timestep), ]

#Pick out the first 20%, de middle 20% and the last 20% from the dataset
## and giving them new group names
ind.fast <- ind3.df[1:nr.x,]$individual                #List of the fastest individuals...
ind.average <- ind3.df[lower.a:upper.a,]$individual    #... the average individuals ...
ind.slow <- ind3.df[lower.s:nrow(ind3.df),]$individual #... and the slowest individuals

#Dividing invividuals into groups depending on how fast they grew up
ind2.df <- ind2.df %>%
  ungroup %>%
  mutate(growth_speed  = case_when(
  .$individual %in% ind.fast ~ "Fast",
  .$individual %in% ind.average ~ "Average",
  .$individual %in% ind.slow ~ "Slow"
))

#Making a group for the individuals that died
ind2.df[ind2.df$end_status == 1,]$growth_speed <- "Dead"

#Making a group for the individual(s) in the static run
ind2.df[ind2.df$run.type=="Static",]$growth_speed <- "Static"

#Changing the order of the growth speed groups, for better plots
gs.order <- c("Fast","Average", "Slow", "Dead","Static")
ind2.df <- ind2.df %>%
            mutate(growth_speed = factor(growth_speed, levels = gs.order)) %>% 
            arrange(growth_speed) 

#Removing ind3.df from memory
rm(ind3.df)


###########

#Making new variable with rounded length categories 
# to nearest number with one decimal
ind2.df$length_rounded <- round(ind2.df$length, 1)

#Making a new variable with rounded environmental categories 
# to the nearest whole number
ind2.df$environment_rounded <- round(ind2.df$environment)

#Replacing all -1000 with NA for better plots
ind2.df[ind2.df==-1000] <- NA

# Plotting
###########

##Plots with mean of growth speed groups and excluding all the invividuals
## that did not fall into a group per TIMESTEP

#Making new subfolder
path <- paste0(main.path, "timestep-growth-speed/")
dir.create(path, showWarnings = TRUE)

#Counter for loop
count = 1 #1 = length
type = "length"

#Plotting (see above)
#~ for (type in ind.var.names){
  
  #Getting the title from the ind.head.names array
  title.plot <- ind.head.names[count]
  
  X11() #Opening plot in new window
  ind.plot <- ggplot(filter(ind2.df, growth_speed!="NA"), aes_string(x="timestep", y=type, colour = "growth_speed", fill="growth_speed")) +
    stat_summary(geom="ribbon", fun.ymin = min, fun.ymax = max, colour=NA, alpha = 0.1) +
    stat_summary(fun.y=median, geom="line", aes(group = growth_speed)) +
    theme_bw(base_size=16) + 
    scale_x_continuous(limits = c(0,max(ind2.df$timestep))) +
    scale_fill_discrete(guide=FALSE) + 
    ggtitle(title.plot) +
    labs(x="Timestep", y=ind.head.names[count], colour="Growth Speed")
  grid.arrange(ind.plot, top=title.forward, nrow=1, ncol=1)
  
  dev.copy(png, paste0(path, "gs_", ind.head.names[count], ".png"), width=800,height=800)
  dev.off()
  
#~   #Counting number of loops
#~   count = count + 1
#~ }

#Syntax for closing all plots on screen
# and removing them from memory
graphics.off()
rm(ind.plot)




##Plots with mean of growth speed groups and excluding all the invividuals
## that did not fall into a group per LENGTH

#Making new subfolder
path <- paste0(main.path, "length-growth-speed/")
dir.create(path, showWarnings = TRUE)

#Counter for loop
count = 1

#Plotting (see above)
for (type in ind.var.names){
  
  #Getting the title from the ind.head.names array
  title.plot <- ind.head.names[count]
  
  X11() #Opening plot in new window
  ind.plot <- ggplot(filter(ind2.df, growth_speed!="NA"), aes_string(x="length_rounded", y=type, colour = "growth_speed", fill="growth_speed")) +
    stat_summary(geom="ribbon", fun.ymin = min, fun.ymax = max, colour=NA, alpha = 0.1) +
    stat_summary(fun.y=median, geom="line", aes(group = growth_speed)) +
    theme_bw(base_size=16) + 
    scale_x_continuous(limits = length.limits) +
    scale_fill_discrete(guide=FALSE) + 
    scale_y_continuous(labels=function(n){format(n, scientific = FALSE)}) + 
    ggtitle(title.plot) +
    labs(x="Length (cm)", y=ind.head.names[count], colour="Growth Speed")
  grid.arrange(ind.plot, top=title.forward, nrow=1, ncol=1)
  
  dev.copy(png, paste0(path, "gs_", ind.head.names[count], ".png"), width=800,height=800)
  dev.off()
  
  #Counting number of loops
  count = count + 1
}


#Syntax for closing all plots on screen
# and removing them from memory
graphics.off()
rm(ind.plot)

