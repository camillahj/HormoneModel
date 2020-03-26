#*******************************************************************************
# SVN Version ID: $Id$
#*******************************************************************************

#Syntax for plotting figures in R

#Clears workspace to get a fresh start (avoids errors if you are plotting several results in different
# folders after each other)
rm(list=ls())

#Install libraries needed (you only need to do this once)
# might take a few minutes to install, so you have time to get a cup of coffee or tea while waiting :)
#install.packages(c("lattice","corrplot","reshape2","tidyverse","gridExtra", "GGally"), dependencies=TRUE)

#Importing libraries
lapply(c("lattice","corrplot","reshape2","tidyverse", "gridExtra", "GGally"), require, character.only = TRUE)

#Setting working directory: Where you keep your result files (add your own here)
setwd("./Results/")
getwd() #Getting working directory

#Header of forward figures, useful for identification
# Just write "" if you don't want one 
title.forward <- "Ny Default (U16-1)"

#Select how many individuals you want to plot in forward (write "max" if you want all individuals)
# if you don't want max, you can write for example 10 000 without " "
max.ind = "max"

#With out without sizeindependent mortality?
Msizeindependent.check = "with"

#Select which individuals you want individual plots for 
# If you want a random individual write "random" !!!NOT AN OPTION AT THE MOMENT!!!
# For example if you want individual 1, 2, 3 and a random individual write:
# singleind.list <- c(1,2,3,"random")
singleind.list <- c(1,2)

#Information from the model where c(min,max)
length.limits   <- c(10,30)
R_cat.limits    <- c(0,100)
OXF.limits      <- c(0,2500)
GHF.limits      <- c(0,200)
THF.limits      <- c(0,5)

#Output folder:
main.path <- "./Figures/"
dir.create(main.path, showWarnings = TRUE) #Creating folder

#Functions for quantiles used in plots
lower.quant = function(z) {quantile(z,0.0)}
upper.quant = function(z) {quantile(z,1.0)}

#Function of importing binary files containing reals
# following http://www.ats.ucla.edu/stat/r/faq/read_binary.htm
import.binary <- function(file_path) {
  bin <- file(file_path, "rb")
  ranking <- readBin(bin,n=1,integer())
  dim_length <- readBin(bin, n = ranking, integer())
  bin.data <- readBin(bin, n = prod(dim_length), numeric())
  bin.data <- array(bin.data,dim=dim_length)
}

#Function of importing binary files containing integers
# following http://www.ats.ucla.edu/stat/r/faq/read_binary.htm
import.binary.int <- function(file_path) {
  bin <- file(file_path, "rb")
  ranking <- readBin(bin,n=1,integer())
  dim_length <- readBin(bin, n = ranking, integer())
  bin.data <- readBin(bin, n = prod(dim_length), integer())
  bin.data <- array(bin.data,dim=dim_length)
}


#---------------------------------------------------------------------------------------------------------


###############
# Environment #
###############

#Creating subfolder for environment figures
path <- paste0(main.path,"Environment/")
dir.create(path, showWarnings = TRUE)

# ProbOfE
##########

#Importing binary file
probe.mx <- import.binary("./ProbOfE.bin")

#Checking that the data is imported correctly
#summary(probe.mx)

#Plotting
X11() #Makes the plot open in a new window
corrplot(probe.mx, method="number", main="ProbOfE")

#Saving plot
dev.copy(png, paste0(path, "ProbOfE.png"))
dev.off()

# ProbOfAllE
############
#See over

#Importing binary file
proballe.mx <- import.binary("./ProbAllE.bin")

#Checking that the data is imported correctly
#summary(proballe.mx)

#Plotting
#~ X11() #Makes the plot open in a new window
plot(proballe.mx, main="ProbAllE")

#Saving plot
dev.copy(png, paste0(path, "ProbAllE.png"))
dev.off()

# ValueOfE
############

valueofe.mx <- import.binary("./ValueOfE.bin")

#summary(valueofe.mx)

#Plotting
#~ X11() #Makes the plot open in a new window
plot(valueofe.mx,main="ValueOfE")

#Saving plot
dev.copy(png, paste0(path, "ValueOfE.png"))

dev.off()

#Syntax for closing all plots on screen
graphics.off()


#---------------------------------------------------------------------------------------------------------


###############
# fitness.bin #
###############

#Creating new subfolder for figures
path <- paste0(main.path, "Fitness/")
dir.create(path, showWarnings = TRUE)

#Timestep that you want to plot
t <- 1

#Importing fitness.bin to R
fit.mx <- import.binary("./fitness.bin")

#Melting Array
fit.df <- melt(fit.mx, varnames=c("Environment","R_categories","L_categories", "Timestep"))
rm(fit.mx) #Removing array from memory 

#Checking that it is done correctly
#summary(fit.df)

#3D
####
envir <- max(fit.df$Environment)
for (E in 1:envir){
  
  #Plotting
  X11() #Opening plot in new window
  fit.3D.plot <- wireframe(value~L_categories*R_categories, data = filter(fit.df, Timestep==t, Environment==E),
            main = "Fitness",
            xlab = "Length (cm)", ylab = "Reserve categories", zlab="",
            scales = list(arrows=FALSE, cex=.5, tick.number = 10, z = list(arrows=T)),
            drape = TRUE,
            colorkey = TRUE)

  grid.arrange(fit.3D.plot, top = paste0("Environment: ", E), nrow=1, ncol=1)

  #Creating a number for each plot file with leading zeros
  if (E < 10) {nr = paste0("000",E)}
  if (E < 100 && E >= 10) {nr = paste0("00",E)}
  if (E >= 100) {nr = paste0("0", E)}
  
  #Saving plot
  dev.copy(png, paste0(path, "Fitness ", nr,".png"))
  dev.off()
}

#Syntax for closing all plots on screen
graphics.off()
rm(fit.3D.plot)
rm(fit.mx, fit.df)


#---------------------------------------------------------------------------------------------------------


################
# strategy.bin #
################

#Making new subfolder for figures
path <- path <- paste0(main.path, "Strategy/")
dir.create(path, showWarnings = TRUE)

#Timestep that you want to plot
t <- 1

#Importing fitness.bin to R
strategy.mx <- import.binary.int("./strategy.bin")

#Melting array
strategy.df <- melt(strategy.mx, varnames=c("Hormonetype","Environment","R_categories","L_categories", "Timestep"))
rm(strategy.mx) #Removing array from memory 

#Checking that it is done correctly
#summary(strategy.df)

#Making subsets for the different hormones, for timestep t
strategy.orexin.df <- filter(strategy.df, Hormonetype==1, Timestep==t)
strategy.growth.df <- filter(strategy.df, Hormonetype==2, Timestep==t)
strategy.thyroid.df <-filter(strategy.df, Hormonetype==3, Timestep==t)

#3D plots
#  from https://www.r-bloggers.com/creating-surface-plots/
##########################################################

#Saving the number of environments
envir <- max(strategy.df$Environment)

#Makes a plot for every environment and saves them, so they can be used for making gifs
for (E in 1:envir){
  #Orexin
  orexin.wire <- wireframe(value~L_categories*R_categories, data = filter(strategy.orexin.df,Environment==E),
            xlab = "Length category", ylab = "Reserve category", zlab="Orexin",
            main = "Optimal Orexin Category",
            drape = TRUE,
            colorkey = TRUE
  )

  #Growth hormone
  growth.wire <- wireframe(value ~ L_categories*R_categories, data = filter(strategy.growth.df,Environment==E),
            xlab = "Length category", ylab = "Reserve category", zlab="IGF-1",
            main = "Optimal IGF Category",
            drape = TRUE,
            colorkey = TRUE
  )

  #Thyroid
  thyroid.wire <- wireframe(value ~ L_categories*R_categories, data = filter(strategy.thyroid.df,Environment==E),
            xlab = "Length category", ylab = "Reserve category", zlab="T3",
            main = "Optimal T3 Category",
            drape = TRUE,
            colorkey = TRUE
  )
  
  #Plotting all in the same window
  X11(width=12,height=12) #Opens a new window and of with width=12" and height = 12"
  grid.arrange(orexin.wire, growth.wire,thyroid.wire, top = paste0("Environment: ", E), nrow=2, ncol=2)

  #Creating a number for each plot file with leading zeros
  if (E < 10) {nr = paste0("000",E)}
  if (E < 100 && E >= 10) {nr = paste0("00",E)}
  if (E >= 100) {nr = paste0("0", E)}

  #Saving plot
  dev.copy(png, paste0(path, "Optimal_Hormones_3D ", nr,".png"), width=800,height=800)
  dev.off()
}

#Syntax for closing all plots on screen
graphics.off()

rm(strategy.orexin.df, strategy.growth.df, strategy.thyroid, orexin.wire, growth.wire, thyroid.wire)
rm(strategy.df)

#---------------------------------------------------------------------------------------------------------


###########
# ind.bin #
###########

#Importing ind.bin to R
ind.mx <- import.binary("./ind.bin")

#Melting 3D array into 2D data frame
ind.df <- melt(ind.mx, varnames=c("type.of.data", "timestep","individual"))
rm(ind.mx) #Removing array from memory 

#Making new folder for forward figures
path.forward <- paste0(main.path, "Forward/")
dir.create(path.forward, showWarnings = TRUE)

#Checking that it is done correctly
#head(ind.df)
#str(ind.df)
#summary(ind.df)

#Turning the ID of individuals and type.of.data from real into a factor
ind.df$type.of.data <- as.factor(ind.df$type.of.data) 

#Selecting how many individuals you want to plot if you have chosen a limit
if(max.ind != "max"){ind.df <- filter(ind.df, individual >= 1 & individual <= max.ind)}

#Variable with the length for halfway to maxlength
halvveis <- round(length.limits[1]+((length.limits[2]-length.limits[1])/2))

#List of all names used in ind array
ind.names <- c("length", "Length (cm)",                             #1
              "R_cat", "Reserve fullness",                          #2      
              "reserves","Reserves (J)",                            #3
              "weight_somatic", "Structural weight (kg)",           #4
              "weight", "Total weight (kg)",                        #5
              "OXF", "OXF (pg ml-1)",                               #6
              "GHF", "GHF (ng ml-1)",                               #7
              "THF", "THF (ng ml-1)",                               #8
              "SMR_std", "SMR standard",                            #9
              "SMR_THF", "SMR adjusted for thyroid hormones",       #10
              "growth", "Growth (kg)",                              #11
              "target_intake", "Target intake given OXF in multiples of standard SMR", #12 
              "foraging_cost", "Cost of foraging for food (J timestep-1)",   #13
              "intake", "Intake (J timestep-1)",                     #14
              "SDA", "SDA - Cost of digestion (J timestep -1)",      #15
              "conversion_costs_via_reserves", "Conversion costs from intake to reserves, and from reserves to metabolism (J)", #16
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
              "M_size_O2_foraging", "Size dependent free-scope and foraging mortality per year", #27
              "environment", "Environment",                                    #28
              "end_status", "Status at the end",                               #29
              "last_timestep", "Last timestep",                                #30
              "foraging_required", "Foraging intensity required in multiples of standard SMR", #31
              "surplus_before_growth", "Surplus before growth (J timestep-1)", #32
              "reserve_diff","Change in reserves (J timestep-1)", #33
              "conversion_cost_to_growth","Conversion cost to growth from reserves (J timestep -1)") #34
              
#Making list with variable names
ind.var.names  <- ind.names[c(TRUE, FALSE)]

#Making list with header names used for plotting
ind.head.names <- ind.names[c(FALSE, TRUE)]

#Making new dataset where each data type has it own column with names based on ind.var.names
ind2.df <- ind.df  %>%
            mutate(type.of.data = ind.var.names[type.of.data]) %>%
            spread(key = type.of.data, value = value) %>%
            filter(length!=-1000)
rm(ind.df) #Removing ind.df from memory 

#Checking that it is done correctly
#head(ind.df)
#str(ind.df)
#summary(ind.df)

#Making a new variable based on the environment they started in 
################################################################

#Add a new column to the ind2.df dataset with the environment the
# environment the individual started in called start_environment 
ind2.df <- ind2.df %>% group_by(individual) %>% mutate(start_environment = first(environment))

#Turns this into an integer for better plots
ind2.df$start_environment <- as.integer(ind2.df$start_environment)


#Making 6 groups based on when they reach max length(growth.speed)
##################################################################

#Making a new variable with where all individuals are put into
# growth speed group = NA
ind2.df$growth_speed <- NA

#Making new dataset for timestep 1 
# without the individuals that died or did not finish
ind3.df <- filter(ind2.df, end_status!=1 & timestep==1)

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

#Changing the order of the growth speed groups, for better plots
gs.order <- c("Fast","Average", "Slow", "Dead")
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

### Control plot to check the individuals for which environments
### they started in, to compare it to the curve of ProbAllE

X11() #Makes the plot open in a new window
test1.plot <-ggplot(filter(ind2.df, timestep==1), aes(x=start_environment)) + 
  geom_histogram(binwidth=1., colour="black", fill="white") +
  scale_x_continuous(limits=c(0,(max((ind2.df$start_environment))+1))) +
  theme_bw(base_size=16) + 
  ggtitle("Start Environment Test") +
  labs(x="Start Environment", y="Number of Individuals")
grid.arrange(test1.plot, top=title.forward, nrow=1, ncol=1)

#Saving plot
dev.copy(png, paste0(path.forward, "test1_environment.png"), width=800,height=800)
dev.off()


###Control plot to check how the individuals change their environments
### using the two first individuals

#Plotting
#~ x11() #Opening a new plotting window
test2.plot <- ggplot(filter(ind2.df, individual >= 1 & individual <= 2),
                    aes(x=timestep, y=environment, colour=as.factor(individual))) +
    geom_line() +
    geom_hline(yintercept = 1) + geom_hline(yintercept = max(ind2.df$environment, na.rm=TRUE)) +
    scale_y_continuous(limits = c(1,max(ind2.df$environment, na.rm=TRUE))) +
    theme_bw(base_size=16) + 
    ggtitle("Environment Test") +
    labs(x="Timestep", y="Environment", colour="Individual")
grid.arrange(test2.plot, top=title.forward, nrow=1, ncol=1)
  
dev.copy(png, paste0(path.forward, "test2_environment.png"), width=800,height=800)
dev.off()
  
#Syntax for closing all plots on screen
#~ graphics.off()



### Control plot to check the intake in different environments

#Making a new dataset with only individual ID, environment, intake and growth speed group
test3.df <- select(ind2.df, individual, environment, intake, growth_speed)

#Rounding environment to nearest whole number
test3.df$environment <- round(test3.df$environment)

#Calculating mean intake for the different environments
test3.mean <- aggregate(intake ~ environment + growth_speed, test3.df, mean)

#Plotting
#~ x11() #Opening a new plotting window
test3.plot <- ggplot(test3.mean, aes(x=environment, y=intake, colour=growth_speed)) +
    geom_line() +
    theme_bw(base_size=16) + 
    ggtitle("Intake Environment Test") +
    labs(x="Environment", y="Intake", colour="Growth Speed")
grid.arrange(test3.plot, top=title.forward, nrow=1, ncol=1)

dev.copy(png, paste0(path.forward, "test3_intake vs environment", ".png"), width=800,height=800)
dev.off()



### Control plot to check max intake and intake for the fast speed group

#Max intake vs intake
#~ X11() #Opening plot in new window
test4.plot <- ggplot(filter(ind2.df, growth_speed=="Fast"), aes(x=timestep, y=intake, colour = growth_speed, fill=growth_speed)) +
  geom_point() +
  stat_summary(geom="ribbon", fun.ymin = lower.quant, fun.ymax = upper.quant, colour=NA, alpha = 0.1) +
  stat_summary(fun.y=mean, geom="line", aes(group = growth_speed)) +
  stat_summary(fun.y=mean, geom="line", aes(y=intake_max, group = growth_speed, colour="max intake")) +
  theme_bw(base_size=16) + 
  scale_linetype_manual(values=c("twodash", "dotted")) +
  scale_x_continuous(limits = c(0,max(ind2.df$timestep))) +
  scale_fill_discrete(guide=FALSE) + 
  ggtitle("Intake vs Max Intake Test") +
  labs(x="Timestep", y="Intake", colour="Growth Speed")
grid.arrange(test4.plot, top=title.forward, nrow=1, ncol=1)
  
dev.copy(png, paste0(path.forward, "test4_intake vs max intake", ".png"), width=800,height=800)
dev.off()




##Plots with mean of growth speed groups and excluding all the invividuals
## that did not fall into a group per TIMESTEP

#Making new subfolder
path <- paste0(path.forward, "timestep-growth-speed/")
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
    stat_summary(geom="ribbon", fun.ymin = lower.quant, fun.ymax = upper.quant, colour=NA, alpha = 0.1) +
    stat_summary(fun.y=mean, geom="line", aes(group = growth_speed)) +
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
graphics.off()




##Plots with mean of growth speed groups and excluding all the invividuals
## that did not fall into a group per LENGTH

#Making new subfolder
path <- paste0(path.forward, "length-growth-speed/")
dir.create(path, showWarnings = TRUE)

#Counter for loop
count = 1

#Plotting (see above)
for (type in ind.var.names){
  
  #Getting the title from the ind.head.names array
  title.plot <- ind.head.names[count]
  
  X11() #Opening plot in new window
  ind.plot <- ggplot(filter(ind2.df, growth_speed!="NA"), aes_string(x="length_rounded", y=type, colour = "growth_speed", fill="growth_speed")) +
    stat_summary(geom="ribbon", fun.ymin = lower.quant, fun.ymax = upper.quant, colour=NA, alpha = 0.1) +
    stat_summary(fun.y=mean, geom="line", aes(group = growth_speed)) +
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
graphics.off()



#Removing plots from memory
rm(list=ls(pattern=".plot"))



#---------------------------------------------------------------------------------------------------------



#######################
# Miscellaneous plots #
#######################

##Making new folder for this plot type
path <- paste0(path.forward, "misc-plots/")
dir.create(path, showWarnings = TRUE)


##Plots for reserve categories in different environments
environment.R_cat.plot <- ggplot(filter(ind2.df, length_rounded==halvveis), aes(x=R_cat, colour=as.factor(environment_rounded))) +
                          geom_freqpoly(bins = max(ind2.df$R_cat, na.rm=TRUE)) +
                          facet_wrap(~environment_rounded) +
                          guides(colour=FALSE) +
                          theme_bw(base_size=16) + 
                          scale_x_continuous(limits = c(0,max(ind2.df$R_cat, na.rm=TRUE))) +
                          ggtitle(paste("Reserve fullness in different environments at",halvveis,"cm")) +
                          labs(x="Reserve fullness (%)", y="# individuals")
grid.arrange(environment.R_cat.plot, top=title.forward)

#Saving and closing plot
dev.copy(png, paste0(path, "environment_R_cat_halvveis.png"), width=800,height=800)
dev.off()


##Plots for hormones given reserve category and length

#Function for plotting with reserve categories on the y axis, length on the x axis
# and hormones as a colour scale
ggplot2_l.r.h <- function(hormone, label, title, hormone.limits) {
  l.r.h.plot <- ggplot(ind2.df, aes_string(y="R_cat", x="length", z=hormone)) + 
                stat_summary_2d(fun=mean) +
                scale_fill_gradient(low="yellow", high="red", limits = hormone.limits) +
                scale_x_continuous(limits = length.limits) +
                scale_y_continuous(limits = R_cat.limits) +
                theme_bw(base_size=16) + 
                ggtitle(title) + labs(y="Reserve fullness (%)", x="Length (cm)", fill=label)
 }

#Plotting
#OXF
for.oxf.plot <- ggplot2_l.r.h("OXF", "OXF", "Hormone levels", OXF.limits)

#GHF
for.ghf.plot <- ggplot2_l.r.h("GHF", "GHF", "", GHF.limits)

#THF
for.thf.plot <- ggplot2_l.r.h("THF", "THF", "", THF.limits)

#Plotting to screen
grid.arrange(for.oxf.plot, for.ghf.plot, for.thf.plot, top=title.forward, nrow=2, ncol=2)

#Saving and closing plot
dev.copy(png, paste0(path, "hormonelevels_length_R_cat.png"), width=800,height=800)
dev.off()
#~ graphics.off()

#Removing plots from memory
rm(for.oxf.plot, for.ghf.plot, for.thf.plot)


##Plots for hormones given length, reserve category and environment

#Function for plotting with hormone levels on the y axis, length on the x axis
# and reservecategory as a colour scale and environment as facet
ggplot2_l.r.h.e <- function(hormone, label, title, hormone.limits) {
  l.r.h.plot <- ggplot(filter(ind2.df, environment_rounded!="NA"), aes_string(x="length", y=hormone, colour="R_cat")) + 
                geom_point() +
                facet_wrap(~environment_rounded) +
                scale_x_continuous(limits = length.limits) +
                scale_y_continuous(limits = hormone.limits) +
                scale_colour_continuous(limits = c(0,max(ind2.df$R_cat, na.rm=TRUE))) +
                theme_bw(base_size=16) + 
                ggtitle(title) + labs(x="Length (cm)", y=label, colour="Reserve \n category")
 }

#Plotting and saving
#OXF
for.en.oxf.plot <- ggplot2_l.r.h.e("OXF", "OXF", "OXF in different environments", OXF.limits)
grid.arrange(for.en.oxf.plot, top=title.forward)
dev.copy(png, paste0(path, "hormonelevels_OXF~length+environment+R_cat.png"), width=800,height=800)
dev.off()

#GHF
for.en.ghf.plot <- ggplot2_l.r.h.e("GHF", "GHF", "GHF in different environments", GHF.limits)
grid.arrange(for.en.ghf.plot, top=title.forward)
dev.copy(png, paste0(path, "hormonelevels_GHF~length+environment+R_cat.png"), width=800,height=800)
dev.off()

#THF
for.en.thf.plot <- ggplot2_l.r.h.e("THF", "THF", "THF in different environments", THF.limits)
grid.arrange(for.en.thf.plot, top=title.forward)
dev.copy(png, paste0(path, "hormonelevels_THF~length+environment+R_cat.png"), width=800,height=800)
dev.off()

#Closing plots
#~ graphics.off()


#~ #Function for plotting with reserve categories on the y axis, length on the x axis
#~ # and hormones as a colour scale and environment as facet
#~ ggplot2_l.r.h.e2 <- function(hormone, label, title, hormone.limits) {
#~   l.r.h.plot2 <- ggplot(filter(ind2.df, environment_rounded!="NA"), aes_string(x="length", y="R_cat", colour=hormone)) + 
#~                 geom_point() +
#~                 facet_wrap(~environment_rounded) +
#~                 scale_x_continuous(limits = length.limits) +
#~                 #scale_y_continuous(limits = R_cat.limits) +
#~                 scale_colour_continuous(limits = hormone.limits) +
#~                 theme_bw(base_size=16) + 
#~                 ggtitle(title) + labs(x="Length (cm)", y="Reserve fullness (%)", z=label)
#~  }

#~ #Plotting and saving
#~ #OXF
#~ for.en.oxf.plot2 <- ggplot2_l.r.h.e2("OXF", "OXF", "OXF and Reserve Categories in different environments", OXF.limits)
#~ grid.arrange(for.en.oxf.plot2, top=title.forward)
#~ dev.copy(png, paste0(path, "R_cat_OXF~length+environment+hormonelevels.png"), width=800,height=800)
#~ dev.off()

#~ #GHF
#~ for.en.ghf.plot2 <- ggplot2_l.r.h.e2("GHF", "GHF", "GHF and Reserve Categories in different environments", GHF.limits)
#~ grid.arrange(for.en.ghf.plot2, top=title.forward)
#~ dev.copy(png, paste0(path, "R_cat_GHF~length+environment+hormonelevels.png"), width=800,height=800)
#~ dev.off()

#~ #THF
#~ for.en.thf.plot2 <- ggplot2_l.r.h.e2("THF", "THF", "THF and Reserve Categories in different environments", THF.limits)
#~ grid.arrange(for.en.thf.plot2, top=title.forward)
#~ dev.copy(png, paste0(path, "R_cat_THF~length+environment+hormonelevels.png"), width=800,height=800)
#~ dev.off()

#Closing plots
#~ graphics.off()



##Histograms for the timestep the individuals died

#Making dataset with only dead individuals
ind2.dead.df <- filter(ind2.df, growth_speed=="Dead") 

#Histogram with the timestep when the individuals died
#~ x11()
timestep.dead.plot <- ggplot(filter(ind2.dead.df, timestep==1), aes(x = last_timestep)) +
                      geom_histogram(bins=max(ind2.df$timestep)) +
                      scale_x_continuous(limits = c(1, max(ind2.df$timestep))) +
                      theme_bw(base_size=16)+
                      labs(title="Timestep individual died", x= "Timestep", y="Count")
grid.arrange(timestep.dead.plot, top=title.forward, nrow=1, ncol=1)

#Saving plot
dev.copy(png, paste0(path, "timestep_of_death.png"), width=800,height=800)
dev.off()

#Closing and resetting plotting window
#~ graphics.off()


#Histogram with length at death
#New dataset with max length for every individual
death.length.df <- aggregate(length ~ individual, data = ind2.dead.df, max)

#Histogram
#~ x11()
length.dead.plot <- ggplot(death.length.df, aes(x = length)) +
                      geom_histogram(bins=(length.limits[2]-length.limits[1])*10) +
                      scale_x_continuous(limits = length.limits) +
                      theme_bw(base_size=16)+
                      labs(title="Length when individual died", x= "Length (cm)", y="Count")
grid.arrange(length.dead.plot, top=title.forward, nrow=1, ncol=1)

#Saving plot
dev.copy(png, paste0(path, "length_of_death.png"), width=800,height=800)
dev.off()

#Closing and resetting plotting window
#~ graphics.off()


#Simple histogram with each individual in each length category
hist(ind2.df$length_rounded, xlim=length.limits, breaks=301)

#Saving and closing plot
dev.copy(png, paste0(path, "length_groups.png"), width=800,height=800)
dev.off()
#~ graphics.off()


##Plot of hormone levels given environment

#Function for making plots
ggplot2_h.e.gs <- function(hormone, hormone.limits, length.to.plot) {
  oxf.e.gs.plot <- ggplot(filter(ind2.df, length_rounded==length.to.plot), aes_string(y=hormone, x="environment", colour="R_cat")) +
                    geom_point() +
                    theme_bw(base_size=16) + 
                    scale_x_continuous(limits = c(1, max(ind2.df$environment, na.rm=TRUE))) +
                    scale_y_continuous(limits = hormone.limits) +
                    scale_colour_continuous(limits = c(0,max(ind2.df$R_cat, na.rm=TRUE))) +
                    ggtitle(paste("Length =",length.to.plot)) + labs(y=hormone, x="Environment", shape="Reserve \n category")
}

#Making plots for OXF using function
OXF_envir1.plot <- ggplot2_h.e.gs("OXF", OXF.limits, length.limits[1])
OXF_envir2.plot <- ggplot2_h.e.gs("OXF", OXF.limits, round(length.limits[1]+((length.limits[2]-length.limits[1])/3)))
OXF_envir3.plot <- ggplot2_h.e.gs("OXF", OXF.limits, round(length.limits[1]+((length.limits[2]-length.limits[1])/3)*2))
OXF_envir4.plot <- ggplot2_h.e.gs("OXF", OXF.limits, length.limits[2])

#Plotting all to screen
grid.arrange(OXF_envir1.plot, OXF_envir2.plot, OXF_envir3.plot, OXF_envir4.plot, top=title.forward, nrow=2, ncol=2)

#Saving and closing plot
dev.copy(png, paste0(path, "OXFlevels_environment_length.png"), width=800,height=800)
dev.off()
#~ graphics.off()

#Making plots for GHF using function
GHF_envir1.plot <- ggplot2_h.e.gs("GHF", GHF.limits, length.limits[1])
GHF_envir2.plot <- ggplot2_h.e.gs("GHF", GHF.limits, round(length.limits[1]+((length.limits[2]-length.limits[1])/3)))
GHF_envir3.plot <- ggplot2_h.e.gs("GHF", GHF.limits, round(length.limits[1]+((length.limits[2]-length.limits[1])/3)*2))
GHF_envir4.plot <- ggplot2_h.e.gs("GHF", GHF.limits, length.limits[2])

#Plotting to screen
grid.arrange(GHF_envir1.plot, GHF_envir2.plot, GHF_envir3.plot, GHF_envir4.plot, top=title.forward, nrow=2, ncol=2)

#Saving and closing plot
dev.copy(png, paste0(path, "GHFlevels_environment_length.png"), width=800,height=800)
dev.off()
#~ graphics.off()

#Making plots for GHF using function
THF_envir1.plot <- ggplot2_h.e.gs("THF", THF.limits, length.limits[1])
THF_envir2.plot <- ggplot2_h.e.gs("THF", THF.limits, round(length.limits[1]+((length.limits[2]-length.limits[1])/3)))
THF_envir3.plot <- ggplot2_h.e.gs("THF", THF.limits, round(length.limits[1]+((length.limits[2]-length.limits[1])/3)*2))
THF_envir4.plot <- ggplot2_h.e.gs("THF", THF.limits, length.limits[2])

#Plotting to screen
grid.arrange(THF_envir1.plot, THF_envir2.plot, THF_envir3.plot, THF_envir4.plot, top=title.forward, nrow=2, ncol=2)

#Saving and closing plot
dev.copy(png, paste0(path, "THFlevels_environment_length.png"), width=800,height=800)
dev.off()
#~ graphics.off()

#Plotting and saving just with hormone values when they are half way to max length

#OXF given environment when they are half way to max length, plotting to screen and saving
OXF_envir_halfway.plot <- ggplot2_h.e.gs("OXF", OXF.limits, halvveis)
grid.arrange(OXF_envir_halfway.plot, top=title.forward, nrow=1, ncol=1)
dev.copy(png, paste0(path, "OXFlevels_environment_halfway_to_maxlength.png"), width=800,height=800)
dev.off()

#GHF given environment when they are half way to max length, plotting to screen and saving
GHF_envir_halfway.plot <- ggplot2_h.e.gs("GHF", GHF.limits, halvveis)
grid.arrange(GHF_envir_halfway.plot, top=title.forward, nrow=1, ncol=1)
dev.copy(png, paste0(path, "GHFlevels_environment_halfway_to_maxlength.png"), width=800,height=800)
dev.off()

#THF given environment when they are half way to max length, plotting to screen and saving
THF_envir_halfway.plot <- ggplot2_h.e.gs("THF", THF.limits, halvveis)
grid.arrange(THF_envir_halfway.plot, top=title.forward, nrow=1, ncol=1)
dev.copy(png, paste0(path, "THFlevels_environment_halfway_to_maxlength.png"), width=800,height=800)
dev.off()

#Closing plots
#~ graphics.off()




#Tables
#######

#~ x11() #Opening new plotting window

#For end_status groups
end_status.table.df <- matrix(c(nrow(filter(ind2.df, end_status==1 & timestep==1)),
              nrow(filter(ind2.df, end_status==2 & timestep==1)),
              nrow(filter(ind2.df, end_status==3 & timestep==1)),
              nrow(filter(ind2.df, timestep==1))),
              ncol=4, byrow=FALSE)
colnames(end_status.table.df) <- c("Dead","To Slow","Reached max length", "Total")
end_status.table <- tableGrob(end_status.table.df)


#For speed groups
growth_speed.table.df <- matrix(c(nrow(filter(ind2.df, growth_speed == "Fast" & timestep==1)),
              nrow(filter(ind2.df, growth_speed == "Average" & timestep==1)),
              nrow(filter(ind2.df, growth_speed == "Slow" & timestep==1)),
              nrow(filter(ind2.df, growth_speed != "NA" & growth_speed != "Dead" & timestep==1))),
              ncol=4, byrow=FALSE)
colnames(growth_speed.table.df) <- c("Fast","Average","Slow", "Total")
growth_speed.table <- tableGrob(growth_speed.table.df)

#Printing them in plotting window
grid.arrange(end_status.table, growth_speed.table, top=title.forward)

#Saving tables
dev.copy(png, paste0(path.forward, "summary_tables.png"), width=800,height=800)
dev.off()

#Syntax for closing all plots on screen
#~ graphics.off()



#Removing plots from memory
rm(list=ls(pattern=".plot"))


#---------------------------------------------------------------------------------------------------------


###################
# Mortality plots #
###################

#Making new folder for figures
path <- paste0(path.forward, "mortality/")
dir.create(path, showWarnings = TRUE)

#Making new dataset with mortality data stacked on top of each other
mort.ind.df <- gather(data=ind2.df, key=mortality_type, value=mortality_data, M_size:M_size_O2_foraging)

#Changing the order of the mortality factors
mort.order <- c("M_size_O2_foraging","M_size_O2", "M_size_foraging", "M_size","M_sizeindependent")
mort.ind.df <- mort.ind.df %>%
               mutate(mortality_type = factor(mortality_type, levels = mort.order)) %>% 
               arrange(mortality_type) 

#Function for plotting mortality plots
ggplot2_mort <- function(df, m.title) {
  m.plot <- ggplot(df, aes(x=length_rounded, y=mortality_data, fill=mortality_type)) +
                    geom_area(colour="black", size=.2, alpha=.4) +
                    theme_bw(base_size=16) + 
                    coord_cartesian(xlim = length.limits, ylim = c(0, 1.5)) +
                    labs(x="Length", y="Mortality (per year)", fill="Mortality type") +
                    ggtitle(m.title)
}


##Mortality plot for the whole population

#~ #Rounding lengths to closest number with two decimal places, and making it into a new variable
#~ mort.ind.df$length_rounded <- round(mort.ind.df$length, 1)

#Calculating the means of the different mortality types based on length and type for
# the whole population
all.mort.mean <- aggregate(mortality_data ~ length_rounded + mortality_type, mort.ind.df, mean)

#Plot with the mean mortality values of the whole population
all_mort_length.plot <- ggplot2_mort(all.mort.mean, "Mean mortality values for the population")
grid.arrange(all_mort_length.plot, top=title.forward, nrow=1, ncol=1)

#Saving and closing plot
dev.copy(png, paste0(path,"all_mortality_length.png"), width=800,height=800)
dev.off()
#~ graphics.off()


##Calculating means and plotting for the 3 speed groups

#Slow
slow.mort.mean <- aggregate(mortality_data ~ length_rounded + mortality_type, filter(mort.ind.df,
    growth_speed=="Slow"), mean)
slow_mort_length.plot <- ggplot2_mort(slow.mort.mean, "Slow")

#Average
average.mort.mean <- aggregate(mortality_data ~ length_rounded + mortality_type, filter(mort.ind.df,
    growth_speed=="Average"), mean)
average_mort_length.plot <- ggplot2_mort(average.mort.mean, "Average")

#Fast
fast.mort.mean <- aggregate(mortality_data ~ length_rounded + mortality_type, filter(mort.ind.df,
    growth_speed=="Fast"), mean)
fast_mort_length.plot <- ggplot2_mort(fast.mort.mean, "Fast")

#Plot window with all 3 plots
grid.arrange(slow_mort_length.plot, average_mort_length.plot, fast_mort_length.plot,
    top=title.forward, nrow=2, ncol=2)

#Saving and closing plot
dev.copy(png, paste0(path,"growth_speed_mortality_length.png"), width=800,height=800)
dev.off()
#~ graphics.off()


## Plots with hormones against mortality

#Function for plotting mortality plots with hormones on xaxis
ggplot2_h_mort <- function(hormone, h.limits) {
  
  #Making plot
  h.m.plot <- ggplot(filter(mort.ind.df, mortality_type!="M_sizeindependent",individual<=1000), aes_string(x=hormone, y="mortality_data", colour="length")) +
                    geom_point() +
                    theme_bw(base_size=16) + 
                    scale_x_continuous(limits = h.limits) +
                    scale_y_continuous(limits = c(0, 1.5)) +
                    facet_wrap(~mortality_type)+
                    labs(x=hormone, y="Mortality per year", colour="Length") +
                    ggtitle(paste("Mortality ~",hormone, "(Individual 1-1000)"))
}

#Plotting for OXF and saving
OXF_h.m.plot <- ggplot2_h_mort("OXF", OXF.limits) #Making plot
grid.arrange(OXF_h.m.plot, top=title.forward) #Plotting to screen 
dev.copy(png, paste0(path,"mortality_OXF.png"), width=800,height=800) #Saving
dev.off()

#Plotting for GHF and saving
GHF_h.m.plot <- ggplot2_h_mort("GHF", GHF.limits)
grid.arrange(GHF_h.m.plot, top=title.forward) 
dev.copy(png, paste0(path,"mortality_GHF.png"), width=800,height=800)
dev.off()

#Plotting for THF and saving
THF_h.m.plot <- ggplot2_h_mort("THF", THF.limits)
grid.arrange(THF_h.m.plot, top=title.forward) 
dev.copy(png, paste0(path,"mortality_THF.png"), width=800,height=800)
dev.off()

#~ #Plotting for R_cat and saving
#~ R_cat.m.plot <- ggplot2_h_mort("R_cat", R_cat.limits)
#~ grid.arrange(R_cat.m.plot, top=title.forward) 
#~ dev.copy(png, paste0(path,"mortality_R_cat.png"), width=800,height=800)
#~ dev.off()

#Plotting for environment and saving
E_h.m.plot <- ggplot2_h_mort("environment", c(min(ind2.df$environment_rounded, na.rm=TRUE),max(ind2.df$environment_rounded, na.rm=TRUE)))
grid.arrange(E_h.m.plot, top=title.forward)   
  dev.copy(png, paste0(path,"mortality_E.png"), width=800,height=800)
dev.off()

#Function for simple plot with x, y and colour for a chosen length
plot.simple <- function(dataset.df,xaxis, xlabel, yaxis, ylabel, colgroup,plot.length){
  easy.plot <- ggplot(filter(dataset.df, length_rounded==plot.length), aes_string(x=xaxis, y=yaxis, colour=colgroup)) +
                geom_point() +
                theme_bw(base_size=16) +
                labs(x=xlabel,y=ylabel) +
                ggtitle(paste("Length =", plot.length))
}


#Making new dataset excluding M_sizeindependent
mort.ind2.df <- filter(mort.ind.df, mortality_type!="M_sizeindependent")

#Making plot with mortality given THF, colour R_cat
THF_mort_Rcat.plot <- plot.simple(mort.ind2.df, "THF","THF (ng ml -1)","mortality_data","Mortality per year","R_cat",halvveis) +
              facet_wrap(~mortality_type)
grid.arrange(THF_mort_Rcat.plot, top=title.forward) 
dev.copy(png, paste0(path,"THF_mortality_halfway_to_maxlength.png"), width=800,height=800) #Saving
dev.off()

#Making plot with mortality given GHF, colour R_cat
GHF_mort_Rcat.plot <- plot.simple(mort.ind2.df, "GHF","GHF (ng ml -1)","mortality_data","Mortality per year","R_cat",halvveis) +
              facet_wrap(~mortality_type)
grid.arrange(GHF_mort_Rcat.plot, top=title.forward) 
dev.copy(png, paste0(path,"GHF_mortality_halfway_to_maxlength.png"), width=800,height=800) #Saving
dev.off()

#Making plot with mortality given OXF, colour R_cat
OXF_mort_Rcat.plot <- plot.simple(mort.ind2.df, "OXF","OXF (pg ml -1)","mortality_data","Mortality per year","R_cat",halvveis) +
              facet_wrap(~mortality_type)
grid.arrange(OXF_mort_Rcat.plot, top=title.forward) 
dev.copy(png, paste0(path,"OXF_mortality_halfway_to_maxlength.png"), width=800,height=800) #Saving
dev.off()

##Plotting mortality given environment, colour R_cat
E_mort_Rcat.plot <- plot.simple(mort.ind2.df, "environment","Environment","mortality_data","Mortality per year","R_cat",halvveis) +
              facet_wrap(~mortality_type)
grid.arrange(E_mort_Rcat.plot, top=title.forward) 
dev.copy(png, paste0(path,"environment_mortality_halfway_to_maxlength.png"), width=800,height=800) #Saving
dev.off()


#~ ##Plotting mortality given R_categories at different length in different environments
#~ E_cat_mort_len_E.plot <- ggplot(filter(ind2.df, environment!="NA"), aes(x=R_cat, y=M_size_foraging, colour=length)) +
#~                 geom_point() + facet_wrap(~as.factor(environment)) +
#~                 theme_bw(base_size=16) +
#~                 #scale_x_continuous(limits = R_cat.limits) +
#~                 labs(x="Reserve fullness (%)",y="Size dependent foraging mortality per year", colour="Length (cm)") +
#~                 ggtitle("Size dependent foraging mortality given reserve categories")
#~ grid.arrange(E_cat_mort_len_E.plot, top=title.forward)
#~ dev.copy(png, paste0(path,"R_cat_foraging_mortality_length_environment.png"), width=800,height=800) #Saving
#~ dev.off()



#Removing plots from memory
rm(list=ls(pattern=".plot"))
#Removing mean datasets from memory
rm(list=ls(pattern=".mean"))


  
#---------------------------------------------------------------------------------------------------------


##############################################
# Trade-off plots (also see mortality above) #
##############################################

#Making new folder for figures
path <- paste0(path.forward, "misc-plots/")
dir.create(path, showWarnings = TRUE)



#Making a new dataset from ind2.df with needed columns 
trade.df <- ind2.df

#Calculating 
trade.df$M_foraging <- trade.df$M_size_foraging/trade.df$M_size

##Plot OXF foraging mortality
trade.fm.fr.plot <- plot.simple(trade.df,"OXF","OXF", "M_foraging", "M_foraging", "environment", halvveis)

##Plot foraging_required intake against foraging mortality
trade.fm.fr.oxf.plot <- plot.simple(trade.df,"foraging_required","foraging_required", "M_foraging", "M_foraging", "OXF", halvveis)

grid.arrange(trade.fm.fr.plot, trade.fm.fr.oxf.plot, ncol=1, top=title.forward) 
dev.copy(png, paste0(path,"OXF-Mforaging-foraging_required.png"), width=800,height=800) #Saving
dev.off()



##Plotting O2 against mortality

#Making a new dataset from ind2.df
trade.df <- ind2.df

#Making a new variable with calculations
trade.df$calcO2 <- trade.df$O2_used / trade.df$O2max_THF

#Plotting and saving
tradeO2.plot <- plot.simple(trade.df,"calcO2","O2/O2max_THF", "M_size_O2", "M_size_O2", "R_cat", halvveis)
grid.arrange(tradeO2.plot, top=title.forward) 
dev.copy(png, paste0(path,"mortality_O2.png"), width=800,height=800) #Saving
dev.off()


##Plotting M_size_O2 given O2_used with THF level as colour
tradeO2_2.plot <- plot.simple(ind2.df,"O2_used","O2_used", "M_size_O2", "M_size_O2", "THF", halvveis)
grid.arrange(tradeO2_2.plot, top=title.forward) 
dev.copy(png, paste0(path,"mortality_O2_used_THF.png"), width=800,height=800) #Saving
dev.off()


##Plotting O2_used given environment with THF as colour
O2_envir.plot <-plot.simple(ind2.df, "environment", "Environment", "O2_used","O2 used", "THF", halvveis)
grid.arrange(O2_envir.plot, top=title.forward) 
dev.copy(png, paste0(path,"environment_O2_used_THF.png"), width=800,height=800) #Saving
dev.off()


##Plotting SMR against mortality

#Making a new dataset from ind2.df
trade.df <- ind2.df

#Calculating mortality 
trade.df$M_O2 <- trade.df$M_size_O2 / trade.df$M_size 

#Plotting and saving
trade.SMR.mort.plot <- plot.simple(trade.df, "SMR_THF", "SMR due to THF","M_O2", "Mortality due to 02 use", "THF", halvveis) 
grid.arrange(trade.SMR.mort.plot, top=title.forward) 
dev.copy(png, paste0(path,"mortality_SMR.png"), width=800,height=800) #Saving
dev.off()

##Plotting THF against SMR

#Making new dataset with SMR stacked on top of each other
trade.df <- gather(data=ind2.df, key=SMR_type, value=SMR, SMR_std:SMR_THF)

#Plotting and saving
tradeSMR1.plot <- plot.simple(trade.df,"THF","THF (ng ml -1)", "SMR", "SMR", "SMR_type", halvveis) #Plot with different colours for SMR_type
tradeSMR2.plot <- plot.simple(trade.df,"THF","THF (ng ml -1)", "SMR", "SMR", "R_cat", halvveis) #Plot with different colours for R_categories

grid.arrange(tradeSMR1.plot, tradeSMR2.plot, top=title.forward, nrow=2, ncol=1) 
dev.copy(png, paste0(path,"THF_SMR.png"), width=800,height=800) #Saving
dev.off()

#Removing trade.df from memory
rm(trade.df)

##Plotting THF against O2 use

#Making new dataset with O2max_THF and O2_used stacked on top of each other
trade.df <- gather(data=ind2.df, key=O2_type, value=O2, O2max_THF:O2_used)

#Plotting and saving
tradeO21.plot <- plot.simple(trade.df,"THF","THF (ng ml -1)", "O2", "Oxygen use", "O2_type", halvveis) #Plot with different colours for O2_type
tradeO22.plot <- plot.simple(trade.df,"THF","THF (ng ml -1)", "O2", "Oxygen use", "R_cat", halvveis) #Plot with different colours for R_categories
grid.arrange(tradeO21.plot,tradeO22.plot, top=title.forward, ncol=1) 
dev.copy(png, paste0(path,"THF_O2.png"), width=800,height=800) #Saving
dev.off()

##Plotting environment against O2 use, mortality and THF
tradeO21e.plot <- plot.simple(trade.df,"environment","Environment", "O2", "Oxygen use", "O2_type", halvveis) #Plot with different colours for O2_type
tradeO22e.plot <- plot.simple(trade.df,"environment","Environment", "O2", "Oxygen use", "THF", halvveis) + facet_wrap(~O2_type,ncol=1) #Plot with different colours for THF
tradeO23e.plot <- plot.simple(ind2.df,"environment","Environment", "M_size_O2", "Size dependent O2 mortality (per year)", "THF", halvveis) #O2 mortality with different colours for THF
tradeO24e.plot <- plot.simple(ind2.df,"environment","Environment", "THF", "THF (ng ml-1)", "R_cat", halvveis) #Plot with different colours for R_categories
grid.arrange(tradeO21e.plot,tradeO22e.plot, tradeO23e.plot, tradeO24e.plot, top=title.forward, nrow=4, ncol=1) 
dev.copy(png, paste0(path,"Environment_O2_THF.png"), width=800,height=1600) #Saving
dev.off()

#Removing trade.df from memory
rm(trade.df)

##Plotting THF against survival 
tradeT_surnext.plot <- plot.simple(ind2.df,"THF","THF (ng ml -1)","survival_next", "Survival", "environment", halvveis)
tradeT_surstep.plot <- plot.simple(ind2.df,"THF","THF (ng ml -1)","survival_step", "Probability of surviving the current timestep", "environment", halvveis)
grid.arrange(tradeT_surnext.plot,tradeT_surstep.plot, top=title.forward, nrow=2, ncol=1) 
dev.copy(png, paste0(path,"THF_survival.png"), width=800,height=800) #Saving
dev.off()


##Plotting THF against intake and foraging_required with environment as colour
tradeT_fr_i.plot <- plot.simple(ind2.df,"THF","THF (ng ml -1)", "foraging_required", "Foraging required (in multiples of SMR)", "environment", halvveis) #Plot with different colours for R_categories
tradeT_i_fr.plot <- plot.simple(ind2.df,"THF","THF (ng ml -1)", "intake", "Intake (J timestep -1)", "environment", halvveis) #Plot with different colours for R_categories
grid.arrange(tradeT_fr_i.plot,tradeT_i_fr.plot, top=title.forward, nrow=2, ncol=1) 
dev.copy(png, paste0(path,"THF_intake_foraging_required.png"), width=800,height=800) #Saving
dev.off()


##Plotting THF against growth with environment as colour
tradeT_grow.plot <- plot.simple(ind2.df,"THF","THF (ng ml -1)", "growth", "Growth (kg)", "environment", halvveis)
grid.arrange(tradeT_grow.plot, top=title.forward) 
dev.copy(png, paste0(path,"THF_growth.png"), width=800,height=800) #Saving
dev.off()


##Plotting environment against intake and foraging_required with THF as colour
tradeT_fr_THF.plot <- plot.simple(ind2.df,"environment","Environment", "foraging_required", "Foraging required (in multiples of SMR)", "THF", halvveis) #Plot with different colours for R_categories
tradeT_i_THF.plot <- plot.simple(ind2.df,"environment","Environment", "intake", "Intake (J timestep -1)", "THF", halvveis) #Plot with different colours for R_categories
grid.arrange(tradeT_fr_THF.plot,tradeT_i_THF.plot, top=title.forward, nrow=2, ncol=1) 
dev.copy(png, paste0(path,"environment_intake_foraging_required_THF.png"), width=800,height=800) #Saving
dev.off()


##Plotting GHF against growth

#Plotting and saving
tradegrowth.plot <- plot.simple(ind2.df,"GHF","GHF (ng ml -1)", "growth", "Growth (kg)", "R_cat", halvveis)
grid.arrange(tradegrowth.plot, top=title.forward) 
dev.copy(png, paste0(path,"GHF_growth.png"), width=800,height=800) #Saving
dev.off()


##Plotting GHF against conversion costs

#Plotting and saving
tradeconversioncost1.plot <- plot.simple(ind2.df,"GHF","GHF (ng ml -1)", "conversion_cost_to_growth", "Conversion cost to growth rescaled to requirement (J)", "R_cat", halvveis)
tradeconversioncost2.plot <- plot.simple(ind2.df,"GHF","GHF (ng ml -1)", "conversion_costs_via_reserves", "Conversion costs from intake to reserves, and from reserves to metabolism (J)", "R_cat", halvveis)
grid.arrange(tradeconversioncost1.plot, tradeconversioncost2.plot, top=title.forward,nrow=2, ncol=1) 
dev.copy(png, paste0(path,"GHF_conversioncosts.png"), width=800,height=800) #Saving
dev.off()


##Plotting GHF against survival 
tradeG_surnext.plot <- plot.simple(ind2.df,"GHF","GHF (ng ml -1)","survival_next", "Survival", "environment", halvveis)
tradeG_surstep.plot <- plot.simple(ind2.df,"GHF","GHF (ng ml -1)","survival_step", "Probability of surviving the current timestep", "environment", halvveis)
grid.arrange(tradeG_surnext.plot,tradeG_surstep.plot, top=title.forward, nrow=2, ncol=1) 
dev.copy(png, paste0(path,"GHF_survival.png"), width=800,height=800) #Saving
dev.off()


##Plotting OXF against intake

#Plotting and saving
tradeintake.plot <- plot.simple(ind2.df,"OXF","OXF (pg ml -1)", "intake", "Intake (J timestep -1)", "R_cat", halvveis)
grid.arrange(tradeintake.plot, top=title.forward) 
dev.copy(png, paste0(path,"OXF_intake.png"), width=800,height=800) #Saving
dev.off()


##Plotting OXF against target_intake

#Plotting and saving
tradeforaging.plot <- plot.simple(ind2.df,"OXF","OXF (pg ml -1)", "target_intake", "Target intake given OXF (Multiples of SMR)", "R_cat", halvveis)
grid.arrange(tradeforaging.plot, top=title.forward) 
dev.copy(png, paste0(path,"OXF_target_intake.png"), width=800,height=800) #Saving
dev.off()


##Plotting OXF against foraging

#Plotting and saving
tradeforaging.plot <- plot.simple(ind2.df,"OXF","OXF (pg ml -1)", "foraging_required", "Foraging Intensity Required (Multiples of SMR)", "R_cat", halvveis)
grid.arrange(tradeforaging.plot, top=title.forward) 
dev.copy(png, paste0(path,"OXF_foraging_required.png"), width=800,height=800) #Saving
dev.off()


##Plotting OXF against foraging costs

#Plotting and saving
tradeforagingcosts.plot <- plot.simple(ind2.df,"OXF","OXF (pg ml -1)", "foraging_cost", "Foraging costs (J timestep -1)", "R_cat", halvveis)
grid.arrange(tradeforagingcosts.plot, top=title.forward) 
dev.copy(png, paste0(path,"OXF_foragingcosts.png"), width=800,height=800) #Saving
dev.off()


##Plotting OXF against survival 
tradeO_surnext.plot <- plot.simple(ind2.df,"OXF","OXF (pg ml -1)","survival_next", "Survival", "environment", halvveis)
tradeO_surstep.plot <- plot.simple(ind2.df,"OXF","OXF (pg ml -1)","survival_step", "Probability of surviving the current timestep", "environment", halvveis)
grid.arrange(tradeO_surnext.plot,tradeO_surstep.plot, top=title.forward, nrow=2, ncol=1) 
dev.copy(png, paste0(path,"OXF_survival.png"), width=800,height=800) #Saving
dev.off()


##Plotting environment against survival 
tradeG_surnext.plot <- plot.simple(ind2.df,"environment","Environment","survival_next", "Survival", "R_cat", halvveis)
tradeG_surstep.plot <- plot.simple(ind2.df,"environment","Environment","survival_step", "Probability of surviving the current timestep", "R_cat", halvveis)
grid.arrange(tradeG_surnext.plot,tradeG_surstep.plot, top=title.forward, nrow=2, ncol=1) 
dev.copy(png, paste0(path,"environment_survival.png"), width=800,height=800) #Saving
dev.off()


##Plotting foraging_required against environment
tradeforagingcosts.plot <- plot.simple(ind2.df,"environment","Environment", "foraging_required", "Foraging Intensity Required (Multiples of SMR)", "target_intake", halvveis)
grid.arrange(tradeforagingcosts.plot, top=title.forward) 
dev.copy(png, paste0(path,"environment_foraging_required.png"), width=800,height=800) #Saving
dev.off()


##Plotting surplus_before_growth / weight * 1000 against environment
tradeE_sbg.plot <- plot.simple(ind2.df,"environment","Environment", "surplus_before_growth/(weight*1000)", "Surplus before growth (J g -1 timestep -1)", "R_cat", halvveis)
grid.arrange(tradeE_sbg.plot, top=title.forward) 
dev.copy(png, paste0(path,"environment_surplus_before_growth.png"), width=800,height=800) #Saving
dev.off()


##Plotting surplus_before_growth / weight * 1000 against environment at different lengths
tradeE_sbg1.plot <- plot.simple(ind2.df,"environment","Environment", "surplus_before_growth/(weight*1000)", "Surplus before growth (J g -1 timestep -1)", "R_cat", length.limits[1]) + 
                    scale_y_continuous(limits = c(min(ind2.df$surplus_before_growth/(ind2.df$weight*1000), na.rm=TRUE),max(ind2.df$surplus_before_growth/(ind2.df$weight*1000), na.rm=TRUE))) +
                    geom_hline(yintercept = 0)
tradeE_sbg2.plot <- plot.simple(ind2.df,"environment","Environment", "surplus_before_growth/(weight*1000)", "Surplus before growth (J g -1 timestep -1)", "R_cat", round(length.limits[1]+((length.limits[2]-length.limits[1])/3))) +
                    scale_y_continuous(limits = c(min(ind2.df$surplus_before_growth/(ind2.df$weight*1000), na.rm=TRUE),max(ind2.df$surplus_before_growth/(ind2.df$weight*1000), na.rm=TRUE))) +
                    geom_hline(yintercept = 0)
tradeE_sbg3.plot <- plot.simple(ind2.df,"environment","Environment", "surplus_before_growth/(weight*1000)", "Surplus before growth (J g -1 timestep -1)", "R_cat", round(length.limits[1]+((length.limits[2]-length.limits[1])/3)*2)) +
                    scale_y_continuous(limits = c(min(ind2.df$surplus_before_growth/(ind2.df$weight*1000), na.rm=TRUE),max(ind2.df$surplus_before_growth/(ind2.df$weight*1000), na.rm=TRUE))) +
                    geom_hline(yintercept = 0)
tradeE_sbg4.plot <- plot.simple(ind2.df,"environment","Environment", "surplus_before_growth/(weight*1000)", "Surplus before growth (J g -1 timestep -1)", "R_cat", length.limits[2]) +
                    scale_y_continuous(limits = c(min(ind2.df$surplus_before_growth/(ind2.df$weight*1000), na.rm=TRUE),max(ind2.df$surplus_before_growth/(ind2.df$weight*1000), na.rm=TRUE))) +
                    geom_hline(yintercept = 0)
grid.arrange(tradeE_sbg1.plot,tradeE_sbg2.plot,tradeE_sbg3.plot,tradeE_sbg4.plot, top=title.forward, nrow=2, ncol=2) 
dev.copy(png, paste0(path,"environment_surplus_before_growth_length.png"), width=800,height=800) #Saving
dev.off()


##Plotting change in reserves against environment
tradeE_rd.plot <- plot.simple(ind2.df,"environment","Environment", "reserve_diff", "Change in reserves (J timestep -1)", "R_cat", halvveis) +
                  scale_y_continuous(labels=function(n){format(n, scientific = FALSE)})
grid.arrange(tradeE_rd.plot, top=title.forward) 
dev.copy(png, paste0(path,"environment_reserve_diff.png"), width=800,height=800) #Saving
dev.off()


##Plotting surplus_before_growth against environment at different lengths
tradeE_rd1.plot <- plot.simple(ind2.df,"environment","Environment", "reserve_diff", "Change in reserves (J timestep -1)", "R_cat", length.limits[1]) +
                   scale_y_continuous(limits = c(min(ind2.df$reserve_diff, na.rm=TRUE), max(ind2.df$reserve_diff, na.rm=TRUE)),labels=function(n){format(n, scientific = FALSE)}) +
                   geom_hline(yintercept = 0)
tradeE_rd2.plot <- plot.simple(ind2.df,"environment","Environment", "reserve_diff", "Change in reserves (J timestep -1)", "R_cat", round(length.limits[1]+((length.limits[2]-length.limits[1])/3))) +
                   scale_y_continuous(limits = c(min(ind2.df$reserve_diff, na.rm=TRUE),max(ind2.df$reserve_diff, na.rm=TRUE)),labels=function(n){format(n, scientific = FALSE)}) +
                  geom_hline(yintercept = 0)
tradeE_rd3.plot <- plot.simple(ind2.df,"environment","Environment", "reserve_diff", "Change in reserves (J timestep -1)", "R_cat", round(length.limits[1]+((length.limits[2]-length.limits[1])/3)*2)) +
                   scale_y_continuous(limits = c(min(ind2.df$reserve_diff, na.rm=TRUE),max(ind2.df$reserve_diff, na.rm=TRUE)),labels=function(n){format(n, scientific = FALSE)}) +
                  geom_hline(yintercept = 0)
tradeE_rd4.plot <- plot.simple(ind2.df,"environment","Environment", "reserve_diff", "Change in reserves (J timestep -1)", "R_cat", length.limits[2]) +
                   scale_y_continuous(limits = c(min(ind2.df$reserve_diff, na.rm=TRUE),max(ind2.df$reserve_diff, na.rm=TRUE)),labels=function(n){format(n, scientific = FALSE)}) +
                  geom_hline(yintercept = 0)
grid.arrange(tradeE_rd1.plot,tradeE_rd2.plot,tradeE_rd3.plot,tradeE_rd4.plot, top=title.forward, nrow=2, ncol=2) 
dev.copy(png, paste0(path,"environment_reserve_diff_length.png"), width=800,height=800) #Saving
dev.off()


##Plot Environment and mean intake and SMR etc.

#Transforming ProbAllE array into a dataset
ProbOfE.df <- melt(proballe.mx, varnames="environment") 
ProbOfE.df2 <- with(ProbOfE.df, approx(environment, value, xout = seq(min(environment),max(environment),length = 1000))) %>%
as.data.frame() %>% select(environment = x, value = y) 

#Calculating mean intake given environment and making it into a new dataset
intake.envir.mean <- aggregate(intake ~ environment_rounded, filter(ind2.df, length_rounded==halvveis), mean)

#Calculating mean SMR_std given environment and making it into a new dataset
SMRstd.envir.mean <- aggregate(SMR_std ~ environment_rounded, filter(ind2.df, length_rounded==halvveis), mean)

#~ x11()
par(mar = c(5,5,2,5), cex.axis=1.5, cex.lab=1.5)
plot(value~environment, type="l", lwd=1.5, col="plum", xlab="",ylab="", axes=FALSE, main="Length=20", data=ProbOfE.df) 
par(new=TRUE)
plot(intake~environment_rounded, cex=1.5, lwd=1.5, col="deeppink", xlab="Rounded Environment", ylab="", axes=FALSE, data=intake.envir.mean) 
abline(abline(h = mean(filter(ind2.df, length_rounded==halvveis)$intake, na.rm=TRUE), lwd=1.5, lty=2, col="deeppink"))
axis(1)
axis(2, col="deeppink")
mtext(side=2, line=3, cex=1.5, "Intake (J timestep -1)")
par(new=TRUE)
plot(SMR_std~environment_rounded, cex=1.5, lwd=1.5, col="darkturquoise", xlab="", ylab="", axes=FALSE,data=SMRstd.envir.mean)
axis(4, col="darkturquoise")
mtext(side=4, line=3, cex=1.5,"SMR standard")
box()
par(new=TRUE, xpd=TRUE)
legend("bottom", inset=0.05, box.lty=2, cex=1,
        legend=c("Probability curve","Total Mean Intake","Mean Intake", "Mean SMR standard"),
        lty=c(1,2,0,0), lwd=1.5, pch=c(NA,NA,1,1), col=c("plum","deeppink","deeppink", "darkturquoise"))
dev.copy(png, paste0(path,"intake_SMR_environment.png"), width=800,height=800) #Saving
dev.off()
#~ graphics.off()

##
#~ test <- ggplot(intake.envir.mean, aes(y=intake, x=environment_rounded)) +
#~         geom_point() +
#~         geom_hline(yintercept=mean(filter(ind2.df, length_rounded==halvveis)$intake)) +
#~         geom_rug(aes(x = environment, alpha = value), data  = ProbOfE.df2, inherit.aes = FALSE)
#~ test
#~ test <- ggplot(intake.envir.mean, aes(y=intake, x=environment_rounded)) +
#~         geom_point() +
#~         geom_hline(yintercept=mean(filter(ind2.df, length_rounded==halvveis)$intake)) +
#~         geom_line(aes(x = environment, y = value), data  = ProbOfE.df2%>%mutate(value = value/max(value) * diff(range(intake.envir.mean$intake, na.rm=TRUE)) + min(intake.envir.mean$intake,na.rm=TRUE)), inherit.aes = FALSE)
#~ test



#Removing plots from memory
rm(list=ls(pattern=".plot"))


  
#---------------------------------------------------------------------------------------------------------


####################
# Individual plots #
####################

#Script for making individual plots (see file names plotting_individuals.r)
#Only run if there is items in the singleind.list
if(length(singleind.list) > 0){
    
  #Making new dataset with the chosen individuals in 
  ind.ind.df <- filter(ind2.df, individual %in% singleind.list)



#~   ##Plots with individuals that were chosen in singleind.list
#~   ## with length against every variable in the dataset

#~   #Making new subfolder
#~   path <- paste0(path.forward , "individual/")
#~   dir.create(path, showWarnings = TRUE)

#~   #Counter for loop
#~   count = 1

#~   for (type in ind.var.names){
    
#~     #Getting the title from the ind.head.names array
#~     title.plot <- ind.head.names[count]
    
#~     #Plotting
#~     X11() #Opening plot in new window
#~     ind.plot <- ggplot(ind.ind.df, aes_string(x="length", y=type, linetype="as.factor(individual)")) +
#~       geom_line() +
#~       geom_point(aes(colour=environment)) + #Points with colour depending on environment
#~       theme_bw(base_size=16) + 
#~       scale_x_continuous(limits = length.limits) +
#~       scale_y_continuous(labels=function(n){format(n, scientific = FALSE)}) + 
#~       ggtitle(title.plot) +
#~       labs(x="Length (cm)", y=ind.head.names[count], linetype="Growth Speed", colour="Environment")
#~     #Printing plot to screen  
#~     grid.arrange(ind.plot, top=title.forward, nrow=1, ncol=1)
    
#~     #Saving plot to file 
#~     dev.copy(png, paste0(path, "ind_", ind.head.names[count], ".png"), width=800,height=800)
#~     dev.off()
    
#~     #Counting number of loops
#~     count = count + 1
#~   }

#~   #Syntax for closing all plots on screen
#~   graphics.off()
  
  
  
  ##Mortality for single individuals
  ###################################
  
  #Making folder
  path <- paste0(path.forward, "mortality/")
  dir.create(path, showWarnings = TRUE)
  
  #Loops that plots for every individual in the singleind.list
  for (singleind in singleind.list){
    #If singleind is set to "random", pick a random individual from the dataset
    #if(singleind == "random"){singleind <- sample(ind2.df$individual, 1, replace=TRUE)}
    if(singleind != "random"){singleind <- as.numeric(singleind)}

    #Making a new dataset for the chosen individual based on mort.ind.df
    mort.single.ind.df <- filter(mort.ind.df, individual == singleind)
    
    #Saving the survival (survival_next) at the end of the individuals life
    survival.var <- filter(ind2.df, individual == singleind, timestep==(max(mort.single.ind.df$last_timestep)+1))$survival_next
    
    #Saving the growth speed group of the individual
    growth_speed.var <- filter(ind2.df, individual==singleind, timestep==1)$growth_speed
    
    #Saving the end status of the individual as a variable
    if(mort.single.ind.df$end_status[1] == 1){end.status.var="Dead"} 
    if(mort.single.ind.df$end_status[1] == 2){end.status.var="To slow"} 
    if(mort.single.ind.df$end_status[1] == 3){end.status.var="Reached max length"} 
    
    #Removing data for sizeindependent if the model is run without size independent mortality
    if (Msizeindependent.check == "without"){
      mort.single.ind.df <- filter(mort.single.ind.df, mortality_type != "M_sizeindependent")
    }
    
    #Calculating mean mortality given length and mortality type for those rare cases 
    # when length do not change that much from one timestep to the next, which ruins the figures
    mort.single.ind.mean <- aggregate(mortality_data ~ length + mortality_type, mort.single.ind.df, mean)
    
    #Calculating mean environment given length and mortality type for those rare cases 
    # when length do not change that much from one timestep to the next, which ruins the figures
    mort.single.ind.mean2 <- aggregate(environment ~ length + mortality_type, mort.single.ind.df, mean)

    #Mortality plot against length
    x11()
    mort_length.plot <- ggplot(mort.single.ind.mean, aes(x=length, y=mortality_data, fill=mortality_type)) +
                        geom_area(colour="black", size=.2, alpha=.4) +
                        geom_line(aes(y = mort.single.ind.mean2$environment/(max(ind2.df$environment,na.rm=TRUE)/1.5)), linetype="dotted", size=.2, alpha=.4) +
                        theme_bw(base_size=16) + 
                        scale_x_continuous(limits = length.limits) +
                        scale_y_continuous(limits = c(0, 1.5), sec.axis = sec_axis(~.*((max(ind2.df$environment,na.rm=TRUE)/1.5)), name = "Environment")) +
                        labs(x="Length", y="Mortality (per year)", fill="Mortality type") +
                        ggtitle(paste0("Individual: ", singleind),
                          subtitle=paste("Last timestep:", max(mort.single.ind.df$last_timestep), " Survival:", round(survival.var,4),
                          "\nGrowth speed:", growth_speed.var, " End status:", end.status.var))
    grid.arrange(mort_length.plot, top=title.forward, nrow=1, ncol=1)

    #Saving and closing plot
    dev.copy(png, paste0(path, singleind,"-mortality_length.png"), width=800,height=800)
    dev.off()
    graphics.off()
    
    #Remove datasets from memory
    rm(mort.single.ind.df, mort.single.ind.mean, mort.single.ind.mean2)
    
    
  }

}



#Removing plots from memory
rm(list=ls(pattern=".plot"))



###Plots with mean of the groups based on which environment they started in

#~ #Making new subfolder
#~ path <- paste0(path.forward, "start environment/")
#~ dir.create(path, showWarnings = TRUE)

#~ ind2.df$start_environment2 = as.factor(ind2.df$start_environment)

#~ #Counter for loop
#~ count = 1

#~ #Plotting (see above)
#~ for (type in ind.var.names){
  
#~   #Setting title if "" is chosen
#~   if(title.forward==""){title.plot <- ind.head.names[count]
#~   }else{title.plot <- title.forward}
  
#~   X11() #Opening plot in new window
#~   ind.plot <- ggplot(ind2.df, aes_string(x="length_rounded", y=type, colour = "start_environment2", fill = "start_environment2")) +
#~     stat_summary(geom="ribbon", fun.ymin = min, fun.ymax = max, colour=NA, alpha = 0.1) +
#~     stat_summary(fun.y=mean, geom="line", aes(group = start_environment2)) +
#~     theme_bw(base_size=16) + 
#~     scale_fill_discrete(guide=FALSE) + 
#~     scale_x_continuous(limits = length.limits) +
#~     ggtitle(title.plot) +
#~     labs(x="Length (cm)", y=ind.head.names[count], colour="Start Environment")
    
#~   grid.arrange(ind.plot, nrow=1, ncol=1)
  
#~   dev.copy(png, paste0(path, "se_", ind.head.names[count], ".png"), width=800,height=800)
#~   dev.off()
  
#~   #Counting number of loops
#~   count = count + 1
#~ }
#~ #Syntax for closing all plots on screen
#~ graphics.off()


#~ #Makes a simple plot for hormone levels given reserve category
#~ par(mfrow=c(2,2))
#~ plot(OXF~R_cat, ylim=OXF.limits, main=title.forward, data=ind2.df) #OXF
#~ plot(GHF~R_cat, ylim=GHF.limits, data=ind2.df)  #GHF
#~ plot(THF~R_cat, ylim=THF.limits, data=ind2.df)  #THF

#~ #Saving plot
#~ dev.copy(png, paste0(path, "hormones~R_cat.png"), width=800,height=800)
#~ dev.off()

#~ #Closing and resetting plotting window
#~ graphics.off()


#~ #Makes a simple plot for hormone levels given reserves
#~ par(mfrow=c(2,2))
#~ plot(OXF~reserves, ylim=OXF.limits, main=title.forward, data=ind2.df) #OXF
#~ plot(GHF~reserves, ylim=GHF.limits, data=ind2.df)  #GHF
#~ plot(THF~reserves, ylim=THF.limits, data=ind2.df)  #THF

#~ #Saving plot
#~ dev.copy(png, paste0(path, "hormones~reserves.png"), width=800,height=800)
#~ dev.off()

#~ #Closing and resetting plotting window
#~ graphics.off()


#~ #Makes a simple plot for hormone levels given length
#~ par(mfrow=c(2,2))
#~ plot(OXF~length, ylim=OXF.limits, main=title.forward, data=ind2.df) #OXF
#~ plot(GHF~length, ylim=GHF.limits, data=ind2.df)  #GHF
#~ plot(THF~length, ylim=THF.limits, data=ind2.df)  #THF

#~ #Saving plot
#~ dev.copy(png, paste0(path, "hormones~length.png"), width=800,height=800)
#~ dev.off()

#~ #Closing and resetting plotting window
#~ graphics.off()


#~ #Function for plotting mortality plots with hormones on xaxis OLD
#~ ggplot2_h_mort <- function(hormone, h.limits, llength) {
  
#~   if (hormone=="OXF"){
#~     hm.mean.df <- aggregate(mortality_data ~ OXF + mortality_type, filter(mort.ind.df, length_rounded==llength), mean)
#~   }else if (hormone=="GHF"){
#~     hm.mean.df <- aggregate(mortality_data ~ GHF + mortality_type, filter(mort.ind.df, length_rounded==llength), mean)
#~   }else if (hormone=="THF"){
#~     hm.mean.df <- aggregate(mortality_data ~ THF + mortality_type, filter(mort.ind.df, length_rounded==llength), mean)
#~   }
  
#~   #Making plot 
#~   h.m.plot <- ggplot(hm.mean.df, aes_string(x=hormone, y="mortality_data", fill="mortality_type")) +
#~                     geom_area(colour="black", size=.2, alpha=.4) +
#~                     theme_bw(base_size=16) + 
#~                     scale_x_continuous(limits = h.limits) +
#~                     scale_y_continuous(limits = c(0, 2.5)) +
#~                     labs(x=hormone, y="Mortality per year", fill="Mortality type") +
#~                     ggtitle(paste("Length =", llength))
                    
#~   #Returning only plot from function
#~   return(h.m.plot)
#~ }


##Simple dataset for testing
#ind2.df <- data.frame(timestep=c(1,1,2,2,3,3,4,4))
#ind2.df$individual <- c(1,2)
#ind2.df$length <- c(2,3,3,4,5,4,0,4)
#ind2.df$end_status <- c(0,0,0,0,3,0,3,2)
#ind2.df$environment <- c(1,2,3,4,5,6,7,8)
#ind2.df

#Printing summary tables again just for an overview in the end:
x11()
grid.arrange(end_status.table, growth_speed.table, top=title.forward)


#Finding the SMR exponent
#~ s_exp = 0.7 #Exponent set in parameter.f90
#~ w_coeff = 0.07 #Weight used in SMR_coeff calculations
#~ SMR_coeff = log((((exp(-5.43)  * 434. * 24.)*(0.001^(0.8)) * 7))*(w_coeff^(s_exp-0.8)))
#~ test.lm <-lm(log(ind2.df$SMR_THF) ~ log(ind2.df$weight))
#~ summary(test.lm)
#~ ggplot(ind2.df, aes(y=log(SMR_THF), x=log(weight), colour=length)) +
#~   geom_point() +
#~   geom_smooth(method="lm", formula = y ~ x, colour="red")

#~ test.df <- filter(ind2.df, environment_rounded==1)
#~ test.lm <-lm(log(test.df$SMR_THF) ~ log(test.df$weight))
#~ summary(test.lm)

