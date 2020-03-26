#*******************************************************************************
# SVN Version ID: $Id$
#*******************************************************************************

#Requires that you have already run the main model and bisection_MO2.f90

#Syntax for plotting figures in R

#Clears workspace to get a fresh start (avoids errors if you are plotting several results in different
# folders after each other)
rm(list=ls())

#Install libraries needed (you only need to do this once)
# might take a few minutes to install, so you have time to get a cup of coffee or tea while waiting :)
#install.packages(c("lattice","corrplot","reshape2","tidyverse","gridExtra","GGally","directlabels","ggcorrplot","RColorBrewer","cowplot", "scales"), dependencies=TRUE)

#Importing libraries
lapply(c("lattice","corrplot","reshape2","tidyverse","gridExtra","GGally","directlabels","ggcorrplot","RColorBrewer"), require, character.only = TRUE)

#Select how many individuals you want to plot in forward (write "max" if you want all individuals)
# if you don't want max, you can write for example 10 000 without " "
max.ind = "max"

#Information from the model where c(min,max)
length.limits   <- c(10,30)
R_cat.limits    <- c(0,100)
OXF.limits      <- c(0,1500)
GHF.limits      <- c(0,200)
THF.limits      <- c(0,5)

#Output folder:
# See under ind.bin

#Functions for quantiles used in plots
lower.quant = function(z) {quantile(z,0.0)}
upper.quant = function(z) {quantile(z,1.0)}

#Function for converting from one min-max scale to another min-max scale
convert.scale <- function(old_value, old_min, old_max, new_min, new_max) {
  new_value <- ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min
}

#Function of importing binary files containing reals
# following http://www.ats.ucla.edu/stat/r/faq/read_binary.htm
import.binary <- function(file_path) {
  bin <- file(file_path, "rb")
  ranking <- readBin(bin,n=1,integer())
  dim_length <- readBin(bin, n = ranking, integer())
  bin.data <- readBin(bin, n = prod(dim_length), numeric())
  bin.data <- array(bin.data,dim=dim_length)
}


#---------------------------------------------------------------------------------------------------------

######################################
#List of all names used in ind array #
######################################

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


#---------------------------------------------------------------------------------------------------------


####################################
# importing data from static model #
####################################

#Script for importing ind.bin files from the static model and get row
# from dataset closest to length = 20

get.static.row <- function(file_path){
  
  ind.tp <- import.binary(file_path)
  ind.tp <- melt(ind.tp, varnames=c("type.of.data", "timestep","individual"))

  ind.tp <- filter(ind.tp, individual==1 & timestep!=max(ind.tp$timestep)) %>%
            mutate(type.of.data = ind.var.names[type.of.data]) %>%
            spread(key = type.of.data, value = value) %>%
            filter(., length==.$length[which.min(abs(.$length-20))])
}                

#Path to folder with static simulations
path.static <- "./AHAsvn/tegsvn/hormones_v2/branches_camilla/testingbranch/Til Paper 1.1"


##Simulation environment 1
##########################
ind.static.df             <- get.static.row(paste0(path.static, "/1 environment_scaling = 0.36 & t_max = 284/Results/ind.bin"))
ind.static.df$environment <- 1

##Simulation environment 2
##########################
ind2.tp                   <- get.static.row(paste0(path.static, "/2 environment_scaling = 0.4763636/Results/ind.bin")) 
ind2.tp$environment       <- 2
ind.static.df             <- bind_rows(ind.static.df, ind2.tp)

##Simulation environment 3
##########################
ind3.tp                   <- get.static.row(paste0(path.static, "/3 environment_scaling = 0.5927273/Results/ind.bin")) 
ind3.tp$environment       <- 3
ind.static.df             <- bind_rows(ind.static.df, ind3.tp)

##Simulation environment 4
##########################
ind4.tp                   <- get.static.row(paste0(path.static, "/4 environment_scaling = 0.7090909/Results/ind.bin")) 
ind4.tp$environment       <- 4
ind.static.df             <- bind_rows(ind.static.df, ind4.tp)

##Simulation environment 5
##########################
ind5.tp                   <- get.static.row(paste0(path.static, "/5 environment_scaling = 0.8254545/Results/ind.bin")) 
ind5.tp$environment       <- 5
ind.static.df             <- bind_rows(ind.static.df, ind5.tp)

##Simulation environment 6
##########################
ind6.tp                   <- get.static.row(paste0(path.static, "/6 environment_scaling = 0.9418182/Results/ind.bin")) 
ind6.tp$environment       <- 6
ind.static.df             <- bind_rows(ind.static.df, ind6.tp)

##Simulation environment 7
##########################
ind7.tp                   <- get.static.row(paste0(path.static, "/7 environment_scaling = 1.0581818/Results/ind.bin")) 
ind7.tp$environment       <- 7
ind.static.df             <- bind_rows(ind.static.df, ind7.tp)

##Simulation environment 8
##########################
ind8.tp                   <- get.static.row(paste0(path.static, "/8 environment_scaling = 1.1745455/Results/ind.bin")) 
ind8.tp$environment       <- 8
ind.static.df             <- bind_rows(ind.static.df, ind8.tp)

##Simulation environment 9
##########################
ind9.tp                   <- get.static.row(paste0(path.static, "/9 environment_scaling = 1.2909091/Results/ind.bin")) 
ind9.tp$environment       <- 9
ind.static.df             <- bind_rows(ind.static.df, ind9.tp)

##Simulation environment 10
###########################
ind10.tp                  <- get.static.row(paste0(path.static, "/10 environment_scaling = 1.4072727/Results/ind.bin")) 
ind10.tp$environment      <- 10
ind.static.df             <- bind_rows(ind.static.df, ind10.tp)

##Simulation environment 11
###########################
ind11.tp                  <- get.static.row(paste0(path.static, "/11 environment_scaling = 1.5236364/Results/ind.bin")) 
ind11.tp$environment      <- 11
ind.static.df             <- bind_rows(ind.static.df, ind11.tp)

##Simulation environment 12
###########################
ind12.tp                  <- get.static.row(paste0(path.static, "/12 environment_scaling = 1.64/Results/ind.bin")) 
ind12.tp$environment      <- 12
ind.static.df             <- bind_rows(ind.static.df, ind12.tp)

#Removing ind.tp from memory
rm(list=ls(pattern=".tp"))


#---------------------------------------------------------------------------------------------------------

# Importing data containing from experimental simulations
ind_ex.mx <- import.binary("./AHAsvn/tegsvn/hormones_v2/stochastic_trunk/Results/Lagrede tester uke 16/1 individuals with fixed environments/Results/ind.bin")
ind_ex.df <- melt(ind_ex.mx, varnames=c("type.of.data", "timestep","individual"))
ind_ex.df <- ind_ex.df %>%
              mutate(type.of.data = ind.var.names[type.of.data]) %>%
              spread(key = type.of.data, value = value) 
              
ind_ex.mx <- import.binary("./AHAsvn/tegsvn/hormones_v2/stochastic_trunk/Results/Lagrede tester uke 16/2 individuals with fixed environments/Results/ind.bin")
ind_ex2.df <- melt(ind_ex.mx, varnames=c("type.of.data", "timestep","individual"))
ind_ex2.df <- ind_ex2.df %>%
              mutate(type.of.data = ind.var.names[type.of.data]) %>%
              spread(key = type.of.data, value = value) 

#Removing matrix from memory
rm(list=ls(pattern="ind_ex.mx"))

#---------------------------------------------------------------------------------------------------------


###########
# ind.bin #
###########

#Setting working directory: Where you keep your result files (add your own here)
setwd("./AHAsvn/tegsvn/hormones_v2/stochastic_trunk/Results/Lagrede tester uke 16/0 Ny Default/Results/")
getwd() #Getting working directory

#Making new folder for forward figures
path <- "./Figures/paper-figures/"
dir.create(path, showWarnings = TRUE)

#Importing ind.bin to R
ind.mx <- import.binary("./ind.bin")

#Melting 3D array into 2D data frame
ind.df <- melt(ind.mx, varnames=c("type.of.data", "timestep","individual"))
rm(ind.mx) #Removing array from memory 

#Turning the ID of individuals and type.of.data from real into a factor
ind.df$type.of.data <- as.factor(ind.df$type.of.data) 

#Selecting how many individuals you want to plot if you have chosen a limit
if(max.ind != "max"){ind.df <- filter(ind.df, individual >= 1 & individual <= max.ind)}

#Variable with the length for halfway to maxlength
halvveis <- round(length.limits[1]+((length.limits[2]-length.limits[1])/2))

#Making new dataset where each data type has it own column with names based on ind.var.names
ind2.df <- ind.df  %>%
            mutate(type.of.data = ind.var.names[type.of.data]) %>%
            spread(key = type.of.data, value = value) %>%
            filter(length!=-1000)
rm(ind.df)

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

#Replacing all -1000 with NA for better plots for both static and stochastic dataset
ind2.df[ind2.df==-1000]             <- NA
ind.static.df[ind.static.df==-1000] <- NA
ind_ex.df[ind_ex.df==-1000]         <- NA
ind_ex2.df[ind_ex2.df==-1000]        <- NA

#Making new variable with rounded length categories 
# to nearest number with one decimal
ind2.df$length_rounded <- round(ind2.df$length, 1)

#Making a new variable with rounded environmental categories 
# to the nearest whole number
ind2.df$environment_rounded <- round(ind2.df$environment)

#Saving old environment format
ind2.df$environment2 <- ind2.df$environment

#Making max and minimum limits for the environment
valueofe.mx <- import.binary("./ValueOfE.bin") #Importing ValueOfE matrix from binary file
E.limits <- c(valueofe.mx[1], valueofe.mx[length(valueofe.mx)]) #Saving min and max

#Transforming environment into a scale that reflects ValueOfE for both static and stochastic dataset
ind2.df$environment       <- convert.scale(ind2.df$environment,       1, length(valueofe.mx), E.limits[1], E.limits[2]) 
ind.static.df$environment <- convert.scale(ind.static.df$environment, 1, length(valueofe.mx), E.limits[1], E.limits[2])
ind_ex.df$environment     <- convert.scale(ind_ex.df$environment,     1, length(valueofe.mx), E.limits[1], E.limits[2]) 
ind_ex2.df$environment     <- convert.scale(ind_ex2.df$environment,   1, length(valueofe.mx), E.limits[1], E.limits[2]) 

#Transforming environment_rounded into a scale that reflects ValueOfE for the stochastic dataset
ind2.df$environment_rounded       <- convert.scale(ind2.df$environment_rounded, 1, length(valueofe.mx), E.limits[1], E.limits[2]) 

#~ #Transforming environment into a scale from 0-1 for both static and stochastic dataset
#~ ind2.df$environment       <- (ind2.df$environment2-1)/(max(ind2.df$start_environment, na.rm=TRUE)-1)
#~ ind.static.df$environment <- (ind.static.df$environment-1)/(max(ind.static.df$environment)-1)

###########

#Calculating total mortality per year
ind2.df$M_total <- ind2.df$M_sizeindependent + ind2.df$M_size + ind2.df$M_size_foraging + ind2.df$M_size_O2 + ind2.df$M_size_O2_foraging 



#Function for simple plot with x, y and colour for a chosen length
plot.simple <- function(dataset.df,xaxis, xlabel, yaxis, ylabel, colgroup,plot.length){
  easy.plot <- ggplot(filter(dataset.df, length_rounded==plot.length), aes_string(x=xaxis, y=yaxis, colour=colgroup,fill=colgroup)) +
                geom_point() +
                theme_classic(base_size=16) +
                labs(x=xlabel,y=ylabel) +
                ggtitle(paste0("Length = ", plot.length," cm"))
}

###########



#Blank plot for filling up grids
blank.plot <- plot(0,type='n',axes=FALSE,ann=FALSE)



####### Figure A
#~ aplot.colours    <- brewer.pal(n = 5, name = "Set2")
aplot.colours    <- c("#F8C52C","#F5B683","#DFD0BB","#FFFAD9","#888D79")
aplot.line.size   <- 1.5
aplot.point.size  <- 2.25
aplot.base.size   <- 16
aplot.timestep.max  <- 100

#valueofe.mx <- import.binary("./ValueOfE.bin")
#proballe.mx <- import.binary("./ProbAllE.bin")

#environment.data <- data.frame(
#                      environment = seq(1,length(valueofe.mx),length.out=length(valueofe.mx)),
#                      valueofe    = valueofe.mx,
#                      proballe    = proballe.mx)
#environment.data



##ProbAllE
#figureA1.plot <- ggplot(environment.data, aes(x=valueofe, y=proballe, colour=proballe, fill=proballe)) +
#                  geom_line(size=aplot.line.size, colour="gray") + geom_point(size=aplot.point.size, colour="#888D79", shape=21) + 
#                  theme_classic(base_size=aplot.base.size) + geom_vline(xintercept=1,linetype="dotted") +
#                  scale_colour_gradient(high="gray",low="gray90", guide=FALSE) +
#                  scale_fill_gradient(  high="gray",low="gray90", guide=FALSE) +
#                  scale_x_continuous(breaks = seq(E.limits[1], E.limits[2], length.out= 5)) +
#                  labs(x="Environment", y="Probability")
#figureA1.plot



##ValueOfE
#figureA2.plot <- ggplot(environment.data, aes(x=environment,y=valueofe, colour=proballe, fill=proballe)) +
#                  geom_line(size=aplot.line.size) + geom_point(size=aplot.point.size, colour="#888D79", shape=21) + 
#                  theme_classic(base_size=aplot.base.size) + geom_vline(xintercept=0.5,linetype="dotted") +
#                  scale_colour_gradient(high="#888D79",low="gray90", guide=FALSE) +
#                  scale_fill_gradient(  high="#888D79",low="gray90", guide=FALSE) +
#                  labs(x="Environment", y="Food Value")
#figureA2.plot



##Environment over time

##Finding an average individual to plot that has a final timestep >= aplot.timestep.max
#ind <- filter(ind2.df, growth_speed=="Average", last_timestep >= aplot.timestep.max)$individual[1]

#figureA2.plot <- ggplot(filter(ind2.df, individual==ind), aes(x=environment, y=timestep, colour=environment, fill=environment)) +
#                  geom_path(size=aplot.line.size, colour="gray") + # geom_point(size=aplot.point.size-1, colour="#888D79", shape=21) + 
#                  theme_classic(base_size=aplot.base.size) + geom_vline(xintercept=1,linetype="dotted") +
##~                   scale_colour_gradient2(high="gray90", mid ="gray", low="gray90", midpoint=1, guide=FALSE) +
##~                   scale_fill_gradient2(  high="gray90", mid ="gray", low="gray90", midpoint=1, guide=FALSE) +
#                  scale_x_continuous(breaks = seq(E.limits[1], E.limits[2], length.out= 5)) +
#                  scale_y_continuous(limits=c(0,aplot.timestep.max)) +
#                  labs(x="Environment", y="Time")
#figureA2.plot



##ProbOfE
#probe.df      <- melt(import.binary("./ProbOfE.bin"))
#probe.df$Var1 <- convert.scale(probe.df$Var1, 1, length(valueofe.mx), valueofe.mx[1], valueofe.mx[length(valueofe.mx)]) 
#probe.df$Var2 <- convert.scale(probe.df$Var2, 1, length(valueofe.mx), valueofe.mx[1], valueofe.mx[length(valueofe.mx)]) 

##~ probe.df$Var1 <- (probe.df$Var1-1)/(max(probe.df$Var1)-1)
##~ probe.df$Var2 <- (probe.df$Var2-1)/(max(probe.df$Var2)-1)

#figureA3.plot <- ggplot(probe.df, aes(y=Var1, x=Var2, fill=value)) +
#                  geom_tile() + geom_text(aes(label=format(round(value, digits=2), nsmall = 2))) + 
#                  theme_classic(base_size=aplot.base.size) +
#                  scale_fill_gradient(high="gray",low="white", guide=FALSE) +
#                  scale_y_continuous(breaks = seq(E.limits[1], E.limits[2], length.out= 5)) +
#                  scale_x_continuous(breaks = seq(E.limits[1], E.limits[2], length.out= 5)) +
#                  labs(y="Next Environment",x="Environment")
#figureA3.plot


##Changing the order of the mortality factors
#mort.order <- c("M_size_O2_foraging","M_size_O2", "M_size_foraging", "M_size","M_sizeindependent")

##Making new dataset with mortality data stacked on top of each other
#figureA4.data <- ind2.df %>%
#            gather(., key=mortality_type, value=mortality_data, M_size:M_size_O2_foraging) %>%
#            mutate(mortality_type = factor(mortality_type, levels = mort.order)) %>% 
#            arrange(mortality_type) %>%
#            aggregate(mortality_data ~ length_rounded + mortality_type, ., mean)
            
#figureA4.plot <- ggplot(figureA4.data,aes(x=length_rounded, y=mortality_data, fill=mortality_type, label=mortality_type)) +
#                  geom_area(colour="black", size=.2, alpha=.4) +
#                  theme_classic(base_size=aplot.base.size) + 
#                  scale_x_continuous(limits = length.limits) +
#                  scale_fill_manual(
#                    labels = c(expression(paste("Size dependent ", "O"[2], " and foraging mortality")),
#                            expression(paste("Size dependent ", "O"[2], " mortality")),
#                            "Size dependent foraging mortality",
#                            "Size dependent mortality",
#                            "Size independent mortality"),
#                    values = aplot.colours) +
#                  labs(x="Length (cm)", y="Mortality (per year)") +
#                  theme(legend.position = c(0.71, 0.84),
##                         legend.background = element_rect(linetype="dotted", colour="grey"),
#                        legend.title=element_blank(),
#                        legend.text.align=0)                  
#figureA4.plot



##Plot with M_O2 ~ O2_used/O2max_THF

##Actual plotting
#figureA5.plot <- ggplot(ind2.df, aes(x=foraging_required, y=(M_size_foraging/M_size))) + 
#                  geom_line(colour=aplot.colours[3], size=aplot.line.size) +
#                  theme_classic(base_size=aplot.base.size) +
#                  labs(x=expression(paste("Foraging activity"," (multiples of SMR"["standard"],")")),y="Foraging mortality")
#figureA5.plot



##Plot with M_O2 ~ O2_used/O2max_THF
#figureA6.plot <- ggplot(ind2.df, aes(x=(O2_used/O2max_THF), y=(M_size_O2/M_size))) + 
#                  geom_line(colour=aplot.colours[2], size=aplot.line.size) +
#                  theme_classic(base_size=aplot.base.size) +
#                  labs(x=expression(paste(over(paste("O"[2]," used"), paste("O"[2]," max"["THF"])))), y=expression(paste("O"[2], " mortality")))
#figureA6.plot


#cowplot::plot_grid(figureA1.plot, figureA4.plot, figureA2.plot, figureA5.plot, blank.plot, figureA6.plot, labels=c("a)", "b)", "c)", "d)", "", "e)"), align="hv", nrow=3, ncol=2)

##Saving plot
#dev.copy(png, paste0(path,"figureA-test.png"), width=1000,height=1000)
#dev.off()

#rm(figureA2.data, figureA4.data)



###### Figure B: Individuals
bplot.colours       <- c("#8DA0CB", "#66C2A5", "#FC8D62")#c("#D55E00","#F0E442","#56B4E9")
bplot.fill          <- aplot.colours #colorRampPalette(rev(brewer.pal(1000, "Spectral")))
bplot.no.colour     <- "darkgrey"
bplot.line.size     <- c(1.6)
bplot.line.alpha    <- 0.5
bplot.point.size    <- 4 #1.6
bplot.base.size     <- aplot.base.size
bplot.timestep      <- c(52,58)
bplot.timesplit     <- 53
bplot.timecut       <- c(40,50)
bplot.strip.alpha   <- 0.1
bplot.strip.colour  <- c("#faf5f1")

#Plot order
plot.order <- c("environment", "GHF", "OXF", "THF", "length", "R_cat")#, "survival2")

y.plot.labels <-  c("environment" = "Food Availability",
                  "GHF"           = "GHF [ng/ml]",
                  "OXF"           = "OXF [pg/ml]",
                  "THF"           = "THF [ng/ml]",
                  "length"        = "Length [cm]",
                  "R_cat"         = "Reserves [%]")
#~                   "survival2"     = "Cumulative \n mortality risk"
#~                   )

colour.plot.labels <- c("Poor","Intermediate","Rich")

x.plot.labels <- c("1" = colour.plot.labels[1],
                   "2" = colour.plot.labels[2],
                   "3" = colour.plot.labels[3]
                  )

#Making dataset containg limits for the different classes
limit.df <- data.frame( environment   = c(min(valueofe.mx),max(valueofe.mx)),
                        GHF           = GHF.limits,
                        OXF           = OXF.limits,
                        THF           = THF.limits,
                        length        = length.limits,
                        R_cat         = c(0, 15),
                        survival2     = c(35, 41),
                        individual    = c(1,1)
                      )
limit.df

#Combining data in one coloumn for plot_grid and sorting individuals
figureB.data <-  filter(ind_ex.df, timestep >= bplot.timestep[1]-1 & timestep <= bplot.timestep[2]+1) %>%
                  mutate(survival2 = log(1/survival_next)) %>%
                  bind_rows(., limit.df) %>% 
                  gather(key=hormone.type, value=hormone.level, OXF, GHF, THF, environment, length, R_cat) %>%
                  mutate(hormone.type = factor(hormone.type, levels = plot.order)) %>% 
                  arrange(hormone.type) %>%
                  mutate(individual = factor(individual))
                  
#Actual plotting
figureB1.plot <- ggplot(figureB.data, aes(x=timestep, y=hormone.level, colour=individual, shape=individual)) +
                  geom_path(size=bplot.line.size[1], lty=3) + geom_point(size=bplot.point.size, stroke = 1.5) +
                  geom_path (data=filter(figureB.data, timestep>=bplot.timesplit+1), size=bplot.line.size[1]) +
                  geom_path (data=filter(figureB.data, timestep<=bplot.timesplit), size=bplot.line.size[1], colour=bplot.no.colour) +
                  geom_point(data=filter(figureB.data, timestep<=bplot.timesplit), size=bplot.point.size,   colour=bplot.no.colour, stroke = 1.5) +
                  theme_classic(base_size=bplot.base.size) +
                  geom_vline(xintercept=bplot.timesplit, colour=bplot.no.colour) +
                  scale_size_manual(values=bplot.line.size) +
                  scale_colour_manual(values=bplot.colours, labels=colour.plot.labels) +
                  scale_shape_manual(values=c(24,21,23), labels=colour.plot.labels) +
                  coord_cartesian(xlim=bplot.timestep) +
#~                   facet_grid(hormone.type~individual, scales="free_y", space="fixed", shrink=FALSE, switch="y", 
#~                             labeller=labeller(hormone.type=y.plot.labels, individual=x.plot.labels)) + 
                  facet_wrap(~hormone.type,ncol=1,scales="free_y", strip.position="left", labeller=as_labeller(y.plot.labels)) +
                  labs(x="Week", y="", colour="Individual", shape="Individual") +
                  #guides(shape=FALSE) + #, colour=FALSE) +
                  theme(legend.position="top")
figureB1.plot



#Mortality plots for the 3 individuals

#Plot order
mort.order <- c("M_size_O2_foraging","M_size_O2", "M_size_foraging", "M_size","M_sizeindependent")

#Plot labels
fill.plot.labels <- rev(c("Size-independent mortality",
                      "Size-dependent mortality", 
                      "Size-dependent foraging mortality",
                      "Size-dependent scope mortality",
                      "Size-dependent active-while- \n vulnerable mortality"))

#Changing the dataset and plots 
figureB2.plot <- filter(ind_ex.df, timestep >= bplot.timestep[1] & timestep <= bplot.timestep[2]) %>%
                  gather(key=mortality.type, value=mortality.data, M_size:M_size_O2_foraging) %>%
                  mutate(mortality.type = factor(mortality.type, levels = mort.order)) %>% 
                  arrange(mortality.type) %>%
                  mutate(individual = factor(individual, levels = c("3", "2", "1"))) %>% 
                  arrange(individual) %>%
                  
                  ggplot(., aes(x=timestep, y=mortality.data, fill=mortality.type)) +
                    geom_area(colour="black", size=.2, alpha=.75) +
                    theme_classic(base_size=bplot.base.size) +
                    geom_vline(xintercept=bplot.timesplit, colour=bplot.no.colour) +
                    scale_fill_manual(values=bplot.fill, labels=fill.plot.labels[1:4],
                                      breaks=mort.order[1:4]) +
                    coord_cartesian(xlim=bplot.timestep, ylim=c(0,1)) +
                    facet_wrap(~individual,ncol=1, strip.position="top", labeller=as_labeller(x.plot.labels)) +
                    labs(x="Week", y=expression(paste("Mortality [",year^-1,"]")), fill="", title="Mortality components") +
                    theme(legend.position="top") + 
                    guides(fill=guide_legend(nrow = 2, byrow = TRUE, reverse=TRUE))
figureB2.plot


cowplot::plot_grid(figureB1.plot, figureB2.plot, labels=c("(a)", "(b)"), align="hv", nrow=1, ncol=2)

#Saving plot
#~ dev.copy(png, paste0(path,"figureB.png"),width=1000,height=1000) # width=500, height=1000)
dev.copy(pdf, paste0(path,"figureB.pdf"),width=2,height=2) # width=500, height=1000)
dev.off()

rm(figureB.data)


#~ #Hormone levels in the different environments


#~ figureB2a.plot <- ind2.df %>%
#~                   gather(key=hormone.type, value=hormone.level, OXF, GHF, THF, environment,survival_next, length) %>%
#~                   filter(., hormone.type %in% c("GHF", "OXF", "THF"), environment_rounded %in% c(min(valueofe.mx), valueofe.mx[round((length(valueofe.mx)/2))], max(valueofe.mx)), growth_speed!="NA") %>%
#~                   aggregate(hormone.level~growth_speed+environment_rounded+hormone.type, ., mean) %>%
#~                   mutate(environment_rounded = round(.$environment_rounded,2)) %>%
#~                   #Plot
#~                   ggplot(., aes(y=hormone.level, x=growth_speed, colour=growth_speed)) +
#~                   geom_boxplot() + theme_classic(base_size=bplot.base.size) + guides(colour=FALSE) +
#~                   facet_grid(hormone.type~environment_rounded,scales="free", switch="y") +
#~                   scale_colour_manual(values=bplot.colours) +
#~                   labs(x="Growth Speed", y="")
#~ figureB2a.plot


#~ figureB2.data <- filter(figureB.data, hormone.type %in% c("environment","GHF", "OXF", "THF"))

#~ figureB2b.plot <- ggplot(filter(figureB2.data, individual!=ind.nr[2]),
#~                               aes(x=timestep, y=hormone.level, colour=growth_speed)) +
#~                   geom_path(alpha=bplot.line.alpha, size=bplot.line.size[1]) +
#~                   geom_path(data=filter(figureB2.data, individual==ind.nr[2]),
#~                               aes(x=timestep, y=hormone.level), colour=bplot.colours[2], size=bplot.line.size[2]) +
#~                   geom_point(data=figureB2.data, aes(x=timestep,y=hormone.level), size=(bplot.point.size-3)) +
#~                   theme_classic(base_size=bplot.base.size) +
#~                   geom_vline(xintercept=(bplot.timestep.max)) +
#~                   scale_size_manual(values=bplot.line.size) +
#~                   scale_colour_manual(values=bplot.colours) +
#~                   facet_wrap(~hormone.type,ncol=1,scales="free_y", strip.position="left", labeller=plot.labels) +
#~                   coord_cartesian(xlim=c((bplot.timecut[1]-0.5),bplot.timecut[2]+0.5), expand=FALSE) +
#~                   labs(x="Week", y="", colour="Individual") +
#~                   guides(size=FALSE, colour=FALSE) +
#~                   theme(panel.background = element_rect(fill = bplot.strip.colour))
#~ figureB2b.plot



#~ first_col  <- cowplot::plot_grid(figureB1.plot, labels=c("a)"))
#~ second_col <- cowplot::plot_grid(figureA1.plot, figureB2a.plot, labels=c("b)", "c)"), align="hv", nrow=2, rel_heights = c(1, 3))
#~ cowplot::plot_grid(first_col, second_col, nrow=1)


#~ #Saving plot
#~ dev.copy(png, paste0(path,"figureB-2-test.png"), width=1000,height=1000)
#~ dev.off()




###### Figure C
cplot.colours     <- colorRampPalette(rev(brewer.pal(1000, "Spectral"))) #"Spectral" "RdYlBu"
cplot.point.size  <- 3
cplot.point.shape <- 23
cplot.arrow.size  <- 0.1
cplot.arrow.length <- unit(0.25,"cm")
cplot.base.size   <- aplot.base.size

#Surplus before growth GHF per environment
figureC1.plot <- plot.simple(ind2.df,"environment", "", "surplus_before_growth", expression(paste(atop("Energy surplus", "[J "~week^-1~"]"))), "GHF", 20) + geom_hline(yintercept=0) + 
#~                 geom_line(data=ind.static.df, aes(x=environment,y=surplus_before_growth),colour="black") +
                geom_point(data=ind.static.df,aes(x=environment,y=surplus_before_growth, fill=GHF), colour="black",shape=cplot.point.shape, size=cplot.point.size) +
                scale_colour_gradientn(colours=cplot.colours(100), limits=GHF.limits) + scale_fill_gradientn(colours=cplot.colours(100), limits=GHF.limits) +
#~                 scale_x_continuous(breaks = seq(E.limits[1], E.limits[2], length.out= 5)) +
                theme(plot.title = element_blank(),
                    axis.title.y = element_text(size=(cplot.base.size-0.8))) +
                labs(colour="GHF \n[ng/ml]", fill="GHF \n[ng/ml]")
#~                 scale_colour_gradient(low="yellow", high="red") + scale_fill_gradient(low="yellow", high="red")



#Foraging mortality and OXF per environment

#Dividing points into groups based on M_size_foraging
figureC2.1.df <- ind2.df %>%
                mutate(axis_class = if_else(.$M_size_foraging <= 0.75, "B","A"))

#Dividing points into groups based on M_size_foraging
figureC2.2.df <- ind.static.df %>%
                mutate(axis_class = if_else(.$M_size_foraging <= 0.75, "B","A"))


figureC2.plot <- plot.simple(figureC2.1.df,"environment", "", "M_size_foraging", expression(paste(atop("Size-dependent Foraging Mortality", "["~year^-1~"]"))), "OXF", 20) +
#~                   geom_line(data=figureC2.2.df, aes(x=environment,y=M_size_foraging),colour="black") +
                  geom_point(data=figureC2.2.df,aes(x=environment,y=M_size_foraging, fill=OXF), colour="black", shape=cplot.point.shape, size=cplot.point.size) +
                  scale_colour_gradientn(colours=cplot.colours(100), limits=OXF.limits) + scale_fill_gradientn(colours=cplot.colours(100), limits=OXF.limits) +
#~                   scale_x_continuous(breaks = seq(E.limits[1], E.limits[2], length.out= 5)) +
                  facet_wrap(~axis_class,ncol=1, scales="free_y") +
                  theme(plot.title  = element_blank(),
                    axis.title.y    = element_text(size=(cplot.base.size-0.8)),
                    strip.background= element_blank(),
                    strip.text      = element_blank()) +
                  labs(colour="OXF \n[pg/ml]", fill="OXF \n[pg/ml]")


#O2 mortality and THF per environment
figureC3.plot <-  plot.simple(ind2.df,"environment", "Food Availability", "M_size_O2", expression(paste(atop("Size-dependent Scope Mortality", "["~year^-1~"]"))), "THF", 20) +
#~                   geom_line(data=ind.static.df, aes(x=environment,y=M_size_O2),colour="black") +
                  geom_point(data=ind.static.df,aes(x=environment,y=M_size_O2, fill=THF), colour="black",shape=cplot.point.shape, size=cplot.point.size) +
                  scale_colour_gradientn(colours=cplot.colours(100), limits=THF.limits) + scale_fill_gradientn(colours=cplot.colours(100), limits=THF.limits) +
#~                   scale_x_continuous(breaks = seq(E.limits[1], E.limits[2], length.out= 5)) +
                  theme(plot.title = element_blank(),
                    axis.title.y = element_text(size=(cplot.base.size-0.8))) +
                  labs(colour="THF \n[ng/ml]", fill="THF \n[ng/ml]")
#~                   scale_colour_gradientn(colours=c("darkblue","darkviolet","violet"),values=c(0,0.2,0.5,5), limits=THF.limits) +
#~                   scale_fill_gradientn(  colours=c("darkblue","darkviolet","violet"),values=c(0,0.2,0.5,5), limits=THF.limits)


cowplot::plot_grid(figureC1.plot, figureC2.plot, figureC3.plot, labels=c("a)", "b)", "c)"), align="v", nrow=3, ncol=1)

#Saving plot
dev.copy(png, paste0(path,"figureC.png"), width=500,height=1000)
dev.off()

rm(figureC2.1.df, figureC2.2.df)



###### Figure D
#~ dplot.colours       <- c("#f1a340", "#998ec3")
#~ dplot.line.size     <- 1.5
#~ dplot.point.size    <- 3
#~ dplot.arrow.length  <- unit(0.40,"cm")
#~ dplot.base.size     <- aplot.base.size

#~ ind.nr     <- as.integer(c(1421,256))
#~ length.nr  <- c(19,21,19.9,20.1)

#~ #Making a dataset containing the two individuals, finding first timestep for each individual and renaming them to A and B
#~ figureD.df <- union(
#~                 filter(ind2.df, individual==ind.nr[1], length>length.nr[1], length<length.nr[2], timestep>55),
#~                 filter(ind2.df, individual==ind.nr[2], length>length.nr[3], length<length.nr[4])
#~                   )

#~ figureD.df <- figureD.df %>%
#~               mutate(timestep_min = case_when(
#~                 .$individual == ind.nr[1] ~ min(filter(., individual == ind.nr[1])$timestep, na.rm=TRUE),
#~                 .$individual == ind.nr[2] ~ min(filter(., individual == ind.nr[2])$timestep, na.rm=TRUE)
#~               )) %>%
#~               mutate(timestep_max = case_when(
#~                 .$individual == ind.nr[1] ~ max(filter(., individual == ind.nr[1])$timestep, na.rm=TRUE),
#~                 .$individual == ind.nr[2] ~ max(filter(., individual == ind.nr[2])$timestep, na.rm=TRUE)
#~               )) %>%
#~               mutate(individual = case_when(
#~                 .$individual == ind.nr[1] ~ "A",
#~                 .$individual == ind.nr[2] ~ "B"
#~               ))


#~ #Actual plotting
#~ figureD1.plot <- ggplot(figureD.df, aes(y=surplus_before_growth, x=environment, colour=as.factor(individual), label=round(length,2), group=individual)) +
#~                   geom_point(data=filter(ind2.df,length>length.nr[1], length<length.nr[2]),aes(x=environment,y=surplus_before_growth),alpha=0.04,colour="lightgrey")+
#~                   geom_hline(yintercept=0, colour="darkgrey") +
#~                   geom_path(size=dplot.line.size, arrow=arrow(ends ="first", length=dplot.arrow.length)) +
#~                   geom_point(data=filter(figureD.df, timestep!=timestep_max), aes(y=surplus_before_growth, x=environment), size=dplot.point.size) +
                  geom_text(colour="black", nudge_y=1500, size=dplot.base.size-12) + #geom_point(size=dplot.point.size) +
                  scale_x_continuous(limits=E.limits, breaks = seq(E.limits[1], E.limits[2], length.out=5)) +
#~                   theme_classic(base_size=dplot.base.size) + 
#~                   scale_colour_manual(values=dplot.colours,
#~                     labels=c(
#~                       paste("A:",round(min(filter(figureD.df, individual=="A")$length, na.rm=TRUE), 2), "-", round(max(filter(figureD.df, individual=="A")$length, na.rm=TRUE), 2), "cm"),
#~                       paste("B:",round(min(filter(figureD.df, individual=="B")$length, na.rm=TRUE), 2), "-", round(max(filter(figureD.df, individual=="B")$length, na.rm=TRUE), 2), "cm")
#~                       )
#~                     )+
#~                   labs(x="", y="Energy surplus [J]", colour="") +
#~                   theme(legend.position="top") 

#~ figureD2.plot <- ggplot(figureD.df, aes(y=R_cat, x=environment, colour=as.factor(individual), label=round(length,2), group=individual)) +
#~                   geom_path(size=dplot.line.size, arrow=arrow(ends ="first", length=dplot.arrow.length)) +
#~                   geom_point(data=filter(figureD.df, timestep!=timestep_max), aes(y=R_cat, x=environment), size=dplot.point.size) +
                  geom_text(colour="black", nudge_y=0.4, size=dplot.base.size-12) + #geom_point(size=dplot.point.size) +
                  scale_x_continuous(limits=E.limits, breaks = seq(E.limits[1], E.limits[2], length.out= 5)) +
#~                   theme_classic(base_size=dplot.base.size) + 
#~                   scale_colour_manual(values=dplot.colours) + guides(colour=FALSE) +
#~                   labs(x="", y="Reserve Fullness [%]")
                
#~ figureD3.plot <- ggplot(figureD.df, aes(y=(growth*1000), x=environment, colour=as.factor(individual), label=round(length,2), group=individual)) +
#~                   geom_path(size=dplot.line.size, arrow=arrow(ends ="first", length=dplot.arrow.length)) +
#~                   geom_point(data=filter(figureD.df, timestep!=timestep_max), aes(y=(growth*1000), x=environment), size=dplot.point.size) +
                  geom_text(colour="black", nudge_y=0., size=dplot.base.size-12) + #geom_point(size=dplot.point.size) +
                  scale_x_continuous(limits=E.limits, breaks = seq(E.limits[1], E.limits[2], length.out= 5)) +
#~                   theme_classic(base_size=dplot.base.size) + 
#~                   scale_colour_manual(values=dplot.colours) + guides(colour=FALSE) +
#~                   labs(x="Food Availability", y=expression(paste("Growth [g ", "timestep"^-1,"]")))

#~ figureD4.plot <- ggplot(figureD.df, aes(y=length, x=environment, colour=as.factor(individual), label=round(length,2), group=individual)) +
#~                   geom_path(size=dplot.line.size, arrow=arrow(ends ="first", length=dplot.arrow.length)) +
#~                   geom_point(data=filter(figureD.df, timestep==timestep_min), aes(y=length, x=environment), size=dplot.point.size) +
                  geom_text(colour="black", nudge_y=0., size=dplot.base.size-12) + #geom_point(size=dplot.point.size) +
                  scale_x_continuous(limits=E.limits, breaks = seq(E.limits[1], E.limits[2], length.out= 5)) +
#~                   theme_classic(base_size=dplot.base.size) + 
#~                   scale_colour_manual(values=dplot.colours) + guides(colour=FALSE) +
#~                   labs(x="", y="Length [cm]")


#~ cowplot::plot_grid(figureD1.plot, figureD2.plot, figureD3.plot, labels=c("a)", "b)", "c)", "d)"), align="hv", ncol=1)

#~ #Saving plot
#~ dev.copy(png, paste0(path,"figureD.png"), width=500,height=1000)
#~ dev.off()

#~ rm(figureD.df)


#~ dplot.colours     <- colorRampPalette(rev(brewer.pal(1000, "Spectral"))) #RdYlBu

#~ figureD1.1.plot <- plot.simple(ind2.df,"environment","Environment","reserves","Reserves","reserve_diff",20) + 
#~                     scale_colour_gradientn(colours=dplot.colours(100)) + scale_fill_gradientn(colours=dplot.colours(100))
#~ figureD1.2.plot <- plot.simple(ind2.df,"environment","Environment","reserve_diff","Change in reserves","reserves",20) + scale_colour_gradientn(colours=cplot.colours(100)) +
#~                     scale_colour_gradientn(colours=dplot.colours(100), limits=c(0,NA)) + scale_fill_gradientn(colours=dplot.colours(100), limits=c(0,NA))

#~ grid.arrange(figureD1.1.plot, figureD1.2.plot)

#~ dev.copy(png, paste0(path,"figureD-test.png"), width=500,height=1000)
#~ dev.off()


#~ figureD2.plot <- plot.simple(ind2.df,"environment","Environment","surplus_before_growth","Energy surplus","reserve_diff",20) +
#~                   scale_colour_gradientn(colours=dplot.colours(100)) + scale_fill_gradientn(colours=dplot.colours(100))
#~ figureD2.plot




###### Figure E
eplot.point.size    <- 1
eplot.point.alpha   <- 1
#~ eplot.point.colour  <- "black"
eplot.line.colour   <- "black"
eplot.scale.colour  <- cplot.colours
eplot.base.size     <- aplot.base.size
eplot.line.size     <- bplot.line.size
eplot.last.timestep <-c(min(ind2.df$last_timestep, na.rm=TRUE), 200)

xtimesteps          <- 1  #Time steps between each "row" used in the figure E3 and E4
length.to.plot      <- 20 #Length to plot individuals in E3 and E4



#Plotte mean total mortality against last timestep
figureE1.plot <- ind2.df %>%
                  mutate(t_t_mortality = (1-survival_next)*100) %>% 
                  group_by(individual) %>%
                  mutate(environment = mean(environment, na.rm=TRUE))%>%
                  ungroup() %>% group_by(individual, last_timestep, environment) %>%
                  filter(timestep==last_timestep) %>%
                  
                  ggplot(., aes(x=20/last_timestep, y=t_t_mortality, colour=environment)) +
                  geom_point(alpha=eplot.point.alpha)+#, colour=eplot.point.colour) + 
                  scale_colour_gradientn(colours=eplot.scale.colour(100),
                    breaks=c(min(valueofe.mx), (min(valueofe.mx)+max(valueofe.mx))/2, max(valueofe.mx)),
                    limits=c(min(valueofe.mx), max(valueofe.mx))) +
                  theme_classic(base_size=eplot.base.size) +
                  #coord_cartesian(xlim=eplot.last.timestep) + 
                  geom_smooth(colour=eplot.line.colour, size=eplot.line.size, 
                    method="gam", formula = y ~ s(x, bs = "cs"), se=FALSE) +
                  labs( x = expression(paste("Mean Growth Rate [cm ", week^-1,"]")),
                        y = "Accumulated \n Mortality Risk [%]") +
                  guides(colour=FALSE)
                  
#Plotte mean total mortality against last timestep
figureE2.plot <- ind2.df %>%
                  mutate(total.mortality = M_sizeindependent + M_size +
                          M_size_foraging + M_size_O2 + M_size_O2_foraging) %>% 
                  group_by(individual) %>%
                  mutate(environment=mean(environment, na.rm=TRUE)) %>%
                  group_by(individual, last_timestep) %>%
                  ungroup() %>% group_by(individual, last_timestep, environment) %>%
                  summarize(m_t_mortality = mean(total.mortality, na.rm=TRUE)) %>% 
                  mutate(m_t_mortality = (1-exp(-m_t_mortality))*100) %>%
                  
                  ggplot(., aes(x=20/last_timestep, y=m_t_mortality, colour=environment)) +
                  geom_point(alpha=eplot.point.alpha) + 
                  geom_smooth(colour=eplot.line.colour, size=eplot.line.size, 
                    method="gam", formula = y ~ s(x, bs = "cs"), se=FALSE) +
                  scale_colour_gradientn(colours=eplot.scale.colour(100),
                    breaks=c(min(valueofe.mx), (min(valueofe.mx)+max(valueofe.mx))/2, max(valueofe.mx)),
                    limits=c(min(valueofe.mx), max(valueofe.mx))) +
#~                  geom_smooth(aes(group=environment_rounded), colour="b
                  theme_classic(base_size=eplot.base.size) +
                  #coord_cartesian(xlim=eplot.last.timestep) + #, ylim=c(0.5,0.9)) +
                  labs( x = "", expression(paste("Mean Growth Rate [cm ", week^-1,"]")),
                        y = expression(paste(atop("Mean Short Term", paste("Mortality Risk [% ", year^-1,"]")))),
                        colour=paste("Mean Food Availability", paste(replicate(8, " "), collapse = ""))) +
                  guides(colour=FALSE)
#~                   theme(legend.position="top")
figureE2.plot

#Function of summarising up to nearest xtimesteps
roundUp <- function(x) ceiling(x/xtimesteps)*xtimesteps

#Dataset with all columns summarised except length and last_timestep, and environment (where we use mean)
sum10.df <- ind2.df %>%
          mutate(rtime = roundUp(timestep)) %>%
          mutate(total.mortality = M_sizeindependent + M_size +
                  M_size_foraging + M_size_O2 + M_size_O2_foraging) %>% 
          mutate(last_timestep = if_else(timestep==rtime, last_timestep, 0),
                 length        = if_else(timestep==rtime, length,        0)) %>%
          group_by(individual, rtime) %>%
          mutate(environment = mean(environment))%>%
          mutate(m.total.mortality = mean(total.mortality))%>%
          mutate(environment        = if_else(timestep==rtime, environment,       0),
                  m.total.mortality = if_else(timestep==rtime, m.total.mortality, 0)) %>%
          select(-growth_speed, -timestep) %>%
          summarise_all(sum) %>%
          mutate( last_timestep     = ifelse((rtime==max(rtime) & last_timestep    ==0), NA, last_timestep),
                  length            = ifelse((rtime==max(rtime) & length           ==0), NA, length),
                  environment       = ifelse((rtime==max(rtime) & environment      ==0), NA, environment),
                  m.total.mortality = ifelse((rtime==max(rtime) & m.total.mortality==0), NA, m.total.mortality)) %>%
          arrange(individual, by_group=rtime) %>%
          #Growth in the xtimesteps time steps
          mutate(growth.rate = (lead(length)-length)/xtimesteps) %>% #Lead one row underneath
          #Selecting out rtime when where where nearest to 20 cm/length.to.plot
          ungroup() %>% group_by(individual) %>%
          filter(length==nth(length, which.min(abs(length-length.to.plot)))) %>%
          distinct(length, .keep_all = TRUE) #Removing duplicated rows and keep the first one 
          
#~ sum10.df$environment_rounded <- round(sum10.df$environment/0.3)*0.3

#Plotte mean total mortality against last time step
figureE3.plot <- ggplot(sum10.df, aes(x=growth.rate, y=(1-exp(-total.mortality))*100, colour=environment)) +
                  geom_point(alpha=eplot.point.alpha) + #, colour=eplot.point.colour) + 
                  scale_colour_gradientn(colours=eplot.scale.colour(100),
                    breaks=c(min(valueofe.mx), (min(valueofe.mx)+max(valueofe.mx))/2, max(valueofe.mx)),
                    limits=c(min(valueofe.mx), max(valueofe.mx))) +
#~                   geom_smooth(aes(group=environment_rounded), colour="black", 
#~                     method="gam", formula = y ~ s(x, bs = "cs"), se=FALSE, size=eplot.line.size+1) +
#~                   geom_smooth(aes(colour=environment_rounded, group=environment_rounded),
#~                     method="gam", formula = y ~ s(x, bs = "cs"), se=FALSE, size=eplot.line.size) +
                  theme_classic(base_size=eplot.base.size) +
#~                   coord_cartesian(y=c(0,1)) +
                  labs( x="", 
                        y="Total \n Mortality Risk for one week [%]",
                        colour="Environment")+
                  guides(colour=FALSE)

#Plotting mean mortality against growth per time step
figureE4.plot <- ggplot(sum10.df, aes(x=growth.rate, y=(1-exp(-m.total.mortality))*100,
                                      colour=environment)) +
                geom_point(alpha=eplot.point.alpha) + #, colour=eplot.point.colour) +
                scale_colour_gradientn(colours=eplot.scale.colour(100),
                  breaks=c(min(valueofe.mx), (min(valueofe.mx)+max(valueofe.mx))/2, max(valueofe.mx)),
                  limits=c(min(valueofe.mx), max(valueofe.mx))) +
#~                  geom_smooth(aes(group=environment_rounded), colour="black", 
#~                    method="gam", formula = y ~ s(x, bs = "cs"), se=FALSE, size=eplot.line.size+1) +
#~                  geom_smooth(aes(colour=environment_rounded, group=environment_rounded),
#~                    method="gam", formula = y ~ s(x, bs = "cs"), se=FALSE, size=eplot.line.size) +
                theme_classic(base_size=eplot.base.size) +
                coord_cartesian(y=c(25,60)) +
                labs( x=expression(paste("Short Term Growth Rate [cm ", week^-1,"]")),
                      y = expression(paste(atop("Short Term", paste("Mortality Risk [% ", year^-1,"]")))),
                      colour=paste("Mean Food Availability", paste(replicate(8, " "), collapse = ""))) +
                theme(legend.position="bottom")
figureE4.plot

#~ legend <- cowplot::get_legend(figureE4.plot)

figureE <- cowplot::plot_grid(figureE4.plot, figureE2.plot, figureE1.plot,
            labels=c("a)", "b)", "c)"), align="hv", ncol=1, nrow=3)
figureE
#~ cowplot::plot_grid(legend,figureE, ncol=1, nrow=2, rel_heights = c(0.1, 1))

#Saving plot
dev.copy(png, paste0(path,"figureE.png"), width=500,height=1000)
dev.off()

rm(sum10.df)




###### Figure F: Individuals
## NB!!! This figure needs Figure B to be plotted first
fplot.colours       <- c(rep("gray40",4),"#FC8D62") #c("#FC8D62", "#66C2A5", "#8DA0CB")
fplot.line.size     <- c(rep(0.8,4),1.2)
fplot.point.size    <- 3
fplot.base.size     <- aplot.base.size
fplot.strip.alpha   <- 0.2
fplot.point.alpha   <- 0.1
fplot.timestep.max  <- 150
fplot.label.hjust   <- -0.3
fplot.label.colour  <- "black"
fplot.quant.individuals   <- c(0.1, 0.3, 0.5, 0.7, 0.9) #The quantiles used for the 5 individuals based on when they reached 30cm

#Ranges used in background in plots (stat_summary(geom="ribbon"...)
fplot.lower.quant1 = function(z) {quantile(z,0   )}
fplot.upper.quant1 = function(z) {quantile(z,1   )}
fplot.lower.quant2 = function(z) {quantile(z,0.25)}
fplot.upper.quant2 = function(z) {quantile(z,0.75)}

#Function for finding individual matching to the quantiles
find.quant.ind <- function(dataset, quant){
  #Filtering out only the last time step from the dataset and arranging individual
  # by when they reached 30 cm
  quant.df  <- filter(dataset, timestep==last_timestep) %>%
                  arrange(timestep)
                  
  #Finding the id of the individual matching the quantile
  quant.ind <- quant.df$individual[floor((length(quant.df$individual)+1)*quant)]
  
  #Returning the id of the individual found
  return(quant.ind)
  }
  
#Finding the IDs of the individuals used in the plots
ind.nr <- c(find.quant.ind(ind2.df, fplot.quant.individuals[1]),
            find.quant.ind(ind2.df, fplot.quant.individuals[2]),
            find.quant.ind(ind2.df, fplot.quant.individuals[3]),
            find.quant.ind(ind2.df, fplot.quant.individuals[4]),
            find.quant.ind(ind2.df, fplot.quant.individuals[5]))
ind.nr 




#Plot of the environment for the 0.5 quantile individual
figureF1.plot <- filter(ind2.df, individual == find.quant.ind(ind2.df,0.5)) %>%
                 
                  ggplot(., aes(x=timestep, y=environment)) +
#~                   stat_summary(data=ind2.df, aes(x=timestep, y=environment), geom="ribbon", fun.ymin = fplot.lower.quant1, fun.ymax = fplot.upper.quant1, size=0, colour="white", alpha = 0.1) +
#~                   stat_summary(data=ind2.df, aes(x=timestep, y=environment), geom="ribbon", fun.ymin = fplot.lower.quant2, fun.ymax = fplot.upper.quant2, size=0, colour="white", alpha = 0.1) +
                  geom_point(data=ind2.df,alpha=fplot.point.alpha, colour="lightgrey")+#, width=0.5) +
                  geom_path(size=fplot.line.size[5], colour=fplot.colours[5]) +
#~                   geom_point(data=filter(ind2.df, individual == find.quant.ind(ind2.df,0.5)),
#~                     aes(fill=environment), size=fplot.point.size, shape=21)+
#~                   scale_fill_gradientn(colours=eplot.scale.colour(100),
#~                     breaks=c(min(valueofe.mx), (min(valueofe.mx)+max(valueofe.mx))/2, max(valueofe.mx)),
#~                     limits=c(min(valueofe.mx), max(valueofe.mx))) +
                  geom_point(
                    data=filter(ind2.df, individual == find.quant.ind(ind2.df,0.5),
                    timestep==last_timestep), size=fplot.point.size, colour=fplot.colours[5])+ 
                  geom_dl(aes(label = "0.5"), method = list("last.points", hjust=fplot.label.hjust), colour=fplot.label.colour) +
                  theme_classic(base_size=fplot.base.size) +
                  coord_cartesian(xlim=c(0,fplot.timestep.max),
                                  ylim=c(min(valueofe.mx), max(valueofe.mx))) +
                  labs(x="", y="Food Availability") + guides(colour=FALSE, fill=FALSE)
#~ figureF1.plot


#Making a dataset containing only the 5 quantile individuals
figureF.data <- filter(ind2.df, individual %in% ind.nr) %>%
                  mutate(survival2 = (1-survival_next)*100) %>% #Calculating cumulative mortality
                  #Giving all a new variable for the quantile they are, and turning this into a factor
                  # also changing the order of the individuals based on this new factor
                  # so that the 0.5 quantile individual is plotted last, and therefore on top
                  # of the others
                  mutate(quant = factor(case_when( 
                        .$individual == ind.nr[1] ~ fplot.quant.individuals[1],
                        .$individual == ind.nr[2] ~ fplot.quant.individuals[2],
                        .$individual == ind.nr[3] ~ fplot.quant.individuals[3],
                        .$individual == ind.nr[4] ~ fplot.quant.individuals[4],
                        .$individual == ind.nr[5] ~ fplot.quant.individuals[5]               
                        ),
                      levels=c( fplot.quant.individuals[1], fplot.quant.individuals[2],
                                fplot.quant.individuals[4], fplot.quant.individuals[5],
                                fplot.quant.individuals[3]))) %>%
                  #Change so last length is always exacly 30cm
                  mutate(length = ifelse(timestep==(last_timestep+1), 30, length))


#Plot with all the quantile individuals with length over time
figureF2.plot <- ggplot(figureF.data, aes(x=timestep, y=length, colour=quant, size=quant)) +
#~                   stat_summary(data=ind2.df, aes(x=timestep, y=length), geom="ribbon", fun.ymin = fplot.lower.quant1, fun.ymax = fplot.upper.quant1, size=0, colour="white", alpha = 0.1) +
#~                   stat_summary(data=ind2.df, aes(x=timestep, y=length), geom="ribbon", fun.ymin = fplot.lower.quant2, fun.ymax = fplot.upper.quant2, size=0, colour="white", alpha = 0.1) +
                  geom_point(data=ind2.df, aes(x=timestep, y=length), alpha=fplot.point.alpha, colour="lightgrey", size=1) + #, width=0.5) +
                  geom_path() +
                  geom_point(
                    data=filter(figureF.data, timestep==last_timestep+1),
                    size=fplot.point.size)+ 
                  geom_dl(aes(label = quant), method = list("last.points", hjust=fplot.label.hjust), colour=fplot.label.colour) +
                  theme_classic(base_size=fplot.base.size) +
                  scale_size_manual(values=fplot.line.size) +
                  scale_colour_manual(values=fplot.colours) +
                  coord_cartesian(xlim=c(0,fplot.timestep.max),
                                  ylim=length.limits) +
                  labs(x="", y="Length [cm]") + guides(size=FALSE, colour=FALSE)
#~ figureF2.plot


#Plot with all the quantile individuals with length over time
figureF3.plot <- ggplot(figureF.data,
#filter(figureF.data, individual %in% c(ind.nr[1], ind.nr[3], ind.nr[5])),
                        aes(x=timestep, y=survival2, colour=quant, size=quant)) +
#~                   stat_summary(data=ind2.df, aes(x=timestep, y=((1-survival_next)*100)), geom="ribbon", fun.ymin = fplot.lower.quant1, fun.ymax = fplot.upper.quant1, size=0, colour="white", alpha = 0.1) +
#~                   stat_summary(data=ind2.df, aes(x=timestep, y=((1-survival_next)*100)), geom="ribbon", fun.ymin = fplot.lower.quant2, fun.ymax = fplot.upper.quant2, size=0, colour="white", alpha = 0.1) +
                  geom_point(data=ind2.df,aes(x=timestep, y=((1-survival_next)*100)), alpha=fplot.point.alpha, colour="lightgrey", size=1) +#, width=0.5) +
                  geom_path() +
                  geom_point(data=filter(figureF.data, #individual %in% c(ind.nr[1], ind.nr[3], ind.nr[5]),
                    timestep==last_timestep+1),
                    size=fplot.point.size)+ 
                  geom_dl(aes(label = quant), method = list("last.points", hjust=fplot.label.hjust), colour=fplot.label.colour) +
                  theme_classic(base_size=fplot.base.size) +
                  scale_size_manual(values=fplot.line.size) +
                  scale_colour_manual(values=fplot.colours) +
#~                   scale_size_manual(values=c(fplot.line.size[1], fplot.line.size[3], fplot.line.size[5])) +
#~                   scale_colour_manual(values=c(fplot.colours[1], fplot.colours[3], fplot.colours[5])) +
                  coord_cartesian(xlim=c(0,fplot.timestep.max)) +
                  labs(x="Week", y=expression(paste("Accumulated Mortality Risk [%]"))) + guides(size=FALSE, colour=FALSE)
#~ figureF3.plot


#Combining plots
x11()
cowplot::plot_grid(figureF1.plot, figureF2.plot, figureF3.plot, labels=c("(a)", "(b)", "(c)"), align="hv", ncol=1)

#Saving plot
dev.copy(png, paste0(path,"figureF.png"), width=(8500/2),height=8500, res=600)
dev.off()

rm(figureF.data)





###### Figure G
#Individuals through environments
gplot.colours       <- c("#f1a340", "#998ec3")
gplot.base.size     <- aplot.base.size   
gplot.line.size     <- 1.5
gplot.point.size    <- 3
gplot.arrow.length  <- unit(0.40,"cm")
gplot.base.size     <- aplot.base.size
gplot.base.timestep <- 56

#Plot
ind_ex2.df$individual <- recode(ind_ex2.df$individual, "1"="A", "2"="B")

figureG.df <- union(
                filter(ind_ex2.df, individual=="A", timestep >= gplot.base.timestep, timestep < gplot.base.timestep+6),
                filter(ind_ex2.df, individual=="B", timestep >= gplot.base.timestep, timestep < gplot.base.timestep+10),
                  )

length.nr <- c(min(filter(figureG.df, individual=="A")$length, na.rm=TRUE), max(filter(figureG.df, individual=="A")$length))

figureG0.plot <- ggplot(figureG.df, aes(y=environment, x=timestep, colour=individual, label=round(length,2), group=individual)) +
                  geom_point(data=union(filter(figureG.df, individual=="A", timestep<gplot.base.timestep+5),filter(figureG.df, individual=="B", timestep<gplot.base.timestep+9)), 
                                aes(y=environment, x=timestep), size=gplot.point.size) +
                  geom_path(size=gplot.line.size, arrow=arrow(ends ="first", length=gplot.arrow.length)) +
                  theme_classic(base_size=gplot.base.size) + 
                  scale_x_continuous(breaks=c(min(figureG.df$timestep):max(figureG.df$timestep))) +
                  scale_colour_manual(values=gplot.colours,
                    labels=c(
                        paste("A:",round(min(filter(figureG.df, individual=="A")$length, na.rm=TRUE), 2), "-", round(max(filter(figureG.df, individual=="A")$length, na.rm=TRUE), 2), "cm"),
                        paste("B:",round(min(filter(figureG.df, individual=="B")$length, na.rm=TRUE), 2), "-", round(max(filter(figureG.df, individual=="B")$length, na.rm=TRUE), 2), "cm")
                        )
                    )+
                  theme(legend.position="bottom") +
                  labs(y="Food Availability", x="Week", colour="Individual")

figureG1.plot <- ggplot(figureG.df, aes(y=surplus_before_growth, x=environment, colour=individual, label=round(length,2), group=individual)) +
                  geom_point(data=filter(ind2.df,length>length.nr[1], length<length.nr[2]),aes(x=environment,y=surplus_before_growth),alpha=0.04,colour="lightgrey")+
                  geom_hline(yintercept=0, colour="darkgrey") +
                  geom_point(data=union(filter(figureG.df, individual=="A", timestep<gplot.base.timestep+5),filter(figureG.df, individual=="B", timestep<gplot.base.timestep+9)), 
                                aes(y=surplus_before_growth, x=environment), size=gplot.point.size) +
                  geom_path(size=gplot.line.size, arrow=arrow(ends ="first", length=gplot.arrow.length)) +
#~                   geom_text(colour="black", nudge_y=1500, size=gplot.base.size-12) + #geom_point(size=gplot.point.size) +
                  theme_classic(base_size=gplot.base.size) + 
                  scale_colour_manual(values=gplot.colours) + guides(colour=FALSE) +
                  labs(x="", y=expression(paste(atop("Energy surplus", "[J "~week^-1~"]"))), colour="")

figureG2.plot <- ggplot(figureG.df, aes(y=R_cat, x=environment, colour=individual, label=round(length,2), group=individual)) +
                  geom_point(data=union(filter(figureG.df, individual=="A", timestep<gplot.base.timestep+5),filter(figureG.df, individual=="B", timestep<gplot.base.timestep+9)), 
                                aes(y=R_cat, x=environment), size=gplot.point.size) +
                  geom_path(size=gplot.line.size, arrow=arrow(ends ="first", length=gplot.arrow.length)) +
#~                   geom_text(colour="black", nudge_y=0.4, size=dplot.base.size-12) + #geom_point(size=dplot.point.size) +
                  theme_classic(base_size=gplot.base.size) + 
                  scale_colour_manual(values=gplot.colours) + guides(colour=FALSE) +
                  labs(x="", y="Reserves \n [%]")

figureG3.plot <- ggplot(figureG.df, aes(y=(growth*1000), x=environment, colour=individual, label=round(length,2), group=individual)) +
                  geom_point(data=union(filter(figureG.df, individual=="A", timestep<gplot.base.timestep+5),filter(figureG.df, individual=="B", timestep<gplot.base.timestep+9)),
                                aes(y=(growth*1000), x=environment), size=gplot.point.size) +
                  geom_path(size=gplot.line.size, arrow=arrow(ends ="first", length=gplot.arrow.length)) +
#~                   geom_text(colour="black", nudge_y=0., size=dplot.base.size-12) + #geom_point(size=dplot.point.size) +
                  theme_classic(base_size=gplot.base.size) + 
                  scale_colour_manual(values=gplot.colours) + guides(colour=FALSE) +
                  labs(x="Food Availability", y=expression(paste(atop("Growth", "[g "~week^-1~"]"))))

figureG4.plot <- ggplot(figureG.df, aes(y=length, x=environment, colour=individual, label=round(length,2), group=individual)) +
                  geom_point(data=union(filter(figureG.df, individual=="A", timestep<gplot.base.timestep+5),filter(figureG.df, individual=="B", timestep<gplot.base.timestep+9)), 
                                aes(y=length, x=environment), size=gplot.point.size) +
                  geom_path(size=gplot.line.size, arrow=arrow(ends ="first", length=gplot.arrow.length)) +
#~                   geom_text(colour="black", nudge_y=0., size=dplot.base.size-12) + #geom_point(size=dplot.point.size) +
                  theme_classic(base_size=gplot.base.size) + 
                  scale_colour_manual(values=gplot.colours) + guides(colour=FALSE) +
                  labs(x="", y="Length (cm)")


cowplot::plot_grid(figureG0.plot, figureG1.plot, figureG2.plot, figureG3.plot, labels=c("(a)", "(b)", "(c)", "(d)"), align="hv", ncol=1)

#Saving plot
dev.copy(png, paste0(path,"figureF.png"), width=(8500/2),height=8500, res=600)
dev.off()

rm(figureG.df)






#------------- Variations of Figure C ----------------------------------------------------------------------

###### Figure C
cplot.colours     <- colorRampPalette(rev(brewer.pal(1000, "Spectral"))) #RdYlBu
cplot.point.size  <- 3
cplot.point.shape <- 21
cplot.arrow.size  <- 0.1
cplot.arrow.length <- unit(0.25,"cm")
cplot.base.size   <- aplot.base.size

#Surplus before growth GHF per environment
figureC1.plot <- plot.simple(ind2.df,"environment", "", "surplus_before_growth","Energy surplus", "GHF", 20) + geom_hline(yintercept=0) + 
                geom_line(data=ind.static.df, aes(x=environment,y=surplus_before_growth),colour="black") +
                geom_point(data=ind.static.df,aes(x=environment,y=surplus_before_growth, fill=GHF), colour="black",shape=cplot.point.shape, size=cplot.point.size) +
                scale_colour_gradientn(colours=cplot.colours(100)) + scale_fill_gradientn(colours=cplot.colours(100)) +
                theme(axis.title.y = element_text(size=(cplot.base.size-0.8)))
#~                 scale_colour_gradient(low="yellow", high="red") + scale_fill_gradient(low="yellow", high="red")



#Foraging mortality and OXF per environment
figureC2.plot <- plot.simple(ind2.df,"environment", "", "M_size_foraging","Size-dependent foraging mortality per year", "OXF", 20) +
                  geom_line(data=ind.static.df, aes(x=environment,y=M_size_foraging),colour="black") +
                  geom_point(data=ind.static.df,aes(x=environment,y=M_size_foraging, fill=OXF), colour="black", shape=cplot.point.shape, size=cplot.point.size) +
                  scale_colour_gradientn(colours=cplot.colours(100)) + scale_fill_gradientn(colours=cplot.colours(100)) +
                  coord_cartesian(ylim=c(0,max(filter(ind2.df, length_rounded==20)$M_size_foraging, na.rm=TRUE)+0.5)) +
#~                   coord_cartesian(ylim=c(0,3.5)) +
                  theme(plot.title  = element_blank(),
                    axis.title.y    = element_text(size=(cplot.base.size-0.8)),
                    axis.line.y     = element_line(arrow = arrow(length = cplot.arrow.length)))
#~                   scale_colour_gradient(low="purple", high="red") + scale_fill_gradient(low="purple", high="red")


        #Point for max in the static.model as this falls outside the plotting area for figureC2.plot
        point.C2.df <- filter(ind.static.df, environment==0) #Making new dataset
        point.C2.df$M_size_foraging2 <- point.C2.df$M_size_foraging #Saving old M_size_foraging value
        point.C2.df$M_size_foraging  <- max(filter(ind2.df, length_rounded==20)$M_size_foraging, na.rm=TRUE)+0.2 #Replacing with max for length_rounded=20 from stochastic simulation
#~         point.C2.df$M_size_foraging  <- 3 #Replacing with max for length_rounded=20 from stochastic simulation


  #Adding point to plot
  figureC2.plot <- figureC2.plot +
#~                     geom_segment(aes(x = point.C2.df$environment, xend = point.C2.df$environment,
#~                       y = point.C2.df$M_size_foraging, yend = (point.C2.df$M_size_foraging + 0.38)),
#~                       arrow = arrow(length = unit(0.2,"cm")), colour="black", size=cplot.arrow.size) +
                    geom_point(data=point.C2.df,aes(x=environment,y=M_size_foraging, fill=OXF), colour="black", shape=cplot.point.shape, size=cplot.point.size) +
#~                     annotate("text", x=point.C2.df$environment, y=(point.C2.df$M_size_foraging + 0.51),
#~                       label=paste0(round(point.C2.df$M_size_foraging2,2)))
                    annotate("text", x=point.C2.df$environment, y=(point.C2.df$M_size_foraging + 0.25),
                      label=paste0(round(point.C2.df$M_size_foraging2,2)))
figureC2.plot



#O2 mortality and THF per environment
figureC3.plot <-  plot.simple(ind2.df,"environment", "Environment", "M_size_O2","Size-dependent oxygen mortality (per year)", "THF", 20) +
                  geom_line(data=ind.static.df, aes(x=environment,y=M_size_O2),colour="black") +
                  geom_point(data=ind.static.df,aes(x=environment,y=M_size_O2, fill=THF), colour="black",shape=cplot.point.shape, size=cplot.point.size) +
                  scale_colour_gradientn(colours=cplot.colours(100), limits=THF.limits) + scale_fill_gradientn(colours=cplot.colours(100), limits=THF.limits) +
                  expand_limits(y=0) +
                  labs(title=element_blank(), x="Environment", y=expression(paste("Size dependent ", "O"[2], " mortality per year"))) +
                  theme(plot.title = element_blank(),
                    axis.title.y = element_text(size=(cplot.base.size-0.8)))
#~                   scale_colour_gradientn(colours=c("darkblue","darkviolet","violet"),values=c(0,0.2,0.5,5), limits=THF.limits) +
#~                   scale_fill_gradientn(  colours=c("darkblue","darkviolet","violet"),values=c(0,0.2,0.5,5), limits=THF.limits)


cowplot::plot_grid(figureC1.plot, figureC2.plot, figureC3.plot, labels=c("a)", "b)", "c)"), align="hv", nrow=3, ncol=1)


#Saving plot
dev.copy(png, paste0(path,"figureF.png"), width=(8500/2),height=8500, res=600)
dev.off()




#------------------------------------------------------------------------------------------

minE20 <- min(filter(ind2.df,length_rounded==20)$environment,na.rm=TRUE)
filter(ind2.df,length_rounded==20,environment>0.999)$individual[1:10]

ind.nr    <- 1421
length.nr <- c(19,21)

my.colour.test <- colorRampPalette(rev(brewer.pal(100, "RdYlBu")))

test1 <- ind2.df %>%
        filter(.,individual==ind.nr, length>length.nr[1], length<length.nr[2]) %>%
        ggplot(aes(y=surplus_before_growth, x=environment, colour=GHF, label=timestep)) +
          geom_point(data=filter(ind2.df,length>length.nr[1], length<length.nr[2]),aes(x=environment,y=surplus_before_growth),alpha=0.04,colour="grey")+
          geom_path(size=1,colour="black")+ geom_point(size=3) + geom_text(colour="black", nudge_y=700) +
          scale_x_continuous(limits=c(0,1)) +
          theme_classic(base_size=16) +
          scale_colour_gradientn(colours=my.colour.test(100), limits=c(0,90))

test2 <- ind2.df %>%
        filter(.,individual==ind.nr, length>length.nr[1], length<length.nr[2]) %>%
        ggplot(aes(y=R_cat, x=environment, colour=GHF, label=timestep)) +
          geom_path(size=1,colour="black")+ geom_point(size=3) + geom_text(colour="black", nudge_y=0.5) +
          scale_y_continuous(limits=c(0,25)) + scale_x_continuous(limits=c(0,1)) +
          theme_classic(base_size=16) +
          scale_colour_gradientn(colours=my.colour.test(100), limits=c(0,90))
          
grid.arrange(test1,test2,ncol=1,nrow=2)

dev.copy(png, paste0(path,"figureF.png"), width=4250,height=8500, res=600)
dev.off()




#Saving the last timestep of the individuals life
lasttime.var <- max(filter(ind2.df, individual==ind.nr)$last_timestep)

#Saving the survival (survival_next) at the end of the individuals life
survival.var <- filter(ind2.df, individual == ind.nr, timestep==(lasttime.var+1))$survival_next
    
#Saving the growth speed group of the individual
growth_speed.var <- filter(ind2.df, individual==ind.nr, timestep==1)$growth_speed

#Saving the end status of the individual as a variable
if(filter(ind2.df, individual==ind.nr)$end_status[1] == 1){end.status.var="Dead"} 
if(filter(ind2.df, individual==ind.nr)$end_status[1] == 2){end.status.var="To slow"} 
if(filter(ind2.df, individual==ind.nr)$end_status[1] == 3){end.status.var="Reached max length"} 


#Ordering levels used in plots
#Changing the order of the growth speed groups, for better plots
plot.order <- c("environment", "GHF", "OXF", "THF", "M_total", "survival_next")



#Plot with all hormones for one individual vs length
test1.plot <- filter(ind2.df, individual %in% ind.nr) %>% 
              gather(key=hormone.type, value=hormone.level, OXF, GHF, THF, environment,M_total,survival_next) %>%
              aggregate(hormone.level ~ length + hormone.type + growth_speed, ., mean) %>%
              mutate(hormone.type = factor(hormone.type, levels = plot.order)) %>% 
              arrange(hormone.type) %>% 
              mutate(hormone.level = case_when(
                hormone.type == "THF" ~ (hormone.level/THF.limits[2])*100,
                hormone.type == "GHF" ~ (hormone.level/GHF.limits[2])*100,
                hormone.type == "OXF" ~ (hormone.level/OXF.limits[2])*100,
                hormone.type == "environment" ~ ((hormone.level-1)/(max(ind2.df$start_environment, na.rm=TRUE)-1)*100),
                hormone.type == "M_total" ~ (hormone.level/2*100),                
                hormone.type == "survival_next" ~ (hormone.level*100)
                )) %>%
              ggplot(aes(x=length, y=hormone.level, colour=growth_speed)) +
              geom_path() + theme_classic(base_size=16) +
              facet_wrap(~hormone.type,ncol=1,nrow=6,scales="free_y") +
              coord_cartesian(ylim=c(0, 100)) +
#~               scale_y_continuous(limits = c(0, 100)) + 
              labs(x="Length", y="% of max", colour="Individual") +
              theme(legend.position="top")
#~               ggtitle(paste0("Individual: ", ind.nr),
#~                 subtitle=paste("Last timestep:", lasttime.var, " Survival:", round(survival.var,4),
#~                 "\nGrowth speed:", growth_speed.var, " End status:", end.status.var))
test1.plot

#Saving plot
dev.copy(png, paste0(path,"test1.png"), width=400,height=800)
dev.off()



#Plot with all hormones for all individuals vs length
test2.plot <- filter(ind2.df, growth_speed!="NA") %>% 
              gather(key=hormone.type, value=hormone.level, OXF, GHF, THF, environment,M_total,survival_next) %>%
              aggregate(hormone.level ~ length_rounded + hormone.type + growth_speed, ., mean) %>%
              mutate(hormone.type = factor(hormone.type, levels = plot.order)) %>% 
              arrange(hormone.type) %>% 
              mutate(hormone.level = case_when(
                hormone.type == "THF" ~ (hormone.level/THF.limits[2])*100,
                hormone.type == "GHF" ~ (hormone.level/GHF.limits[2])*100,
                hormone.type == "OXF" ~ (hormone.level/OXF.limits[2])*100,
                hormone.type == "environment" ~ ((hormone.level-1)/(max(ind2.df$start_environment, na.rm=TRUE)-1)*100),
                hormone.type == "M_total" ~ (hormone.level/2*100),                
                hormone.type == "survival_next" ~ (hormone.level*100)
                )) %>%             
              ggplot(aes(x=length_rounded,y=hormone.level, colour=growth_speed)) +
              geom_line() + theme_classic(base_size=16) +
              facet_wrap(~hormone.type,ncol=1,nrow=6,scales="free_y") +
              coord_cartesian(ylim=c(0, 100)) +
#~               scale_y_continuous(limits = c(0, 100)) + 
              labs(x="Length", y="% of max", colour="Growth speed") +
              theme(legend.position="top")
test2.plot

#Saving plot
dev.copy(png, paste0(path,"test2.png"), width=400,height=800)
dev.off()



#Individual mortality for 2 individuals at once with facet_wrap
mort.order <- c("M_size_O2_foraging","M_size_O2", "M_size_foraging", "M_size","M_sizeindependent")

test3.plot <- filter(ind2.df, individual %in% ind.nr) %>% 
              gather(key=mortality_type, value=mortality_data, M_size:M_size_O2_foraging) %>%
              mutate(mortality_type = factor(mortality_type, levels = mort.order)) %>% 
              arrange(mortality_type) %>%
              aggregate(mortality_data ~ length+mortality_type+individual, ., mean) %>%
              ggplot(aes(x=length, y=mortality_data, fill=mortality_type)) +
                    geom_area(colour="black", size=.2, alpha=.4) +
                    theme_classic(base_size=16) + 
                    coord_cartesian(ylim=c(0, 1.5)) +
                    facet_wrap(~as.factor(individual), ncol=1) +
                    labs(x="Length", y="Mortality (per year)", fill="Mortality type")
test3.plot

#Saving plot
dev.copy(png, paste0(path,"test3.png"), width=800,height=800)
dev.off()




#Function for simple plot with x, y and colour for a chosen length
plot.simple <- function(dataset.df,xaxis, xlabel, yaxis, ylabel, colgroup,plot.length){
  easy.plot <- ggplot(filter(dataset.df, length_rounded==plot.length), aes_string(x=xaxis, y=yaxis, colour=colgroup)) +
                geom_point() +
                theme_classic(base_size=16) +
                labs(x=xlabel,y=ylabel) +
                ggtitle(paste("Length =", plot.length))
}




test4.plot <- plot.simple(ind2.df, "environment", "Environment", "R_cat", "Reserve Fullness (%)", "surplus_before_growth", 20)
test4.plot

#Saving plot
dev.copy(png, paste0(path,"test4.png"), width=800,height=800)
dev.off()

#Foraging mortality and OXF per environment
test5.plot <- plot.simple(ind2.df,"environment", "Environment", "M_size_foraging","Size dependent foraging mortality (per year)", "OXF", 20) + scale_colour_gradient(low="purple", high="red") + coord_cartesian(ylim=c(0, 2))
test5.plot

#Saving plot
dev.copy(png, paste0(path,"test5.png"), width=800,height=800)
dev.off()

#O2 mortality and THF per environment
test6.plot <- plot.simple(ind2.df,"environment", "Environment", "M_size_O2","Size dependent oxygen mortality (per year)", "THF", 20) + scale_colour_gradient(low="blue", high="violet")
test6.plot

#Saving plot
dev.copy(png, paste0(path,"test6.png"), width=800,height=800)
dev.off()

#O2 mortality and THF per environment
test7.plot <- plot.simple(ind2.df,"environment", "Environment", "surplus_before_growth","Surplus before growth", "intake", 20) + geom_hline(yintercept=0) + scale_colour_gradient(low="yellow", high="red")
test7.plot

#Saving plot
dev.copy(png, paste0(path,"test7.png"), width=800,height=800)
dev.off()
