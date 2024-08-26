## Define path to your CoMOLA folder
path <- 'C:/+PAPER_WORK+/Opti-Tool/CoMOLA_CS1_240503'

## Define names of objectives (fit1, fit2, etc. must correspond with your 
## specifications in models/SWAT.R)
fit1 <- 'Objective 1' # give a more meaningful name for objective 1
fit2 <- 'Objective 2' # give a more meaningful name for objective 2
fit3 <- 'Objective 3' # give a more meaningful name for objective 3
fit4 <- 'Objective 4' # give a more meaningful name for objective 4

## Execute the code below (do not modify)

# get functions
setwd(paste0(path,'/output_analysis'))
source('functions_postprocessing.R')

foo1(c('mco', 'dplyr', 'tidyverse', 'ggplot2', 'viridis'))

## extract results
pareto <- get_pareto()

## calculate hypervolume development
HV <- hv_generations()

## plot hypervolumes for each generation
ggplot(HV, aes(Generation, HV)) +
  geom_line() +
  geom_point()

## plot Pareto solutions
plot_2D(mode=2)

# Use mode to modify the assignment of objectives to axis
# for 3-dimensional plots
# mode = 1: x = fit1, y = fit2, color = fit3
# mode = 2: x = fit2, y = fit3, color = fit1
# mode = 3: x = fit3, y = fit1, color = fit2
#
# for 4-dimensional plots
# mode = 1: x = fit1, y = fit2, color = fit3, size = fit4
# mode = 2: x = fit2, y = fit3, color = fit4, size = fit1
# mode = 3: x = fit3, y = fit4, color = fit1, size = fit2
# mode = 4: x = fit4, y = fit1, color = fit2, size = fit3