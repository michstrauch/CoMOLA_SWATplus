## Define path to your CoMOLA folder
path <- 'F:/bonato/CoMOLA_CS1_m'

## Define names of objectives (fit1, fit2, etc. must correspond with SWAT.R)
fit1 <- 'Habitat connectivity'
fit2 <- 'Habitat quality'
fit3 <- 'P load (kg/a)'
fit4 <- 'Agr. production (grain units)'

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
plot_2D(mode=3)

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