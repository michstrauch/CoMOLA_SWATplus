## Define path to your CoMOLA folder
path <- 'F:/mstrauch/CoMOLA_OPTAIN_CS1'

## Define names of objectives (fit1, fit2, etc. must correspond with your 
## specifications in models/SWAT.R)
## Think about including a line break (\n) for long names when indicator should be displayed as color or point size
fit1 <- 'P load at outlet [kg/a]' # meaningful name for objective 1
fit2 <- 'Lowflow\n[fraction of days\nbelow threshold]' # meaningful name for objective 2
fit3 <- 'Crop production\n[sum of grain units/a]' # meaningful name for objective 3
fit4 <- 'NSWRM costs [???/a]' # meaningful name for objective 4

## Execute the code below (do not modify)

# get functions
setwd(paste0(path,'/output_analysis'))
source('functions_postprocessing.R')

foo1(c('mco', 'dplyr', 'tidyverse', 'ggplot2', 'viridis'))

## extract pareto results
pareto <- get_pareto()

## write status quo fitness values
write_sq()

## calculate hypervolume development
HV <- hv_generations()

## plot hypervolumes for each generation
ggplot(HV, aes(Generation, HV)) +
  geom_line() +
  geom_point()

## plot Pareto solutions
plot_2D(mode=19)

# Use mode to modify the assignment of objectives to axis
# for 2-dimensional plots
# mode = 1: x = fit1, y = fit2
# mode = 2: x = fit2, y = fit1
#
# Use mode to modify the assignment of objectives to axis
# for 3-dimensional plots
# mode = 1: x = fit1, y = fit2, color = fit3
# mode = 2: x = fit1, y = fit3, color = fit2
# mode = 3: x = fit2, y = fit1, color = fit3
# mode = 4: x = fit2, y = fit3, color = fit1
# mode = 5: x = fit3, y = fit1, color = fit2
# mode = 6: x = fit3, y = fit2, color = fit1
#
# for 4-dimensional plots
# mode = 1: x = fit1, y = fit2, color = fit3, size = fit4
# mode = 2: x = fit1, y = fit2, color = fit4, size = fit3
# mode = 3: x = fit1, y = fit3, color = fit2, size = fit4
# mode = 4: x = fit1, y = fit3, color = fit4, size = fit2
# mode = 5: x = fit1, y = fit4, color = fit2, size = fit3
# mode = 6: x = fit1, y = fit4, color = fit3, size = fit2
# mode = 7: x = fit2, y = fit1, color = fit3, size = fit4
# mode = 8: x = fit2, y = fit1, color = fit4, size = fit3
# mode = 9: x = fit2, y = fit3, color = fit1, size = fit4
# mode = 10: x = fit2, y = fit3, color = fit4, size = fit1
# mode = 11: x = fit2, y = fit4, color = fit1, size = fit3
# mode = 12: x = fit2, y = fit4, color = fit3, size = fit1
# mode = 13: x = fit3, y = fit1, color = fit2, size = fit4
# mode = 14: x = fit3, y = fit1, color = fit4, size = fit2
# mode = 15: x = fit3, y = fit2, color = fit1, size = fit4
# mode = 16: x = fit3, y = fit2, color = fit4, size = fit1
# mode = 17: x = fit3, y = fit4, color = fit1, size = fit2
# mode = 18: x = fit3, y = fit4, color = fit2, size = fit1
# mode = 19: x = fit4, y = fit1, color = fit2, size = fit3
# mode = 20: x = fit4, y = fit1, color = fit3, size = fit2
# mode = 21: x = fit4, y = fit2, color = fit1, size = fit3
# mode = 22: x = fit4, y = fit2, color = fit3, size = fit1
# mode = 23: x = fit4, y = fit3, color = fit1, size = fit2
# mode = 24: x = fit4, y = fit3, color = fit2, size = fit1