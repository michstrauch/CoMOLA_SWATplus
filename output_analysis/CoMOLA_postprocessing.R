if(! require('here' , character.only = TRUE)){
  #  If package was not able to be loaded then re-install
  install.packages('here' , dependencies = TRUE)
  #  Load package after installing
  require('here' , character.only = TRUE)
}

## Define paths
# to optimization results
opt_path <- here('output')

# to some post-processing folder
post_path <- here('output_analysis')

# get functions
setwd(post_path)
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
plot_2D()

