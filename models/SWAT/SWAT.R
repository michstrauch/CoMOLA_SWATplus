setwd("C:/+PAPER_WORK+/Opti-Tool/CoMOLA_CS1_240502/models/SWAT")
sink("C:/+PAPER_WORK+/Opti-Tool/CoMOLA_CS1_240502/models/SWAT/console.txt", append=FALSE)

#-------------------------------------------------------------------------------
### Script to run SWAT+ within CoMOLA
#-------------------------------------------------------------------------------

### 1 - Load libraries and functions -------------------------------------------
source('calc_opt_indis.R')

library(data.table)
library(dplyr)
library(purrr)
library(readr)
library(stringr)

library(SWATmeasR)

### 2 - Define paths -----------------------------------------------------------
wd <- getwd()
txt_path <- paste0(wd,'/txt') # path to SWAT+ model txt folder

### 3 - Implement measures and run SWAT ---------------------------------
# Load the measR project which is located in the project path.
load_measr(paste0(txt_path, '/', measr_file))
# assign the data of the measr project with a specific name to the generic 
# variable with the name 'measr'
assign('measr', get(gsub('.measr$', '', measr_file)))

# Reset SWAT files
measr$reset()

# Read genome
genome <- read.csv('genom.csv', header=T)

# Define HRUs subject of NSWRM implementation
idx <- which(genome == 2)

# Implement NSWRMs
if(is.integer0(idx) == F){
  measr$implement_nswrm(nswrm_id = idx)
  measr$write_swat_inputs()
}

# Run SWAT
set2wd(txt_path)

swat_exe <- list.files(txt_path, '.exe$')
system(swat_exe)

### 4 - Calculate indicators ---------------------------------------------------
##
#
# If you want to consider crop yield as an optimization objective, specify
# grain units to normalize the basin wide sum of crop yields by crop-specific
# nutritional values, please specify grain units for relevant crops
# The grain units must be applicable to dry mass!!!
grain_units <- data.frame('wbar' = 1.163, 
                          'csil' = 1.071, 
                          'wwht' = 1.209, 
                          'wira' = 1.429,
                          'barl' = 1.163,
                          'akgs' = 0.682, 
                          'wiry' = 1.174, 
                          'sgbt' = 1)

# Define objectives. Please keep the naming syntax with fit1, fit2, ...
fit1 <- ind_cha_aa(project_path, 'cha0926')[3] * -1 #loads should be minimized
fit2 <- ind_cha_day(project_path, 'cha0926', 'Q_p05')[7]
fit3 <- ind_bsn_aa_crp(project_path, names(grain_units), out_type = "yield", grain_units)[1]

# Add the fit variables here. Please do not rename 'out'.
out <- t(cbind.data.frame(fit1, fit2, fit3))

write.table(out, paste0(wd,'/SWAT_output.csv'), 
            row.names = F, quote= F, col.names = F)







sink()
