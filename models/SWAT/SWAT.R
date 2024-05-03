setwd("C:/+PAPER_WORK+/Opti-Tool/CoMOLA_CS1_240502/models/SWAT")
sink("C:/+PAPER_WORK+/Opti-Tool/CoMOLA_CS1_240502/models/SWAT/console.txt", append=FALSE)

#-------------------------------------------------------------------------------
### Script to run SWAT+ within CoMOLA
#-------------------------------------------------------------------------------


### 1 - Load libraries and functions -------------------------------------------

source('calc_opt_indis.R')
foo1(c('dplyr' , 'readr' , 'tidyverse', 'data.table', 'remotes'))
foo2('SWATmeasR')

if (file.exists('error_message.txt')) {
  file.remove('error_message.txt')
}

# parent_path <- dirname(dirname(getwd()))
# if(!'input' %in% list.dirs(parent_path, recursive = FALSE)) {
#   input_path <- paste0(parent_path, 'input')
#   dir.create(input_path)
#   writeLines(rep(0, 3), paste0(input_path, 'worst_fitness_values_maximize.txt'))
#   writeLines(rep(0, 3), paste0(input_path, 'worst_fitness_values_maximize.txt'))
# }

### 2 - Define paths -----------------------------------------------------------

wd <- getwd()
project_path <- paste0(wd,'/txt') # path to SWAT+ model txt folder

### 3 - Implement measures and run SWAT ---------------------------------

measr_file <- list.files(project_path, '.measr$')

if (length(measr_file) > 1) {
  write('Multiple SWATmeasR projects found. Project folder must contain only ONE measR project!', 
        'error_message.txt', append = TRUE)
  stop('Multiple SWATmeasR projects found. Project folder must contain only ONE measR project!')
}

# Load the measR project which is located in the project path.
load_measr(paste0(project_path, '/', measr_file))
# assign the data of the measr project with a specific name to the generic 
# variable with the name 'measr'
assign('measr', get(gsub('.measr$', '', measr_file)))

# Check measR project version
# Due to some updates in the SWATmeasR code a CoMOLA project requires at least
# version 0.7.0
measr_version <- measr$.data$meta$measr_version
if(is.null(measr_version)) {
  write(c('SWATmeasR which was used to build current measR project had a version <= 0.7.0!',
          'The CoMOLA workflow however requires a version >= 0.7.0'), 
          'error_message.txt', append = TRUE)
  stop('SWATmeasR which was used to build current measR project <= 0.7.0!\n',
       'The CoMOLA workflow however requires a version >= 0.7.0')
} else if (measr_version <= '0.7.0') {
  write(c('SWATmeasR which was used to build current measR project had a version <= 0.7.0!',
          'The CoMOLA workflow however requires a version >= 0.7.0'), 
        'error_message.txt', append = TRUE)
  stop('SWATmeasR which was used to build current measR project <= 0.7.0!\n',
       'The CoMOLA workflow however requires a version >= 0.7.0')
}

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
# Find executable file and check if only one exe exists.
swat_exe <- list.files(project_path, '.exe$')

if (length(swat_exe) > 1) {
  write(c('Multiple executable files found in project folder.', 
          'Project folder must contain only ONE executable file!'), 
        'error_message.txt', append = TRUE)
  stop('Multiple executable files found in project folder')
}

set2wd(project_path)
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

fit1 <- ind_cha_aa(project_path, 'cha0926')[3] * -1 #loads should be minimized
fit2 <- ind_cha_day(project_path, 'cha0926', 'Q_p05')[7]
fit3 <- ind_bsn_aa_crp(project_path, names(grain_units), out_type = "yield", grain_units)[1]

out <- t(cbind.data.frame(fit1, fit2, fit3))

write.table(out, paste0(wd,'/SWAT_output.csv'), row.names = F, quote= F, col.names = F)







sink()
