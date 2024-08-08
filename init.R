# init.R ------------------------------------------------------------------
# 
# Version 0.1.0 
# Date: 2024-05-17 
# Author: Christoph Schürz christoph.schuerz@ufz.de
#
# init.R runs before starting the CoMOLA routine. It performs checks, to 
# catch potential issues which would cause problems in the CoMOLA routine.
# The bottom part generates and writes the 
# -------------------------------------------------------------------------

# Load required R packages ------------------------------------------------
# Check if any of the later required R packages is not installed.
if(!'SWATmeasR' %in% installed.packages()) {
  write(paste0("SWATmeasR is not installed!", 
               "Please install before re-running."),
        'error_log.txt', append = TRUE)
  stop()
}
if(!'data.table' %in% installed.packages()) {
  install.packages('data.table')
}
if(!'dplyr' %in% installed.packages()) {
  install.packages('dplyr')
}
if(!'purrr' %in% installed.packages()) {
  install.packages('purrr')
}
if(!'readr' %in% installed.packages()) {
  install.packages('readr')
}
if(!'stringr' %in% installed.packages()) {
  install.packages('stringr')
}
if(!'readxl' %in% installed.packages()) {
  install.packages('readxl')
}

library(SWATmeasR)
library(stringr)
library(readxl)

# Source economic function to access the function create_spi_pars_df()
source(paste0(getwd(), '/models/SWAT/calc_spi_indis.R'))

# Paths and parameters ----------------------------------------------------
# Path where the SWAT txt folder and the SWAT.R script are located
txt_path <- file.path(getwd(),'models/SWAT/txt')

# Path to the SWAT.R script
script_path <- file.path(getwd(), 'models/SWAT/SWAT.R')

# Minimum required version of SWATmeasR
measr_version_min <- '0.8.0'

# -------------------------------------------------------------------------
# Perform checks
# -------------------------------------------------------------------------
# Check if SWAT txt folder is available -----------------------------------
if(!dir.exists(txt_path)) {
  write(paste("SWAT model folder 'txt' is missing in",
              file.path(getwd(),'models/SWAT')),
        'error_log.txt', append = TRUE)
  stop()
}

# Checks on SWATmeasR projects and versions -------------------------------
# Find all files with file extension .measr in txt folder.
measr_file <- list.files(txt_path, '.measr$')

# Project must contain only one measR project. Report if otherwise.
if (length(measr_file) > 1) {
  write(paste('Multiple SWATmeasR projects found.', 
              'Project folder must contain only ONE measR project!'), 
        'error_log.txt', append = TRUE)
}

# Check measR project version
# Due to some updates in the SWATmeasR code a CoMOLA project requires at least
# a version defined in measr_version_min.
# 
# Load the measR project.
load_measr(paste0(txt_path, '/', measr_file))

# assign the data of the measr project with a specific name to the generic 
# variable with the name 'measr'
assign('measr', get(gsub('.measr$', '', measr_file)))

# Get the measR version of the measR project.
measr_version <- measr$.data$meta$measr_version

# Check if the version of the measR project is at least measr_version_min
if(is.null(measr_version)) {
  write(paste0('SWATmeasR which was used to build current measR project ',
               'had a version <= ', measr_version_min, '. ',
               'The CoMOLA workflow however requires a version >= ', 
               measr_version_min, '!'), 
        'error_log.txt', append = TRUE)
} else if (measr_version < measr_version_min) {
  write(paste0('SWATmeasR which was used to build current measR project ', 
               'had the version ',  measr_version, '. ',
               'The CoMOLA workflow however requires a version >= ', 
               measr_version_min, '!'), 
        'error_log.txt', append = TRUE)
}

# Check if the version of the installed SWATmeasR is at least 
# measr_version_min
# 
measr_version_installed <- as.character(packageVersion('SWATmeasR'))

if(measr_version_installed < measr_version_min) {
  write(paste0('The installed version of SWATmeasR is ', 
               measr_version_installed, '. ',
               'The CoMOLA workflow however requires a version >= ', 
               measr_version_min, '!'), 
        'error_log.txt', append = TRUE)
}

# Check if only one SWAT exe is in txt folder -----------------------------
# Find executable file.
swat_exe <- list.files(txt_path, '.exe$')

# Check if only one exe exists.
if (length(swat_exe) > 1) {
  write(paste('Multiple executable files found in the project folder.', 
              'The project folder must contain only ONE executable file!'), 
        'error_log.txt', append = TRUE)
}

# -------------------------------------------------------------------------
# Build CoMOLA input folder
# -------------------------------------------------------------------------
# Build input files -------------------------------------------------------
# Build file_HRU input file
n_gene <- nrow(measr$.data$nswrm_definition$nswrm_locations)

file_hru <- data.frame(HRU = rep(1, n_gene),
                       Area_ha = 0,
                       Area_rel = 1/n_gene)

# Build worst_fitness_values_maximize
r_script <- readLines(script_path)

pos_out <- str_which(str_trim(r_script), '^out')
n_crit  <- str_count(r_script[pos_out], 'fit')

w_fit <- rep('0', n_crit)

# Write input files -------------------------------------------------------
if(!dir.exists(file.path(getwd(), 'input'))) {
  dir.create(file.path(getwd(), 'input'))
}
write.csv(file_hru,
          file = file.path(getwd(), 'input', 'file_HRU.csv'),
          row.names = FALSE)
writeLines(w_fit,
           file.path(getwd(), 'input', 'worst_fitness_values_maximize.txt'))

# -------------------------------------------------------------------------
# Build economic model input file
# -------------------------------------------------------------------------
# Find all files with file extension .measr in txt folder.
econ_file <- list.files(paste0(getwd(), '/models/SWAT/economic_model'), 
                        '.xlsx$', full.names = TRUE)

# Project must contain only one measR project. Report if otherwise.
if (length(econ_file) > 1) {
  write(paste('Multiple Excel files were found in the folder economic_model.', 
              'The folder must contain only ONE Excel input file!'), 
        'error_log.txt', append = TRUE)
}

# Create the spi input data.frames and write them into the script CS_input_data.R
cs_input_path <- paste0(getwd(), '/models/SWAT/economic_model/CS_input_data.R')

create_spi_pars_df(pars_xlsx_path = econ_file, 
                   out_fn = cs_input_path)
  

