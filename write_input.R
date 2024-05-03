# Load R packages ---------------------------------------------------------
library(SWATmeasR)
library(tibble)

# Paths -------------------------------------------------------------------
project_path <- 'C:/Users/schuerz/Documents/optain/optimization/240502_CoMOLA_CS1'
txt_path <- file.path(project_path, 'models', 'SWAT', 'txt')

n_crit <- 3

# Find and load SWATmeasR file --------------------------------------------
measr_file <- list.files(txt_path, '.measr$')

# Check if only one measr file contained in project folder
if (length(measr_file) > 1) {
  write('Multiple SWATmeasR projects found. Project folder must contain only ONE measR project!',
        'error_message.txt', append = TRUE)
  stop('Multiple SWATmeasR projects found. Project folder must contain only ONE measR project!')
}

# Load the measR project which is located in the project path.
load_measr(paste0(txt_path, '/', measr_file))
# assign the data of the measr project with a specific name to the generic
# variable with the name 'measr'
assign('measr', get(gsub('.measr$', '', measr_file)))

# Build input files -------------------------------------------------------
# Build file_HRU input file
n_gene <- nrow(measr$.data$nswrm_definition$nswrm_locations)

file_hru <- tibble(HRU = rep(1, n_gene),
                   Area_ha = 0,
                   Area_rel = 1/n_gene)

# Build worst_fitness_values_maximize
w_fit <- rep('0', n_crit)

# Write input files -------------------------------------------------------
if(!dir.exists(file.path(project_path, 'input'))) {
  dir.create(file.path(project_path, 'input'))
  write.csv(file_hru,
            file = file.path(project_path, 'input', 'file_HRU.csv'),
            row.names = FALSE)
  writeLines(w_fit,
             file.path(project_path, 'input', 'worst_fitness_values_maximize.txt'))
}
