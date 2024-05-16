library(SWATmeasR)

txt_path <- file.path(getwd(),'models/SWAT')
measr_version_min <- '0.8.0'

if(!dir.exists('./models/SWAT/txt')) {
  write(paste0("SWAT model folder 'txt' is missing in ", txt_path),
        'error_log.txt', append = TRUE)
  stop()
}

measr_file <- list.files(txt_path, '.measr$')

if (length(measr_file) > 1) {
  write('Multiple SWATmeasR projects found. Project folder must contain only ONE measR project!', 
        'error_message.txt', append = TRUE)
}

# Check measR project version
# Due to some updates in the SWATmeasR code a CoMOLA project requires at least
# measr_version
load_measr(paste0(txt_path, '/', measr_file))
# assign the data of the measr project with a specific name to the generic 
# variable with the name 'measr'
assign('measr', get(gsub('.measr$', '', measr_file)))

measr_version <- measr$.data$meta$measr_version
if(is.null(measr_version)) {
  write(c('SWATmeasR which was used to build current measR project had a version <=', 
          measr_version_min, '!',
          'The CoMOLA workflow however requires a version >=', measr_version_min), 
        'error_message.txt', append = TRUE)
} else if (measr_version <= '0.7.0') {
  write(c('SWATmeasR which was used to build current measR project had the version ', 
          measr_version, '!',
          'The CoMOLA workflow however requires a version >=', measr_version_min), 
        'error_message.txt', append = TRUE)
}

# Find executable file and check if only one exe exists.
swat_exe <- list.files(txt_path, '.exe$')

if (length(swat_exe) > 1) {
  write(c('Multiple executable files found in project folder.', 
          'Project folder must contain only ONE executable file!'), 
        'error_message.txt', append = TRUE)
}