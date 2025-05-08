setwd("C:/enricodata/progetto_OPTAIN/COMOLA/models/SWAT")
sink("C:/enricodata/progetto_OPTAIN/COMOLA/models/SWAT/console.txt", append=FALSE)

#-------------------------------------------------------------------------------
### Script to run SWAT+ within CoMOLA
#-------------------------------------------------------------------------------


### 1 - Load libraries and functions -------------------------------------------

source('calc_opt_indis.R')

source('calc_spi_indis.R') # add here spi library

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
msr_impl_costs_df = data.frame(
  "name" = c("buffer", "channres", "droughtplt", "pond", "terrace"),
  "value" = c(9973, 3314.68, 0, 20000, 5223)
)

msr_mnt_costs_df = data.frame(
  "name" = c("buffer", "channres", "droughtplt", "pond", "terrace"),
  "value" = c(NA, 570, NA, 3986, 3252.6)
)

msr_subsidies_df = data.frame(
  "name" = c("buffer", "channres", "droughtplt", "pond", "terrace"),
  "value" = c(0, 0, 0, 0, 392.5102)
)

## measure indicators ----
# can be run without SWAT.exe

fit1 <- calc_measr_impl_cost_ind(measr_obj = measr, # a measr object or the complete path to the measr file
                                 impl_nswrms_list = idx, # = idx from SWAT.R
                                 measr_IC = msr_impl_costs_df)      # implementation costs: a dataframe with two column: name, value
                                 
fit2 <- calc_measr_maint_cost_ind(measr_obj = measr, # a measr object or the complete path to the measr file
                                  impl_nswrms_list = idx, # = idx from SWAT.R
                                  measr_MC  = msr_mnt_costs_df)     # maintenance costs: a dataframe with two column: name, value
                                              
fit3 <- calc_measr_subsidies_ind(measr_obj = measr, # a measr object or the complete path to the measr file
                                impl_nswrms_list = idx, # = idx from SWAT.R 
                                measr_SB  = msr_subsidies_df)     # subsidies: a dataframe with two column: name, value

fit4 <- calc_measr_all_ind(measr_obj = measr, # a measr object or the complete path to the measr file
                          impl_nswrms_list = idx, # = idx from SWAT.R
                          measr_IC = msr_impl_costs_df,     # implementation costs: a dataframe with two column: name, value
                          measr_MC = msr_mnt_costs_df,     # maintenance costs: a dataframe with two column: name, value
                          measr_FP = msr_subsidies_df)       # subsidies: a dataframe with two column: name, value

## carbon sequestration ----

# Put the following code before calling the function: 
measr_hru_con_df = measr[[".data"]][["model_setup"]][["modified_inputs"]][["hru.con"]]
measr_hru_data_df = measr[[".data"]][["model_setup"]][["modified_inputs"]][["hru_data.hru"]]
measr_hru_data_df = measr_hru_data_df %>% left_join(measr_hru_con_df %>% select(id, area), by = "id")
measr_plant_ini_df = measr[[".data"]][["model_setup"]][["modified_inputs"]][["plant.ini"]]

tot_init_biom = calc_biom_init(measr_hru_data_df,measr_plant_ini_df)

env_prices_df = data.frame(
  "name" = c("N", "sed", "C", "P"),
  "value" = c(14, 12, 86.53, 15.67)
)

c_seq_performance_df = data.frame(
  "name" = c("urld", "past", "frst", "orcd", "rnge", "urml", "ucom", "field"),
  "value" = c(0, 0.9, 0.9, 0.9, 0.9, 0, 0, 0.1)
)

fit5 = calc_carbon_ind(path = project_path,
                       measr_hru_data_df= measr_hru_data_df,
                       env_prices_df = env_prices_df,
                       c_seq_per_df  = c_seq_performance_df,
                       tot_init_biom = tot_init_biom)

## water quality ----

fit6 = calc_water_quality_ind(path = project_path,
                              hru_data_df = measr_hru_data_df,
                              env_prices_df = env_prices_df)


## farm indicator ----

crop_prod_costs_mgt_df = data.frame(
  "op_typ" = c("harv", "kill", "fert", "till", "till", "plnt", "till", "fert", "till", "harv", "kill", "till", "plnt", "skip", "harv", "till", "irrm", "kill", "plnt", "fert", "plnt", "harv", "kill", "plnt", "harv", "kill", "fert", "harv", "kill", "fert", "till", "plnt", "fert", "harv", "kill", "plnt", "plnt", "harv", "kill", "plnt", "harv", "kill", "fert", "plnt", "harv", "fert", "kill", "harv", "harv", "harv"),
  "op_data1" = c("mustard", "mustard", "liq_man_it", "fallplow", "harrow", "grsg", "rolpkrat", "urea", "rodweedr", "grsg", "grsg", "rollhrrw", "alfa", "null", "alfa", "packer", "irrig", "alfa", "mustard", "00_46_00", "soy1", "soy1", "soy1", "csi1", "csi1", "csi1", "urea", "wwht", "wwht", "liq_man_it", "furwdike", "ryeg", "26_00_00", "ryeg", "ryeg", "wwht", "soy2", "soy2", "soy2", "csi2", "csi2", "csi2", "00_46_00", "hay", "hay", "urea", "hay", "frsdit", "orcd", "hayd"),
  "op_data2" = c("grass_mulch", "null", "aerial_liquid", "null", "null", "null", "null", "broadcast", "null", "silage", "null", "null", "null", "null", "hay_cut_high", "null", "null", "null", "null", "broadcast", "null", "grain", "null", "null", "silage", "null", "broadcast", "silage", "null", "aerial_liquid", "null", "null", "broadcast", "hay_cut_high", "null", "null", "null", "grain", "null", "null", "silage", "null", "broadcast", "null", "hay_cut_high", "broadcast", "null", "forest_cut", "orchard", "hay_cut_low"
  ),
  "op_data3" = c(0, 0, 4800, 0, 0, 0, 0, 423, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 249, 0, 0, 0, 0, 0, 0, 120, 0, 0, 40000, 0, 0, 150, 0, 0, 0, 0, 0, 0, 0, 0, 0, 150, 0, 0, 110, 0, 0, 0, 0),
  "value" = c(110, NA, NA, NA, NA, NA, NA, NA, NA, 675, NA, NA, NA, NA, 816, NA, NA, NA, NA, NA, NA, 847.1, NA, NA, 1420, NA, NA, 557, NA, NA, NA, NA, NA, 377, NA, NA, NA, 656, NA, NA, 1086, NA, NA, NA, 717, NA, NA, 136.666666666667, 1467.89, 280)
)


# crop cost calculated directly from management
fit7 = calc_mgt_cost_ind(path = project_path,
                         measr_obj = measr,
                         mgt_costs_df = crop_prod_costs_mgt_df)


crop_prod_costs_df = data.frame(
  "name" = c("frsd", "orcd", "rnge", "must", "grsg", "alfa", "soy1", "soy2", "wwht", "ryeg", "hay", "hayd", "csi1", "csi2", "soyb", "csil", "grsd", "grsw"),
  "value" = c(136.666666666667, 1467.89, 0, 110, 675, 816, 929, 656, 557, 377, 717, 280, 1420, 1086, 847.1, 1319.8, 675, 675)
)

fit8 = calc_crop_prod_cost_ind(path = project_path,
                              crp_prd_costs_df = crop_prod_costs_df)

crop_farm_payms_df = data.frame(
  "name" = c("frsd", "orcd", "rnge", "buff", "must", "grsg", "alfa", "soy1", "soy2", "wwht", "ryeg", "hay", "hayd", "csi1", "csi2", "soyb", "csil", "grsd", "grsw"),
  "value" = c(0, 261, 0, 600, 300, 0, 0, 0, 0, 0, 0, 110, 110, 0, 0, 0, 0, 0, 0)
)


fit9 = calc_crop_paymets_ind(path = project_path,
                            crp_frm_pays_df = crop_farm_payms_df)

crop_prices_df = data.frame(
  "name" = c("frsd", "orcd", "rnge", "rngb", "must", "grsg", "alfa", "soy1", "soy2", "wwht", "ryeg", "hay", "hayd", "csi1", "csi2", "soyb", "csil", "grsd", "grsw"),
  "value" = c(207.778, 5295.24, 31.94, 283.26, 0, 126.58, 222.75, 488.31, 488.31, 116.62, 205.52, 185.63, 168.75, 137.12, 137.12, 488.31, 137.12, 126.58, 126.58)
)

fit10 = calc_crop_values_ind(path = project_path,
                             crp_prc_df = crop_prices_df)


# set crp_prd_costs_df to empty dataframe and offset to fit7 to include management costs
# instead of harvest area costs
fit11 = calc_agr_gross_marg_ind(path = project_path,
                                crp_prd_costs_df = data.frame(name = character(),
                                                              value = double()),
                                crp_frm_pays_df = crop_farm_payms_df,
                                crp_prc_df = crop_prices_df,
                                offset = fit7)



out <- t(cbind.data.frame(fit1, fit2, fit3, fit10))

write.table(out, paste0(wd,'/SWAT_output.csv'), row.names = F, quote= F, col.names = F)


sink()
