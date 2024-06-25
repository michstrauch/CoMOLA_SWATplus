# dear user,
# define here the parameters necessary to the model
# each parameter table is a dataframe
# use create_spi_pars_df to create the initialization script from the spreadsheet

# REQUIREMENTS ----
# read_table, %nin% from calc_opt_indis.R
#library('readxl')

# ANCILLARY ----
# create_spi_pars_df 
# the function reads the parameter required by spi indicator from the spreadsheet file
# and create an R script that can be used for initialization
# with all the parameters ready to be included in data.frame
# the name of the data.frames are the same as the defaults in the indicator function 
# NOTE: use once or every time the spreadsheet is updated
# USAGE
# create_spi_pars_df('./output_to_send_CS9/CS9_input data.xlsx','./output_to_send_CS9/CS9_input data.r')

create_spi_pars_df <- function(pars_xlsx_path, out_fn='spi_init.r',
                               cols_to_use=c('name','op_typ','op_data1','op_data2','op_data3','value')){
  
  # read all the sheets in the in_fn (spreadsheet type)
  sheet_name_list = c('crop_prices','crop_farm_payms','crop_prod_costs','crop_prod_costs_mgt',
                      'msr_impl_costs','msr_mnt_costs','msr_subsidies','env_prices','c_seq_performance','TIC','ST')
  
  df_txt_list = c()
  
  for (sheet_name in sheet_name_list){
    df = read_excel(pars_xlsx_path,sheet=sheet_name)
    # get col names
    col_name_list = names(df)
    col_txt_list = c()
    for (col_name in col_name_list){
      if (any(cols_to_use==col_name)){
        col_txt = paste0('"',col_name, '" = ',paste0(df[col_name]),'')
        col_txt_list = c(col_txt_list,col_txt)
        }
      }
    df_txt = paste0(sheet_name,'_df = data.frame(\n',paste0(col_txt_list, collapse = ',\n'),'\n)')
    df_txt_list = c(df_txt_list,df_txt)
  }
  
  txt = paste(df_txt_list, collapse = '\n\n')
  txt = paste0('# spi parameters\n# origin: ',pars_xlsx_path,'\n\n',txt)
  write(txt, file = out_fn,
        append = FALSE, sep = "")
}

get_parametric_value <-function(par_name, # the name of the parametric value
                                par_table, # a dataframe with two column: name, value
                                par_default = 0.){
  # try with the entire name
  par_name_gen = gsub('_\\d+', '', par_name)
  par_name_list = c(par_name, par_name_gen)
  par_value = par_default
  for (par_name_test in par_name_list){
    par_value = par_table$value[sub(par_name_test, '',par_table$name)=='']
    if (length(par_value)>0)  break
    else par_value = par_default
  }
  
  # empty parameter values are set to default 
  if (is.na(par_value)) par_value = par_default
  
  return(par_value)
}

calc_mgt_cost <-function(hru_df,sim_df,n_skip_yrs,mgt_df,crp_cst_df){
  # keep used fields only
  hru_df = hru_df[,c('id','name', 'area', 'lu_mgt')]
  
  # fix mgt code
  hru_df$lu_mgt=sub(pattern='_lum', replacement='_mgt', hru_df$lu_mgt)
  
  # add years based on skip operations
  yrs_seq = seq(sim_df$yrc_start,sim_df$yrc_end)
  
  # add years
  n_skip_ops = length(mgt_df$name[mgt_df$op_typ=='skip'])
  n_flds = n_skip_ops/length(yrs_seq) # number of fields in management schedule
  
  mgt_df['yr'] = NA
  
  only_skip = subset(mgt_df, mgt_df$op_typ=='skip')
  count_skip = aggregate(only_skip$name,list(only_skip$name),length)
  names(count_skip)=c('name','x')
  
  min_count = min(count_skip$x)
  max_count = max(count_skip$x)
  
  if (min_count != max_count) {warning('N. of skip operations differs between managements')}
  
  for (cc in 1:length(count_skip$name)){
    mgt_name = count_skip$name[cc] 
    idx = (mgt_df$name == mgt_name) & (mgt_df$op_typ =='skip')
    # TODO: check here
    mgt_df$yr[idx][1:min_count] = yrs_seq[1:min_count]#yrs_seq[1:count_skip$x[cc]]
  }
  #mgt_df$yr[mgt_df$op_typ=='skip'] = rep(yrs_seq,n_flds)
  
  mgt_df$yr = repeat_before(mgt_df$yr,T)
  
  print(paste('n. of not allocated operation:',length(mgt_df$yr[is.na(mgt_df$yr)])))
  mgt_df = mgt_df[!is.na(mgt_df$yr),]
  
  # join with operation cost
  crp_cst_df['code']=paste(crp_cst_df$op_typ,
                           crp_cst_df$op_data1,
                           crp_cst_df$op_data2,
                           crp_cst_df$op_data3, sep='')
  
  mgt_df['code']=paste(mgt_df$op_typ,
                       mgt_df$op_data1,
                       mgt_df$op_data2,
                       mgt_df$op_data3, sep='')
  
  mgt_cost_df = merge(mgt_df[,c('code','name','mon','yr','op_typ','op_data1','op_data2','op_data3')],
                      crp_cst_df[,c('code','value')],
                      by='code',all.x = T)
  
  names(mgt_cost_df)[names(mgt_cost_df) == 'name'] <- 'name_mgt'
  names(mgt_cost_df)[names(mgt_cost_df) == 'value'] <- 'op_cost'
  
  # merge hru data
  mgt_cost_df = merge(mgt_cost_df,
                      hru_df[,c('lu_mgt','area','name')],
                      by.x='name_mgt', by.y='lu_mgt')
  
  ### CALCULATE AVERAGE TOTAL COSTS ###
  mgt_cost_df['fld_op_cost'] = mgt_cost_df$area*mgt_cost_df$op_cost
  sum_mgt_cost = sum(mgt_cost_df$fld_op_cost[mgt_cost_df$yr>=(sim_df$yrc_start+n_skip_yrs)], na.rm=T)/
    (sim_df$yrc_end-sim_df$yrc_start+1-n_skip_yrs)
  return(sum_mgt_cost)
}

# calc_biom_init #
# Calculate the total biomass at the beginning of the simulation
# To be use AFTER measr$implement_nswrm(nswrm_id = idx) but BEFORE calc_carbon_ind 
# Put the following code before calling the function: 
# measr_hru_con_df = measr[[".data"]][["model_setup"]][["modified_inputs"]][["hru.con"]]
# measr_hru_data_df = measr[[".data"]][["model_setup"]][["modified_inputs"]][["hru_data.hru"]]
# measr_hru_data_df = measr_hru_data_df %>% left_join(measr_hru_con_df %>% select(id, area), by = "id")
# measr_plant_ini_df = measr[[".data"]][["model_setup"]][["modified_inputs"]][["plant.ini"]]
calc_biom_init <- function(measr_hru_data_df,measr_plant_ini_df){
  # fix mgt code
  measr_hru_data_df$lu_mgt=sub(pattern='_lum', replacement='_com', measr_hru_data_df$lu_mgt)
  
  # select only records with bm_init > 0
  measr_plant_ini_df = measr_plant_ini_df[measr_plant_ini_df$bm_init>0,]
  
  # merge plant.ini
  measr_plant_ini_df = merge(measr_hru_data_df[,c('lu_mgt','area','name')],
                             measr_plant_ini_df[,c('pcom_name','bm_init')],
                             by.x='lu_mgt', by.y='pcom_name')
  
  
  measr_plant_ini_df <- measr_plant_ini_df %>%
    mutate(weighted_bm_init = bm_init * area)
  
  tot_init_biom = sum(measr_plant_ini_df$weighted_bm_init, na.rm= T)
  
  return(tot_init_biom)
}

# read time.sim
read_time_sim <- function(file_path){
  colNames = c('day_start','yrc_start','day_end','yrc_end','step')
  
  # read text files
  df = read_swat_table_file(file_path,colNames,2,'var',F)
  df
}

# read general text file as table
read_swat_table_file <-function(file_path, col_names = NULL,
                                skip=1,dummy = 'dummy',
                                force_col= F, header = F){
  
  if (is.null(col_names)) force_col= FALSE
  
  if (force_col==TRUE){
    df = read.delim(file_path, header = header, sep = "",
                    dec = ".", skip =skip, col.names = col_names)
    
  }else{
    df = read.delim(file_path, header = header, sep = "",
                    dec = ".", skip =skip)
    if (!is.null(col_names)){
      n_cols = ncol(df)
      n_names = length(col_names)
      n_missing = n_cols-n_names 
      cat('n_cols',n_cols,'n_names',n_names,'n_missing',n_missing,'\n')
      if (n_missing>0){
        col_names = c(col_names,
                      paste0(dummy,seq(1,n_missing))
        )  
      }
      cat('col_names',col_names,'\n')
      names(df) = col_names
    }
  }
  
  df
}

# read skip yr from print.prt
read_skip_yr <- function(file_path){
  lines = read_lines(file_path, n_max = 3)
  toks = strsplit(lines[3], split = "")
  skip_yr = as.numeric(toks[[1]][1])
  skip_yr
}

# credits: https://stackoverflow.com/questions/7735647/replacing-nas-with-latest-non-na-value
repeat_before <- function(xx, before = F) {   # repeats the last non NA value. Keeps leading NA
  if (before) x = rev(xx)# flip
  else x = xx
  
  ind = which(!is.na(x))      # get positions of nonmissing values
  if(is.na(x[1]))             # if it begins with a missing, add the 
    ind = c(1,ind)        # first position to the indices
  
  res = rep(x[ind],
            times = diff(   # repeat the values at these indices
              c(ind, length(x) + 1) )) # diffing the indices + length yields how often   
  # they need to be repeated
  
  if (before) res = rev(res)# flip
  res
}                               


# FARM RELATED INDICATORs ----

## MANAGEMENT COST ##

calc_mgt_cost_ind <-function(path,
                             measr_obj,
                             mgt_costs_df = data.frame(op_typ = character(),
                                                       op_data1 = character(),
                                                       op_data2 = character(),
                                                       op_data3 = character(),
                                                       value = double())){
  # get time extreme
  sim_df = read_time_sim(paste0(path,'/time.sim'))
  n_skip_yrs = read_skip_yr(paste0(path,'/print.prt'))
  
  # get hru info updated after measr
  measr_hru_con_df = measr_obj[[".data"]][["model_setup"]][["modified_inputs"]][["hru.con"]]
  measr_hru_data_df = measr_obj[[".data"]][["model_setup"]][["modified_inputs"]][["hru_data.hru"]]
  
  # add area to hru data
  measr_hru_data_df = measr_hru_data_df %>%
    left_join(measr_hru_con_df %>% select(id, area), by = "id")
  
  mgt_df = measr_obj[[".data"]][["model_setup"]][["modified_inputs"]][["management.sch"]]
  tot_crop_prod_costs_mgt <- calc_mgt_cost(measr_hru_data_df,
                                           sim_df,n_skip_yrs,
                                           mgt_df, mgt_costs_df)
  
  tot_crop_prod_costs_mgt= - tot_crop_prod_costs_mgt # cost are always negative so the higher is the better for optimization
}

## CROP PRODUCTION COST ##

calc_crop_prod_cost_ind <-function(path,
                                   crp_prd_costs_df = data.frame(name = character(),
                                                                 value = double())
                                   ){
  
  crop_aa_df <- read_tbl('basin_crop_yld_aa.txt', path, 2)
  
  # calculate crop production costs by simplified approach
  temp_df = merge(crop_aa_df, crp_prd_costs_df, by = intersect(names(x), names(y)),
                  by.x ='plant_name', by.y = 'name', all.x=T)
  if (anyNA(temp_df$value)) message('Missing crops in prod. costs: ',paste(temp_df$plant_name[is.na(temp_df$value)], collapse=' '),'\n')
  tot_crop_prod_costs <- sum(temp_df$`harv_area(ha)` * temp_df$value, na.rm = T)
  return(-tot_crop_prod_costs) # cost are always negative so the higher is the better for optimization
}

## CROP PAYMENTS ##

calc_crop_paymets_ind <-function(path,
                                 crp_frm_pays_df = data.frame(name = character(),
                                                                    value = double())){
  
  crop_aa_df <- read_tbl('basin_crop_yld_aa.txt', path, 2)
  
  # calculate crop farmer payments
  temp_df = merge(crop_aa_df, crp_frm_pays_df, by = intersect(names(x), names(y)),
                  by.x ='plant_name', by.y = 'name', all.x=T)
  if (anyNA(temp_df$value)) message('Missing crops in farm payments: ',paste(temp_df$plant_name[is.na(temp_df$value)], collapse=' '),'\n')
  tot_crop_frm_pays <- sum(temp_df$`harv_area(ha)` * temp_df$value, na.rm = T)
  tot_crop_frm_pays
}

## CROP VALUES ##

calc_crop_values_ind <-function(path,
                                crp_prc_df = data.frame(name = character(),
                                                        value = double())){
  
  crop_aa_df <- read_tbl('basin_crop_yld_aa.txt', path, 2)
  
  # calculate crop farmer payments
  temp_df = merge(crop_aa_df, crp_prc_df, by = intersect(names(x), names(y)),
                  by.x ='plant_name', by.y = 'name', all.x=T)
  if (anyNA(temp_df$value)) message('Missing crops in crop prices: ',paste(temp_df$plant_name[is.na(temp_df$value)], collapse=' '),'\n')
  tot_crop_values <- sum(temp_df$`yld(t)` * temp_df$value, na.rm = T)
  tot_crop_values
}

## AGRICULTURAL GROSS MARGIN ##

# NOTE: use offset to add/remove values externally calculated
#       e.g. use external function to calculate production costs
#            set crp_prd_costs_df to empty dataframe and negative offset
calc_agr_gross_marg_ind <- function(path,
                            crp_prd_costs_df = data.frame(name = character(),
                                                          value = double()),
                            crp_frm_pays_df = data.frame(name = character(),
                                                         value = double()),
                            crp_prc_df = data.frame(name = character(),
                                                    value = double()),
                            offset = 0.0){
  
  crop_aa_df <- read_tbl('basin_crop_yld_aa.txt', path, 2)
  
  # calculate crop production costs by simplified approach
  temp_df = merge(crop_aa_df, crp_prd_costs_df, by = intersect(names(x), names(y)),
                  by.x ='plant_name', by.y = 'name', all.x=T)
  if (anyNA(temp_df$value)) message('Missing crops in prod. costs: ',paste(temp_df$plant_name[is.na(temp_df$value)], collapse=' '),'\n')
  tot_crop_prod_costs <- round(sum(temp_df$`harv_area(ha)` * temp_df$value,
                                   na.rm = T),3)
  # calculate crop farmer payments
  temp_df = merge(crop_aa_df, crp_frm_pays_df, by = intersect(names(x), names(y)),
                  by.x ='plant_name', by.y = 'name', all.x=T)
  if (anyNA(temp_df$value)) message('Missing crops in farm payments: ',paste(temp_df$plant_name[is.na(temp_df$value)], collapse=' '),'\n')
  tot_crop_frm_pays <- round(sum(temp_df$`harv_area(ha)` * temp_df$value,
                                 na.rm = T),3)
  # calculate crop prices
  temp_df = merge(crop_aa_df, crp_prc_df, by = intersect(names(x), names(y)),
                  by.x ='plant_name', by.y = 'name', all.x=T)
  if (anyNA(temp_df$value)) message('Missing crops in crop prices: ',paste(temp_df$plant_name[is.na(temp_df$value)], collapse=' '),'\n')
  tot_crop_values <- round(sum(temp_df$`yld(t)` * temp_df$value,
                               na.rm = T),3)
  
  net_marg_farm = tot_crop_values + tot_crop_frm_pays - tot_crop_prod_costs + offset
  net_marg_farm
  
}

# WATER QUALITY INDICATORs----

# calc_water_quality_ind #
# hru_data_df should be set externally one time at the beginning of the COMOLA
# Put the following code before calling the function: 
# hru_con_df = measr[[".data"]][["model_setup"]][["original_inputs"]][["hru.con"]]
# hru_data_df = measr[[".data"]][["model_setup"]][["original_inputs"]][["hru_data.hru"]]
# hru_data_df = hru_data_df %>% left_join(hru_con_df %>% select(id, area), by = "id") # add area to hru data
calc_water_quality_ind <- function(path, hru_data_df,
                                   env_prices_df = data.frame(name = character(),
                                                              value = double())){
  
  hru_pw <- read_tbl('hru_pw_aa.txt', path, 3)
  
  hru_ls <- read_tbl('hru_ls_aa.txt', path, 3)
  
  # Add the HRU area to each HRU
  hru_ls <- hru_ls %>%
    left_join(hru_data_df %>% select(id, area), by = c("unit" = "id"))
  
  hru_pw <- hru_pw %>%
    left_join(hru_data_df %>% select(id, area), by = c("unit" = "id"))
  
  # Calculate the total weighted sum for N and P
  # Calculate N losses 
  # Add N_losses as a new column to df_selected_hru_nb
  
  hru_ls <- hru_ls %>%
    mutate(N_losses = sedorgn + surqno3 + lat3no3 + tileno3 + hru_pw$percn)
  
  # The P input is only fertp from hru_nb
  # Calculate the sum of P losses
  hru_ls <- hru_ls %>%
    mutate(P_losses = sedorgp + surqsolp + sedminp + tilelabp + lchlabp)
  
  # Calculate the weighted values for N, P and sediment losses
  hru_ls <- hru_ls %>%
    mutate(weighted_N_losses = N_losses * area,
           weighted_P_losses = P_losses * area,
           weighted_sed_losses = sedyld * area)
  
  # water quality improvements: Nitrogen
  tot_weighted_N_losses = sum(hru_ls$weighted_N_losses)/1000
  N_cost = get_parametric_value('N',env_prices_df)
  tot_N_cost = tot_weighted_N_losses * N_cost
  
  # water quality improvements: Phosphorous
  tot_weighted_P_losses = sum(hru_ls$weighted_P_losses)/1000
  P_cost = get_parametric_value('P',env_prices_df)
  tot_P_cost = tot_weighted_P_losses * P_cost
  
  # water quality improvements: sediments
  tot_weighted_sed_losses = sum(hru_ls$weighted_sed_losses)
  sed_cost = get_parametric_value('sed',env_prices_df)
  tot_sed_cost = tot_weighted_sed_losses * sed_cost
  
  tot_wat_quality_cost = tot_N_cost + tot_P_cost + tot_sed_cost
  
  return(-tot_wat_quality_cost) # cost are always negative so the higher is the better for optimization
}

# CARBON SEQUESTRATION INDICATORs ----

# To be use AFTER measr$implement_nswrm(nswrm_id = idx)
# Put the following code before calling the function: 
# measr_hru_con_df = measr[[".data"]][["model_setup"]][["modified_inputs"]][["hru.con"]]
# measr_hru_data_df = measr[[".data"]][["model_setup"]][["modified_inputs"]][["hru_data.hru"]]
# measr_hru_data_df = measr_hru_data_df %>% left_join(measr_hru_con_df %>% select(id, area), by = "id")
# measr_plant_ini_df = measr[[".data"]][["model_setup"]][["modified_inputs"]][["plant.ini"]]
# tot_init_biom = calc_biom_init(measr_hru_data_df, measr_plant_ini_df)
# NOTE: use "calc_biom_init" function to calculate "tot_init_biom"

calc_carbon_ind <- function(path, measr_hru_data_df,
                            env_prices_df = data.frame(name = character(),
                                                       value = double()),
                            c_seq_per_df  = data.frame(name = character(),
                                                       value = double()),
                            tot_init_biom = 0.0){

  hru_pw <- read_tbl('hru_pw_aa.txt', path, 3)
  
  # update measr_hru_data_df in order to include also C seq performance
  measr_hru_data_df['c_seq_perf'] = 0.
  for (k in seq_len(nrow(c_seq_per_df))){
    measr_hru_data_df$c_seq_perf[grep(paste0('^',c_seq_per_df$name[k]),
                                      measr_hru_data_df$lu_mgt)] = c_seq_per_df$value[k]
  }
  
  # join area
  hru_pw <- hru_pw %>%
    left_join(measr_hru_data_df %>% select(id, area,c_seq_perf), by = c("unit" = "id")) %>%
    mutate(weighted_bioms = bioms * area) %>%
    mutate(weighted_yield = yield * area) %>%
    mutate(weighted_c_seq_perf = c_seq_perf * area)
  
  weighted_c_seq_perf = sum(hru_pw$weighted_c_seq_perf)/sum(hru_pw$area)
  
  # make total biomass
  # +0.34 considers also roots
  # see average values from 
  # Mathew, I., Shimelis, H., Mutema, M., & Chaplot, V. (2017).
  # What crop type for atmospheric carbon sequestration: Results from a global data analysis.
  # Agriculture, Ecosystems & Environment, 243, 34–46. https://doi.org/10.1016/j.agee.2017.04.008
  
  tot_weighted_bioms = 1.34*sum(hru_pw$weighted_bioms)/1000
  tot_weighted_yield = sum(hru_pw$weighted_yield)/1000
  tot_init_biom = tot_init_biom/1000
  
  # remove init biomass and add performance index
  net_weighted_bioms = weighted_c_seq_perf*(tot_weighted_bioms - tot_weighted_yield - tot_init_biom)
  
  # carbon sequestration
  # Ma et al., 2018
  # Variations and determinants of carbon content in plants: a global synthesis.
  # Biogeosciences, 15, 693–702, 2018
  tot_C = net_weighted_bioms * 0.45
  tot_CO2 = tot_C * 3.67
  C_price = get_parametric_value('C',env_prices_df)
  tot_C_value = tot_CO2 * C_price
  
  return(tot_C_value)                         
}

# MEASURE IMPLEMENTATION INDICATORs ----

## NSWRMs implementation costs ##

# To be use AFTER measr$implement_nswrm(nswrm_id = idx)
calc_measr_impl_cost_ind <- function(measr_obj, # a measr object or the complete path to the measr file
                                     impl_nswrms_list, # = idx from SWAT.R
                                     measr_IC = data.frame(name = character(),
                                                           value = double())      # implementation costs: a dataframe with two column: name, value
                                     ){
  # get the hru data table
  hru_con = measr_obj$.data$model_setup$original_inputs$hru.con
  # get the ids of the hru involved by the measr scenario
  nswrms_loca <- measr_obj$.data$nswrm_definition$nswrm_locations 
  # get the ids of the hru involved by the measr scenario
  tot_hru_area = 0.
  tot_IC = 0.
  
  # Loop over all measures and calculate implementation cost
  for (impl_nswrm in impl_nswrms_list) {
    # get the name of the measure
    nswrm_name = nswrms_loca$name[nswrms_loca$id == impl_nswrm]
    # get the list of hru involved by the measure
    hru_ids = unique(unlist(nswrms_loca$obj_id[nswrms_loca$id == impl_nswrm]))
    # get the areas in hectares from each hrus
    hru_area = hru_con$area[hru_con$id %in% hru_ids]
    
    sum_hru_area = sum(hru_area)
    tot_hru_area = tot_hru_area+sum_hru_area
    
    # get the cost/payment per hectare for each measure, find the best matching label
    
    # Measure implementation costs
    IC = get_parametric_value(nswrm_name,measr_IC)
    sum_IC = sum_hru_area*IC
    tot_IC = tot_IC + sum_IC
  }
  
  return(-tot_IC) # cost are always negative so the higher is the better for optimization
}

## NSWRMs maintenance costs ##

# To be use AFTER measr$implement_nswrm(nswrm_id = idx)
calc_measr_maint_cost_ind <- function(measr_obj, # a measr object or the complete path to the measr file
                                      impl_nswrms_list, # = idx from SWAT.R
                                      measr_MC  = data.frame(name = character(),
                                                            value = double())     # maintenance costs: a dataframe with two column: name, value
  ){
  # get the hru data table
  hru_con = measr_obj$.data$model_setup$original_inputs$hru.con
  # get the ids of the hru involved by the measr scenario
  nswrms_loca <- measr_obj$.data$nswrm_definition$nswrm_locations 
  tot_hru_area = 0.
  tot_MC = 0.
  
  # Loop over all measures and calculate implementation cost
  for (impl_nswrm in impl_nswrms_list) {
    # get the name of the measure
    nswrm_name = nswrms_loca$name[nswrms_loca$id == impl_nswrm]
    # get the list of hru involved by the measure
    hru_ids = unique(unlist(nswrms_loca$obj_id[nswrms_loca$id == impl_nswrm]))
    # get the areas in hectares from each hrus
    hru_area = hru_con$area[hru_con$id %in% hru_ids]
    
    sum_hru_area = sum(hru_area)
    tot_hru_area = tot_hru_area+sum_hru_area
    
    # get the cost/payment per hectare for each measure, find the best matching label
    
    # Measure maintanance costs
    MC = get_parametric_value(nswrm_name,measr_MC)
    sum_MC = sum_hru_area*MC
    tot_MC = tot_MC + sum_MC
  }
  
  return(-tot_MC) # cost are always negative so the higher is the better for optimization
}

## NSWRMs subsidies ##

# To be use AFTER measr$implement_nswrm(nswrm_id = idx)
calc_measr_subsidies_ind <- function(measr_obj, # a measr object or the complete path to the measr file
                                     impl_nswrms_list, # = idx from SWAT.R 
                                     measr_SB  = data.frame(name = character(),
                                                             value = double())     # subsidies: a dataframe with two column: name, value
){
  # get the hru data table
  hru_con = measr_obj$.data$model_setup$original_inputs$hru.con
  # get the ids of the hru involved by the measr scenario
  nswrms_loca <- measr_obj$.data$nswrm_definition$nswrm_locations 
  # get the ids of the hru involved by the measr scenario
  tot_hru_area = 0.
  tot_SB = 0.
  
  # Loop over all measures and calculate implementation cost
  for (impl_nswrm in impl_nswrms_list) {
    # get the name of the measure
    nswrm_name = nswrms_loca$name[nswrms_loca$id == impl_nswrm]
    # get the list of hru involved by the measure
    hru_ids = unique(unlist(nswrms_loca$obj_id[nswrms_loca$id == impl_nswrm]))
    # get the areas in hectares from each hrus
    hru_area = hru_con$area[hru_con$id %in% hru_ids]
    
    sum_hru_area = sum(hru_area)
    tot_hru_area = tot_hru_area+sum_hru_area
    
    # get the cost/payment per hectare for each measure, find the best matching label
    
    # Measure maintanance costs
    SB = get_parametric_value(nswrm_name,measr_SB)
    sum_SB = sum_hru_area*SB
    tot_SB = tot_SB + sum_SB
  }
  
  return(tot_SB)
}


## OVERALL MEASURE INDICATOR ##

# To be use AFTER measr$implement_nswrm(nswrm_id = idx)
calc_measr_all_ind <- function(measr_obj, # a measr object or the complete path to the measr file
                               impl_nswrms_list, # = idx from SWAT.R
                               measr_IC = data.frame(name = character(),
                                                      value = double()),     # implementation costs: a dataframe with two column: name, value
                               measr_MC = data.frame(name = character(),
                                                      value = double()),     # maintenance costs: a dataframe with two column: name, value
                               measr_FP = data.frame(name = character(),
                                                     value = double())       # subsidies: a dataframe with two column: name, value
                              ){
  # get the hru data table
  hru_con = measr_obj$.data$model_setup$original_inputs$hru.con
  # get the ids of the hru involved by the measr scenario
  nswrms_loca <- measr_obj$.data$nswrm_definition$nswrm_locations 
  # get the ids of the hru involved by the measr scenario
  tot_hru_area = 0.
  tot_FP = 0.
  tot_MC = 0.
  tot_IC = 0.
  
  # Loop over all measures and calculate implementation cost
  for (impl_nswrm in impl_nswrms_list) {
    # get the name of the measure
    nswrm_name = nswrms_loca$name[nswrms_loca$id == impl_nswrm]
    # get the list of hru involved by the measure
    hru_ids = unique(unlist(nswrms_loca$obj_id[nswrms_loca$id == impl_nswrm]))
    # get the areas in hectares from each hrus
    hru_area = hru_con$area[hru_con$id %in% hru_ids]
    
    sum_hru_area = sum(hru_area)
    tot_hru_area = tot_hru_area+sum_hru_area
    
    # get the cost/payment per hectare for each measure, find the best matching label
    
    # Measure implementation costs
    IC = get_parametric_value(nswrm_name,measr_IC)
    sum_IC = sum_hru_area*IC
    tot_IC = tot_IC + sum_IC
    
    # Measure maintenance costs
    MC = get_parametric_value(nswrm_name,measr_MC)
    sum_MC = sum_hru_area*MC
    tot_MC = tot_MC+sum_MC
    
    # Measure farmer payments
    FP = get_parametric_value(nswrm_name,measr_FP)
    sum_FP = sum_hru_area*FP
    tot_FP = tot_FP+sum_FP
  }
  
  # make totals
  measr_all = tot_FP - tot_MC - tot_IC
  return(measr_all)
}
