#############################################################################################
##
##~~~~~~~~~~~~~~~ Functions calculating performance indicators ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Indicator     | Description                                                   | Function(parameter)[output]                    | Files (to be defined in print.prt)
##~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Q_mean        | mean discharge [m?/s]                                         | ind_cha_aa(path, channel)[1]                   | channel_sd_aa.txt
## Nload         | total N load [kg/yr]                                          | ind_cha_aa(path, channel)[2]                   | channel_sd_aa.txt
## Pload         | total P load [kg/yr]                                          | ind_cha_aa(path, channel)[3]                   | channel_sd_aa.txt
## Sedload       | total sediment load [tons/yr]                                 | ind_cha_aa(path, channel)[4]                   | channel_sd_aa.txt
##---------------|---------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------
## Q_max         | maximum daily discharge [m?/s]                                | ind_cha_day(path, channel, 'Q_max')[1]         | channel_sd_day.txt
## Q_max_aa      | average maximum daily discharge of each year [m?/s]           | ind_cha_day(path, channel, 'Q_max')[2]         | channel_sd_day.txt
## Q_p95         | 95 percentile daily discharge [m?/s]                          | ind_cha_day(path, channel, 'Q_p95')[3]         | channel_sd_day.txt
## Q_p90         | 90 percentile daily discharge [m?/s]                          | ind_cha_day(path, channel, 'Q_p90')[4]         | channel_sd_day.txt
## Q_p50         | 50 percentile daily discharge [m?/s]                          | ind_cha_day(path, channel, 'Q_p50')[5]         | channel_sd_day.txt
## Q_p10         | 10 percentile daily discharge [m?/s]                          | ind_cha_day(path, channel, 'Q_p10')[6]         | channel_sd_day.txt
## Q_p05         | 5 percentile daily discharge [m?/s]                           | ind_cha_day(path, channel, 'Q_p05')[7]         | channel_sd_day.txt
## Q_min         | minimum daily discharge [m?/s]                                | ind_cha_day(path, channel, 'Q_min')[8]         | channel_sd_day.txt
## Q_min_aa      | average minimum daily discharge of each year  [m?/s]          | ind_cha_day(path, channel, 'Q_min')[9]         | channel_sd_day.txt
## Q_maxmin      | Q_max/Q_min ratio []                                          | ind_cha_day(path, channel, 'Q_maxmin')[8]      | channel_sd_day.txt
## Q_maxmin_aa   | Q_max_aa/Q_min_aa ratio []                                    | ind_cha_day(path, channel, 'Q_maxmin')[8]      | channel_sd_day.txt
## Q_low_days    | frequency daily discharge is below low flow threshold []      | ind_cha_day(path, channel, 'Q_low_days', threshold_lowQ)[11]       | channel_sd_day.txt
## Q_high_days   | frequency daily discharge is below high flow threshold []     | ind_cha_day(path, channel, 'Q_high_days', threshold_highQ)[12]     | channel_sd_day.txt
## Nconc_days    | frequency total N concentrations is below threshold []        | ind_cha_day(path, channel, 'Nconc_days', threshold_N)[13]          | channel_sd_day.txt
## Pconc_days    | frequency total P concentrations is below threshold []        | ind_cha_day(path, channel, 'Pconc_days', threshold_P)[14]          | channel_sd_day.txt
## Sedconc_days  | frequency total sediment concentrations is below threshold [] | ind_cha_day(path, channel, 'Sedconc_days', threshold_Sed)[15]      | channel_sd_day.txt
##---------------|---------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------
## N_loss        | average annual Nitrogen loss from land objects [kg N/ha,yr]   | ind_hru_aa_nb(path, area)[1] # area='basin' by default, area='agr' for cropland only (hru_agr.txt must be provided!) | hru_nb_aa.txt, hru_ls_aa.txt, hru_pw_aa.txt
## P_loss        | average annual Phosphorus loss from land objects [kg P/ha,yr] | ind_hru_aa_nb(path, area)[3] # area='basin' by default, area='agr' for cropland only (hru_agr.txt must be provided!) | hru_nb_aa.txt, hru_ls_aa.txt, hru_pw_aa.txt
## Sed_loss      | average annual Sediment loss from land objects [tons/ha,yr]   | ind_hru_aa_nb(path, area)[5] # area='basin' by default, area='agr' for cropland only (hru_agr.txt must be provided!) | hru_nb_aa.txt, hru_ls_aa.txt, hru_pw_aa.txt
## N_loss_ratio  | average annual Nitrogen loss/input ratio []                   | ind_hru_aa_nb(path, area)[2] # area='basin' by default, area='agr' for cropland only (hru_agr.txt must be provided!) | hru_nb_aa.txt, hru_ls_aa.txt, hru_pw_aa.txt
## P_loss_ratio  | average annual Phosphorus loss/input ratio []                 | ind_hru_aa_nb(path, area)[4] # area='basin' by default, area='agr' for cropland only (hru_agr.txt must be provided!) | hru_nb_aa.txt, hru_ls_aa.txt, hru_pw_aa.txt
##---------------|---------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------
## sw            | average annual total soil moisture [mm]                       | ind_hru_aa_wb(path, area)[1] # area='basin' by default, area='agr' for cropland only (hru_agr.txt must be provided!) | hru_wb_aa.txt
## sw300         | average annual soil moisture in top 30 cm [mm]                | ind_hru_aa_wb(path, area)[2] # area='basin' by default, area='agr' for cropland only (hru_agr.txt must be provided!) | hru_wb_aa.txt
## perc          | average annual percolation [mm]                               | ind_hru_aa_wb(path, area)[3] # area='basin' by default, area='agr' for cropland only (hru_agr.txt must be provided!) | hru_wb_aa.txt
##---------------|---------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------
## sw            | average soil moisture in period of interest [mm]              | ind_hru_mon_wb(path, sw, period, area) # e.g. period=c(5:9), area='basin' by default, area='agr' for cropland only (hru_agr.txt must be provided!) | hru_wb_mon.txt
## sw300         | average soil moisture (top 30cm) in period of interest [mm]   | ind_hru_mon_wb(path, sw300, period, area) # e.g. period=c(5:9), area='basin' by default, area='agr' for cropland only (hru_agr.txt must be provided!) | hru_wb_mon.txt
##---------------|---------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------
## grain_units   | average annual sum of grain units in whole 'basin' []         | ind_bsn_aa_crp(path, grain_units, 'grain_units', crops_sel)        | 'basin'_crop_yld_aa.txt
## cropland      | average annual area of cropland in whole 'basin' [ha]         | ind_bsn_aa_crp(path, grain_units, 'cropland', crops_sel)           | 'basin'_crop_yld_aa.txt

foo1 <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}

foo2 <- function(x){
  #  require returns TRUE invisibly if it was able to load package
  if( ! require( x , character.only = TRUE ) ){
    #  If package was not able to be loaded then re-install
    remotes::install_github("chrisschuerz/SWATfarmR")
    remotes::install_git('https://git.ufz.de/schuerz/SWATmeasR')
    #  Load package after installing
    require( x , character.only = TRUE )
  }
}

is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}

`%nin%` = Negate(`%in%`)

set2wd <- function(x) setwd(x)

# Read table function
read_tbl <- function(file, run_path, n_skip) {
  file_path <- paste0(run_path, '/', file)
  
  col_names <- read_lines(file = file_path, skip = 1, n_max = 1, lazy = FALSE) %>%
    str_trim(.) %>%
    str_split(., '[:space:]+') %>%
    unlist()
  
  name_duplicate <- table(col_names) %>%
    .[. > 1]
  if(length(name_duplicate) > 0) {
    for (i in 1:length(name_duplicate)) {
      col_names[col_names == names(name_duplicate[i])] <-
        paste0(names(name_duplicate[i]), 1:name_duplicate[i])
    }
  }
  
  fread(file_path, skip = n_skip, header = FALSE) %>%
    set_names(., col_names) %>%
    tibble(.)
}

# Indicators based on annual average channel output
ind_cha_aa <- function(path, channel){
  
  df_out <- data.frame(Q_mean=NA, 
                       Nload=NA, 
                       Pload=NA, 
                       Sedload=NA)
  # Read file
  channel_sd_aa <- read_tbl('channel_sd_aa.txt', path, 3)
  
  # Specify the columns you want to keep
  columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                       "flo_out", "sed_out", "orgn_out", "sedp_out", "no3_out", 
                       "solp_out", "nh3_out", "no2_out")
  
  # Create a new data frame with only the selected columns
  df_selected <- channel_sd_aa[, columns_to_keep]
  
  df_selected <- df_selected %>%
    mutate(total_N = no3_out + orgn_out + nh3_out + no2_out) %>% 
    mutate(total_P = solp_out + sedp_out)
  
  Q_mean <- round(df_selected$flo_out[which(df_selected$name==channel)],3)
  Nload <- round(df_selected$total_N[which(df_selected$name==channel)],3)
  Pload <- round(df_selected$total_P[which(df_selected$name==channel)],3)
  Sedload <- round(df_selected$sed_out[which(df_selected$name==channel)],3)
  
  df_out[1,1:4] <- c(Q_mean, Nload, Pload, Sedload)
  
  return(df_out)
}

# Indicators based on daily channel output in cha_day.out (to be defined in object.prt file,
# see example file, you have to specify this file in file.cio)
# !!! And don't forget to deactivate daily channel_sd printing in print.prt !!!
ind_cha_day <- function(path,
                        channel, 
                        ind='all', 
                        threshold_lowQ=0.0344,
                        threshold_highQ=2.7911,
                        threshold_N=2.3, 
                        threshold_P=0.082, 
                        threshold_Sed=50){
  
  df_out <- data.frame(Q_max=NA,
                       Q_max_aa = NA,
                       Q_p95=NA, 
                       Q_p90=NA,
                       Q_p50=NA, 
                       Q_p10=NA, 
                       Q_p05=NA, 
                       Q_min=NA,
                       Q_min_aa = NA,
                       Q_maxmin=NA,
                       Q_maxmin_aa = NA,
                       Q_low_days=NA,
                       Q_high_days=NA,
                       Nconc_days=NA, 
                       Pconc_days=NA,
                       Sedconc_days=NA)
  
  channel_sd_day <- read.table(paste0(path,'/cha_day.out'), h=T)
  names(channel_sd_day)[c(5,6)] <- c('type','name')
  
  if('Q_maxmin' %in% ind | 'Q_max' %in% ind  | 'Q_min' %in% ind | 'all' %in% ind){
    
    # convert m?/day to m?/s
    channel_sd_day$flo <- channel_sd_day$flo/86400
    
    # Handle zero flow (defining 0.00000001 cm?/s)
    channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
    
    # Group data by channel (replace "name" with the actual column name for channel)
    max_min_ratio <- channel_sd_day %>%
      group_by(name) %>% 
      summarise(
        max_discharge = max(flo, na.rm = TRUE),
        min_discharge = min(flo, na.rm = TRUE),
        extreme_streamflow_ratio = max_discharge / min_discharge
      )
    
    max_min_ratio_aa <- channel_sd_day %>%
      group_by(name,yr) %>% 
      summarise(
        max_discharge_yr = max(flo, na.rm = TRUE),
        min_discharge_yr = min(flo, na.rm = TRUE),
      ) %>% 
      group_by(name) %>% 
      summarise(
        max_discharge_aa = mean(max_discharge_yr, na.rm = TRUE),
        min_discharge_aa = mean(min_discharge_yr, na.rm = TRUE),
        extreme_streamflow_ratio_aa = max_discharge_aa / min_discharge_aa
      )
    
    df_out[1,1] <- round(max_min_ratio$max_discharge[which(max_min_ratio$name==channel)],5)
    df_out[1,2] <- round(max_min_ratio_aa$max_discharge_aa[which(max_min_ratio_aa$name==channel)],5)
    df_out[1,8] <- round(max_min_ratio$min_discharge[which(max_min_ratio$name==channel)],5)
    df_out[1,9] <- round(max_min_ratio_aa$min_discharge_aa[which(max_min_ratio_aa$name==channel)],5)
    df_out[1,10] <- round(max_min_ratio$extreme_streamflow_ratio[which(max_min_ratio$name==channel)],5)
    df_out[1,11] <- round(max_min_ratio_aa$extreme_streamflow_ratio_aa[which(max_min_ratio_aa$name==channel)],5)
  }
  if('Q_p50' %in% ind | 'all' %in% ind){
    
    # convert m?/day to m?/s
    if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
    
    # Handle zero flow (defining 0.00000001 cm?/s)
    channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
    
    # Group data by channel (replace "name" with the actual column name for channel)
    p50 <- channel_sd_day %>%
      group_by(name) %>% 
      summarise(
        p50_discharge = quantile(flo, probs = 0.50, na.rm = TRUE)
      )
    
    df_out[1,5] <- round(as.numeric(p50$p50_discharge[which(p50$name==channel)]),5)
  }
  if('Q_p95p05' %in% ind | 'Q_p95' %in% ind | 'Q_p05' %in% ind | 'all' %in% ind){
    
    # convert m?/day to m?/s
    if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
    
    # Handle zero flow (defining 0.00000001 cm?/s)
    channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
    
    # Group data by channel (replace "name" with the actual column name for channel)
    Q_p95p05 <- channel_sd_day %>%
      group_by(name) %>%
      summarise(
        p05_discharge = quantile(flo, probs = 0.05, na.rm = TRUE),
        p95_discharge = quantile(flo, probs = 0.95, na.rm = TRUE),
      )
    df_out[1,3] <- round(as.numeric(Q_p95p05$p95_discharge[which(Q_p95p05$name==channel)]),5)
    df_out[1,7] <- round(as.numeric(Q_p95p05$p05_discharge[which(Q_p95p05$name==channel)]),5)
  }
  if('Q_p90p10' %in% ind | 'Q_p90' %in% ind | 'Q_p10' %in% ind | 'all' %in% ind){
    
    # convert m?/day to m?/s
    if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
    
    # Handle zero flow (defining 0.00000001 cm?/s)
    channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
    
    # Group data by channel (replace "name" with the actual column name for channel)
    Q_p90p10 <- channel_sd_day %>%
      group_by(name) %>%
      summarise(
        p10_discharge = quantile(flo, probs = 0.10, na.rm = TRUE),
        p90_discharge = quantile(flo, probs = 0.90, na.rm = TRUE),
      )
    df_out[1,4] <- round(as.numeric(Q_p90p10$p90_discharge[which(Q_p90p10$name==channel)]),5)
    df_out[1,6] <- round(as.numeric(Q_p90p10$p10_discharge[which(Q_p90p10$name==channel)]),5)
  }
  if('Q_low_days' %in% ind | 'Q_high_days' %in% ind | 'Nconc_days' %in% ind | 'Pconc_days' %in% ind | 'Sedconc_days' %in% ind | 'all' %in% ind){
    
    # convert m?/day to m?/s
    if('all' %in% ind == F) channel_sd_day$flo <- channel_sd_day$flo/86400
    
    # Handle zero flow (defining 0.00000001 cm?/s)
    channel_sd_day$flo[which(channel_sd_day$flo==0)] <- 1e-40
    
    # Calculate the sum of no3, orgn, nh3, and no2
    channel_sd_day <- channel_sd_day %>%
      mutate(total_N = no3 + orgn + nh3 + no2) %>% 
      mutate(N_conc_mgl = ifelse(flo == 0, 0, (total_N * 1000) / (flo * 86400))) %>% 
      mutate(total_P = solp + sedp) %>% 
      mutate(P_conc_mgl = ifelse(flo == 0, 0, (total_P * 1000) / (flo * 86400))) %>% 
      mutate(sed_conc_mgl = ifelse(flo == 0, 0, (sed * 1e6) / (flo * 86400)))
    
    # Calculate the frequency of exceeding the thresholds for each "unit" (channel)
    frequency_summary_mean <- channel_sd_day %>%
      group_by(name) %>%
      summarize(
        freq_below_threshold_lowQ = mean(flo <= threshold_lowQ, na.rm = TRUE),
        freq_beyond_threshold_highQ = mean(flo > threshold_highQ, na.rm = TRUE),
        freq_beyond_threshold_N = mean(N_conc_mgl > threshold_N, na.rm = TRUE),
        freq_beyond_threshold_P = mean(P_conc_mgl > threshold_P, na.rm = TRUE),
        freq_beyond_threshold_Sed = mean(sed_conc_mgl > threshold_Sed, na.rm = TRUE)
      )
    df_out[1,14] <- round(as.numeric(frequency_summary_mean$freq_beyond_threshold_N[which(frequency_summary_mean$name==channel)]),5)
    df_out[1,15] <- round(as.numeric(frequency_summary_mean$freq_beyond_threshold_P[which(frequency_summary_mean$name==channel)]),5)
    df_out[1,16] <- round(as.numeric(frequency_summary_mean$freq_beyond_threshold_Sed[which(frequency_summary_mean$name==channel)]),5)
    df_out[1,12] <- round(as.numeric(frequency_summary_mean$freq_below_threshold_lowQ[which(frequency_summary_mean$name==channel)]),5)
    df_out[1,13] <- round(as.numeric(frequency_summary_mean$freq_beyond_threshold_highQ[which(frequency_summary_mean$name==channel)]),5)
  }
  
  return(df_out)
}
                          

# water balance related indicators based on annual average hru output
ind_hru_aa_wb <- function(path, area = 'basin'){
  
  df_out <- data.frame(sw=NA,
                       sw300=NA,
                       perc=NA)
  
  # Read file
  hru_wb <- read_tbl('hru_wb_aa.txt', path, 3)
  hru_area <- read_tbl ('hru.con', path, 2)
  
  # Specify the columns you want to keep
  # Keep all hru_ls columns
  
  columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                       "sw_ave", "sw_300","perc")
  
  # Create a new data frame with only the selected columns
  if(area == 'basin'){
    df_selected_hru_wb <- hru_wb[, columns_to_keep]
    idx <- c(1:dim(hru_wb)[1])
  }else{
    # Read in vector for agricultural area
    hru_agr <- read.table(paste0(path[i],'/hru_agr.txt'), h=T)
    idx <- hru_agr$hru_id
    df_selected_hru_wb <- hru_wb[idx, columns_to_keep]
  }
  
  # Add the HRU area to each HRU
  hru_wb <- hru_wb[idx,] %>% left_join(hru_area[idx,] %>% select(id, area), by = c("unit" = "id"))
  
  # Calculate the weighted values
  df_selected_hru_wb  <- df_selected_hru_wb  %>%
    mutate(weighted_sw = sw_ave * hru_wb$area,
           weighted_sw300 = sw_300 * hru_wb$area,
           weighted_perc = perc * hru_wb$area)
  
  # Calculate the total weighted sum
  total_weighted_sw <- sum(df_selected_hru_wb$weighted_sw, na.rm = TRUE)
  total_weighted_sw300 <- sum(df_selected_hru_wb$weighted_sw300, na.rm = TRUE)
  total_weighted_perc <- sum(df_selected_hru_wb$weighted_perc, na.rm = TRUE)
  
  # Calculate the total area across all HRUs
  total_area <- sum(hru_wb$area, na.rm = TRUE)
  
  # Calculate the area-weighted averages
  df_out[1,1] <- round(total_weighted_sw / total_area,3)
  df_out[1,2] <- round(total_weighted_sw300 / total_area,3)
  df_out[1,3] <- round(total_weighted_perc / total_area,3)
  
  return(df_out)
}

# nutrient and sediment related indicators based on annual average hru output
ind_hru_aa_nb <- function(path, area = 'basin'){
  
  df_out <- data.frame(N_loss=NA, 
                       P_loss=NA,
                       Sed_loss=NA,
                       N_loss_ratio=NA, 
                       P_loss_ratio=NA)
  
  # Read file
  hru_ls <- read_tbl('hru_ls_aa.txt', path, 3)
  hru_nb <- read_tbl('hru_nb_aa.txt', path, 3)
  hru_pw <- read_tbl('hru_pw_aa.txt', path, 3)
  hru_area <- read_tbl ('hru.con', path, 2)
  
  # Specify the columns you want to keep
  # Keep all hru_ls columns
  
  columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                       "fertn", "fixn", "no3atmo", "nh4atmo", "fertp", "denit")
  
  
  # Create a new data frame with only the selected columns
  if(area == 'basin'){
    df_selected_hru_nb <- hru_nb[, columns_to_keep]
    idx <- c(1:dim(hru_nb)[1])
  }else{
    # Read in vector for agricultural area
    hru_agr <- read.table(paste0(path[i],'/hru_agr.txt'), h=T)
    idx <- hru_agr$hru_id
    df_selected_hru_nb <- hru_nb[idx, columns_to_keep]
  }
  
  # Calculate the sum of N inputs
  df_selected_hru_nb <- df_selected_hru_nb %>%
    mutate(N_inputs = fertn + fixn + no3atmo + nh4atmo)
  
  # Calculate N losses 
  # Add N_losses as a new column to df_selected_hru_nb
  
  hru_ls <- hru_ls[idx,] %>%
    mutate(N_losses = sedorgn + surqno3 + lat3no3 + tileno3 + hru_pw$percn[idx])
  
  # The P input is only fertp from hru_nb
  # Calculate the sum of P losses
  hru_ls <- hru_ls %>%
    mutate(P_losses = sedorgp + surqsolp + sedminp + tilelabp + lchlabp)
  
  # Add the HRU area to each HRU
  hru_ls <- hru_ls %>% left_join(hru_area[idx,] %>% select(id, area), by = c("unit" = "id"))
  
  # Calculate the weighted values for N and P inputs
  df_selected_hru_nb  <- df_selected_hru_nb  %>%
    mutate(weighted_N_inputs = N_inputs * hru_ls$area,
           weighted_P_inputs = fertp * hru_ls$area)
  # Calculate the weighted values for N, P and sediment losses
  hru_ls <- hru_ls %>%
    mutate(weighted_N_losses = N_losses * area,
           weighted_P_losses = P_losses * area,
           weighted_Sed_losses = sedyld * area)
  
  # Calculate the total weighted sum for N and P
  total_weighted_N_inputs <- sum(df_selected_hru_nb$weighted_N_inputs, na.rm = TRUE)
  total_weighted_N_losses <- sum(hru_ls$weighted_N_losses, na.rm = TRUE)
  total_weighted_P_inputs <- sum(df_selected_hru_nb$weighted_P_inputs, na.rm = TRUE)
  total_weighted_P_losses <- sum(hru_ls$weighted_P_losses, na.rm = TRUE)
  total_weighted_Sed_losses <- sum(hru_ls$weighted_Sed_losses, na.rm = TRUE)
  
  # Calculate the total area across all HRUs
  total_area <- sum(hru_ls$area, na.rm = TRUE)
  
  # Calculate the area-weighted averages for N and P
  area_weighted_average_N_inputs <- total_weighted_N_inputs / total_area
  area_weighted_average_P_inputs <- total_weighted_P_inputs / total_area
  area_weighted_average_N_losses <- total_weighted_N_losses / total_area
  area_weighted_average_P_losses <- total_weighted_P_losses / total_area
  area_weighted_average_Sed_losses <- total_weighted_Sed_losses / total_area
  
  df_out[1,1] <- round(area_weighted_average_N_losses,3)
  df_out[1,2] <- round(area_weighted_average_P_losses,3)
  df_out[1,3] <- round(area_weighted_average_Sed_losses,3)
  df_out[1,4] <- round(area_weighted_average_N_losses/area_weighted_average_N_inputs,3)
  df_out[1,5] <- round(area_weighted_average_P_losses/area_weighted_average_P_inputs,3)
  
  return(df_out)
}

# water balance related indicators based on monthly hru outputs
ind_hru_mon_wb <- function(path, ind = 'sw', period = c(5:9), area = 'basin'){
  
  if(ind == 'sw'){
    sw_colnames <- paste0('sw_', map_chr(period, ~paste(.x, collapse = "_")))
    df_out <- data.frame(matrix(data=NA, nrow=1, ncol=length(sw_colnames)))
    names(df_out) <- sw_colnames
    
    # Read file
    hru_wb <- read_tbl('hru_wb_mon.txt', path, 3)
    hru_area <- read_tbl ('hru.con', path, 2)
    
    # Specify the columns you want to keep for hru_nb
    # Keep all hru_ls columns
    
    columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                         "sw_ave")
    
    for(k in 1:length(period)){
      # Create a new data frame with only the selected rows and columns
      if(area == 'basin'){
        df_selected_hru_wb <- hru_wb[which(hru_wb$mon %in% period[[k]]), columns_to_keep]
        idx <- c(1:dim(hru_area)[1])
      }else{
        # Read in vector for agricultural area
        hru_agr <- read.table(paste0(path,'/hru_agr.txt'), h=T)
        idx <- hru_agr$hru_id
        df_selected_hru_wb <- hru_wb[which(hru_wb$unit %in% idx & hru_wb$mon %in% period[[k]]), columns_to_keep]
      }
      
      # Calculate the area weighted values
      df_selected_hru_wb  <- df_selected_hru_wb  %>%
        left_join(hru_area %>% select(id, area), by = c("unit" = "id")) %>% 
        mutate(weighted_sw = sw_ave * area)
      
      # Number of months and years
      n_mon <- length(period[[k]])
      n_years <- length(unique(df_selected_hru_wb$yr))
      
      # Calculate the total weighted sum
      total_weighted_sw <- sum(df_selected_hru_wb$weighted_sw, na.rm = TRUE)/n_mon/n_years
      
      # Calculate the total area across all HRUs
      total_area <- sum(hru_area$area[idx], na.rm = TRUE)
      
      # Calculate the area-weighted averages for N and P
      df_out[1,k] <- round(total_weighted_sw / total_area,3)
    }
  }
  
  if(ind == 'sw300'){
    sw_colnames <- paste0('sw300_', map_chr(period, ~paste(.x, collapse = "_")))
    df_out <- data.frame(matrix(data=NA, nrow=1, ncol=length(sw_colnames)))
    names(df_out) <- sw_colnames
    
    # Read file
    hru_wb <- read_tbl('hru_wb_mon.txt', path, 3)
    hru_area <- read_tbl ('hru.con', path, 2)
    
    # Specify the columns you want to keep for hru_nb
    # Keep all hru_ls columns
    
    columns_to_keep <- c("jday", "mon", "day", "yr", "unit", "gis_id", "name", 
                         "sw_300")
    
    for(k in 1:length(period)){
      # Create a new data frame with only the selected rows and columns
      if(area == 'basin'){
        df_selected_hru_wb <- hru_wb[which(hru_wb$mon %in% period[[k]]), columns_to_keep]
        idx <- c(1:dim(hru_area)[1])
      }else{
        # Read in vector for agricultural area
        hru_agr <- read.table(paste0(path,'/hru_agr.txt'), h=T)
        idx <- hru_agr$hru_id
        df_selected_hru_wb <- hru_wb[which(hru_wb$unit %in% idx & hru_wb$mon %in% period[[k]]), columns_to_keep]
      }
      
      # Calculate the area weighted values
      df_selected_hru_wb  <- df_selected_hru_wb  %>%
        left_join(hru_area %>% select(id, area), by = c("unit" = "id")) %>% 
        mutate(weighted_sw = sw_300 * area)
      
      # Number of months and years
      n_mon <- length(period[[k]])
      n_years <- length(unique(df_selected_hru_wb$yr))
      
      # Calculate the total weighted sum
      total_weighted_sw <- sum(df_selected_hru_wb$weighted_sw, na.rm = TRUE)/n_mon/n_years
      
      # Calculate the total area across all HRUs
      total_area <- sum(hru_area$area[idx], na.rm = TRUE)
      
      # Calculate the area-weighted averages for N and P
      df_out[1,k] <- round(total_weighted_sw / total_area,3)
    }
  }
  
  return(df_out)
}

# crop yield related indicators (grain units) based on annual average hru output
ind_bsn_aa_crp <- function(path, crop_sel, ind, grain_units){
  
  if (ind == "grain_units") {
    df_out <- data.frame(grain_units=NA)
  } else {
    df_out <- data.frame(crops_ha=NA)
  }
  
  # read output file
  crop_aa <- read_tbl('basin_crop_yld_aa.txt', path, 2)
  # Index for reading
  if (ind == "grain_units") { 
    crop_sel <- names(grain_units)[match(crop_sel, names(grain_units))]
    idx <- pmatch(crop_sel, crop_aa$plant_name, duplicates.ok = T)
    # Convert yield to grain units and sum it up
    crop_yld_gu <- round(sum(crop_aa$`yld(t)`[idx] * grain_units, na.rm = T),3)
    # collect grain unit values in df_out
    df_out[1,1] <- crop_yld_gu
  } else {
    idx <- pmatch(crop_sel, crop_aa$plant_name, duplicates.ok = T)
    # sum up hectare values and collect in df_out
    df_out[1,1] <- round(sum(crop_aa$`harv_area(ha)`[idx], na.rm = T),2)
  }
  
  return(df_out)
}

