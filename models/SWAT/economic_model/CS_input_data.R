# spi parameters
# origin: F:/mstrauch/CoMOLA_OPTAIN_CS1/models/SWAT/economic_model/CS_input_data.xlsx

crop_prices_df = data.frame(
"name" = c("rnge", "wwht", "wira", "wbar", "wiry", "akgs", "csil", "barl", "radi", "sgbt", "frst", "fesc", "orcd", "rngb", "wetl", "frsd"),
"value" = c(0, 242.597, 2435, 205.275, 196.078, 30.303, 117.284, 255.734, 0, 152.174, 0, 81.369, 4199.5, 0, 0, 0)
)

crop_farm_payms_df = data.frame(
"name" = c("rnge", "wwht", "wira", "wbar", "wiry", "akgs", "csil", "barl", "radi", "sgbt", "frst", "fesc", "orcd", "rngb", "wetl", "frsd"),
"value" = c(540, 0, 0, 0, 0, 0, 0, 0, 78, 0, 0, 0, 0, 0, 0, 0)
)

crop_prod_costs_df = data.frame(
"name" = c("rnge", "wwht", "wira", "wbar", "wiry", "akgs", "csil", "barl", "radi", "sgbt", "frst", "fesc", "orcd", "rngb", "wetl", "frsd"),
"value" = c(524, 891, 898, 799, 750, 1159, 1064, 646, 847, 1317, 136.67, 524, 2446.9, 524, 524, 137)
)

crop_prod_costs_mgt_df = data.frame(
"op_typ" = c("fert", "fert", "fert", "harv", "kill", "till", "fert", "till", "till", "plnt", "skip", "fert", "fert", "harv", "kill", "fert", "till", "till", "plnt", "fert", "till", "plnt", "fert", "fert", "harv", "kill", "till", "fert", "fert", "harv", "kill", "till", "plnt", "harv", "fert", "kill", "plnt", "fert", "fert", "fert", "fert", "till", "plnt", "harv", "kill", "fert", "fert", "plnt", "harv", "kill", "fert", "fert", "fert", "fert", "fert", "fert", "fert", "fert", "fert", "fert", "plnt", "fert", "harv", 
"kill", "harv", "kill", "plnt", "fert", "fert", "fert", "fert", "fert", "fert", "fert", "plnt", "harv", "kill", "fert", "fert", "fert", "fert", "fert", "fert", "fert", "harv", "fert", "fert", "fert", "harv"),
"op_data1" = c("elem_n", "elem_n", "elem_n", "wwht", "wwht", "fldcul10", "elem_p", "cultiv20", "harrow5", "wira", "null", "elem_n", "elem_n", "wira", "wira", "elem_p", "cultiv25", "harrow7", "wwht", "elem_p", "fldcul12", "wbar", "elem_n", "elem_n", "wbar", "wbar", "fldcul15", "beefg_fl", "elem_n", "wiry", "wiry", "cultiv30", "akgs", "akgs", "elem_n", "akgs", "wiry", "beefg_fl", "elem_n", "elem_n", "elem_p", "harrow8", "csil", "csil", "csil", "beefg_fl", "elem_n", "radi", "radi", "radi", "beefg_fl", "elem_n", 
"elem_p", "beefg_fl", "beefg_fl", "elem_n", "elem_n", "elem_n", "elem_n", "beefg_fl", "barl", "elem_n", "barl", "barl", "rnge_test", "rnge_test", "rnge_test", "elem_p", "elem_n", "elem_n", "elem_p", "elem_n", "elem_p", "elem_n", "sgbt", "sgbt", "sgbt", "beefg_fs", "elem_n", "elem_n", "beefg_fs", "elem_n", "elem_n", "elem_p", "fesc_cs1", "elem_n", "elem_n", "beefg_fl", "orcd"),
"op_data2" = c("broadcast", "broadcast", "broadcast", "grain", "null", "null", "broadcast", "null", "null", "null", "null", "broadcast", "broadcast", "grain", "null", "broadcast", "null", "null", "null", "broadcast", "null", "null", "broadcast", "broadcast", "grain", "null", "null", "aerial_liquid", "broadcast", "grain", "null", "null", "null", "hay_cut_low", "broadcast", "null", "null", "aerial_liquid", "broadcast", "broadcast", "broadcast", "null", "null", "silage", "null", "aerial_liquid", "broadcast", "null", 
"grass_mulch", "null", "aerial_liquid", "broadcast", "broadcast", "aerial_liquid", "aerial_liquid", "broadcast", "broadcast", "broadcast", "broadcast", "aerial_liquid", "null", "broadcast", "grain", "null", "grass_mulch", "null", "null", "broadcast", "broadcast", "broadcast", "broadcast", "broadcast", "broadcast", "broadcast", "null", "tuber", "null", "broadcast", "broadcast", "broadcast", "broadcast", "broadcast", "broadcast", "broadcast", "hay_cut_low", "broadcast", "broadcast", "aerial_liquid", 
"orchard"),
"op_data3" = c(85.47, 55, 16.5, 0, 0, 0, 27.8, 0, 0, 0, 0, 105, 59.01, 0, 0, 25.4, 0, 0, 0, 22.9, 0, 0, 71.82, 42, 0, 0, 0, 22000, 43417, 0, 0, 0, 0, 0, 110, 0, 0, 27500, 88.99, 132, 32.9, 0, 0, 0, 0, 21000, 68145, 0, 0, 0, 40000, 15, 15, 15750, 15000, 95, 44.3, 88, 77, 20000, 0, 40.25, 0, 0, 0, 0, 0, 20.2, 64995, 21, 16.7, 62.7, 45.6, 93, 0, 0, 0, 20000, 54.58, 102.3, 22000, 60038, 60, 25, 0, 40, 70, 25000, 0),
"value" = c(NA, NA, NA, 891, NA, -58, NA, -58, -58, NA, NA, NA, NA, 898, NA, NA, -58, -58, NA, NA, -58, NA, NA, NA, 799, NA, -58, NA, NA, 750, NA, -58, NA, 1159, NA, NA, NA, NA, NA, NA, NA, -58, NA, 1064, NA, NA, NA, NA, 847, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 646, NA, 524, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1317, NA, NA, NA, NA, NA, NA, NA, NA, 524, NA, NA, NA, 2447)
)

msr_impl_costs_df = data.frame(
"name" = c("buffer", "hedge", "grassslope", "pond", "lowtillcc"),
"value" = c(0, 2500, 27, 7407, 134)
)

msr_mnt_costs_df = data.frame(
"name" = c("buffer", "hedge", "grassslope", "pond", "lowtillcc"),
"value" = c(-896, 22539, -896, 26545, -53)
)

msr_subsidies_df = data.frame(
"name" = c("buffer", "hedge", "grassslope", "contr_drn", "pond", "lowtillcc"),
"value" = c(0, 0, 299, NA, 0, 80)
)

env_prices_df = data.frame(
"name" = c("N", "sed", "C", "P"),
"value" = c(14, 12, 86.53, 15.67)
)

c_seq_performance_df = data.frame(
"name" = c("field", "rnge", "fld", "frst", "fd", "mdw_cts", "orcd", "rngb", "mw_ts_drn", "mdw_cts_drn", "utrn", "urld", "wetl", "urmd", "bsvg", "buffer", "grassslope", "hedge"),
"value" = c(0.1, 0.9, 0.1, 0.9, 0.1, 0.9, 0.9, 0.9, 0.9, 0.9, 0, 0, 0.9, 0, 0.9, 0.9, 0.9, 0.9)
)

TIC_df = data.frame(
"name" = logical(0),
"value" = logical(0)
)

ST_df = data.frame(
"name" = logical(0),
"value" = logical(0)
)
