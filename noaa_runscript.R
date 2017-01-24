# This is a script calling the noaa mapping function (noaa.mapping.function.R). Inputs need to be changed as needed
# and cegwas outputs need to be adjusted (e.g. making the phenotype dataframe) based on what columns are available.

library(dplyr)
library(cegwas)
library(biomaRt)

setwd("~/Dropbox/AndersenLab/LabFolders/Katie/projects/noaa/data/")

#Download data for 3 months
wi <- read.csv('~/Dropbox/AndersenLab/LabFolders/Katie/projects/noaa/data/now/20160715.wild.isolates.adjusted.csv', header=T) %>% 
  dplyr::select(-X)

wi_weather_1yr <- noaa_mappings(df = wi, time_period = 12, important_trait = "temp", additional_data = "AA1")
full_wi_weather_1yr <- wi_weather_1yr

#Remove traits and keep ones with enough data (90%)
traits_to_keep <- NULL
for(traits in colnames(wi_weather_1yr)) {
  colnum <- which(colnames(wi_weather_1yr) == traits)
  if(length(which(!is.na(wi_weather_1yr[colnum]))) > 0.9*nrow(wi_weather_1yr)) { 
    traits_to_keep <- cbind(traits_to_keep, traits)
  }
}

wi_weather_1yr <- wi_weather_1yr %>% 
  dplyr::select(one_of(traits_to_keep))

write.csv(full_wi_weather_1yr, "20161103_wiw_1yr_full.csv")
write.csv(wi_weather_1yr, "20161103_wiw_1yr.csv")

wi_weather_3yr <- noaa_mappings(df = wi, time_period = 36, important_trait = "temp", additional_data = "AA1")
full_wi_weather_3yr <- wi_weather_3yr

#Remove traits and keep ones with enough data
traits_to_keep <- NULL
for(traits in colnames(wi_weather_3yr)) {
  colnum <- which(colnames(wi_weather_3yr) == traits)
  if(length(which(!is.na(wi_weather_3yr[colnum]))) > 0.9*nrow(wi_weather_3yr)) { 
    traits_to_keep <- cbind(traits_to_keep, traits)
  }
}

wi_weather_3yr <- wi_weather_3yr %>% 
  dplyr::select(one_of(traits_to_keep))

write.csv(full_wi_weather_3yr, "20161103_wiw_3yr_full.csv")
write.csv(wi_weather_3yr, "20161103_wiw_3yr.csv")

wi_weather_1yr <- wi_weather_1yr %>% mutate(time_period = "1yr")
wi_weather_3mo <- wi_weather_3mo %>% mutate(time_period = "3mo", ceil_hgt_min = c(NA), ceil_hgt_max = c(NA), ceil_hgt_avg = c(NA))
wi_weather_3yr <- wi_weather_3yr %>% mutate(time_period = "3yr")
wi_weather <- rbind(wi_weather_3mo, wi_weather_1yr, wi_weather_3yr)


#Use cegwas to make phenotypes and mapp
phenotypes_3mo <- wi_weather_3mo %>%
  dplyr::select(-X, -isolation_date, -time_period,-ceil_hgt_avg, -ceil_hgt_min, -ceil_hgt_max, -nearest_station, -station_distance, -start_date, -end_date, -latitude, -longitude, -abslat, -elevation, -abslong) %>%
gather(trait, value, wd_min:rh_avg) %>% ####NEEDS TO BE DYNAMIC!!!#######
spread(isotype, value)

#Plot map with specific strains used in competition assays noted
processed_phenotypes_3mo <- process_pheno(phenotypes_3mo)
mapping_df_3mo <- gwas_mappings(processed_phenotypes_3mo, cores = 4)
processed_mapping_df_3mo <- process_mappings(mapping_df_3mo, phenotype_df = processed_phenotypes_3mo,
                                                 CI_size = 50, snp_grouping = 200)
#saveRDS(processed_mapping_df_3mo, "20161103_pmd_3mo_minmax.rds")

# variant_cor_3mo <- variant_correlation(processed_mapping_df_3mo, condition_trait = F)
# processed_variant_cor_3mo <- process_correlations(variant_cor_3mo)

#wi_1yr <- wi %>% filter(time_period == "1yr")

#Use cegwas to make phenotypes and mapp
phenotypes_1yr <- wi_weather_1yr %>%
  dplyr::select(-X,-time_period,  -isolation_date, -nearest_station, -station_distance, -start_date, -end_date, -latitude, -longitude, -abslat, -elevation, -abslong) %>%
  gather(trait, value, wd_min:rh_avg) %>% ####NEEDS TO BE DYNAMIC!!!#######
  spread(isotype, value)

#Plot map with specific strains used in competition assays noted
processed_phenotypes_1yr <- process_pheno(phenotypes_1yr)
mapping_df_1yr <- gwas_mappings(processed_phenotypes_1yr, cores = 4)
processed_mapping_df_1yr <- process_mappings(mapping_df_1yr, phenotype_df = processed_phenotypes_1yr,
                                             CI_size = 50, snp_grouping = 200)
#saveRDS(processed_mapping_df_1yr, "20161103_pmd_1yr_minmax.rds")

# variant_cor_1yr <- variant_correlation(processed_mapping_df_1yr, condition_trait = F)
# processed_variant_cor_1yr <- process_correlations(variant_cor_1yr)

#wi_3yr <- wi %>% filter(time_period == "3yr")

#Use cegwas to make phenotypes and mapp
phenotypes_3yr <- wi_weather_3yr %>%
  dplyr::select(-X, -time_period, -isolation_date, -nearest_station, -station_distance, -start_date, -end_date, -latitude, -longitude, -abslat, -elevation, -abslong) %>%
  gather(trait, value, wd_min:rh_avg) %>% ####NEEDS TO BE DYNAMIC!!!#######
  spread(isotype, value)

#Plot map with specific strains used in competition assays noted
processed_phenotypes_3yr <- process_pheno(phenotypes_3yr)
mapping_df_3yr <- gwas_mappings(processed_phenotypes_3yr, cores = 4)
processed_mapping_df_3yr <- process_mappings(mapping_df_3yr, phenotype_df = processed_phenotypes_3yr,
                                             CI_size = 50, snp_grouping = 200)

variant_cor_3yr <- variant_correlation(processed_mapping_df_3yr, condition_trait = F)
processed_variant_cor_3yr <- process_correlations(variant_cor_3yr)

processed_mapping_df_1yr <- processed_mapping_df_1yr %>% mutate(time_period = "1yr")
processed_mapping_df_3mo <- processed_mapping_df_3mo %>% mutate(time_period = "3mo")
processed_mapping_df_3yr <- processed_mapping_df_3yr %>% mutate(time_period = "3yr")
processed_mapping_df <- rbind(processed_mapping_df_3mo, processed_mapping_df_1yr, processed_mapping_df_3yr)

varcor <- variant_correlation(processed_mapping_df, condition_trait = F)
pvarcor <- process_correlations(varcor)

#Save R file to upload into figures file
save.image("20161103.noaa.pmd.all.RData")
