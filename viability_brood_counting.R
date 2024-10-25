## ----setup, include=FALSE---------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---------------------------------------------------------------------------------
library(ggplot2); library(dplyr); library(reshape); library(Hmisc); library(tidyverse); library(stats);library(lubridate);library(plotrix); library(sjPlot); library(cowplot); library(svglite)

## ---------------------------------------------------------------------------------
get_file_info = function(trailingOnly = TRUE, asValues = TRUE) {
  data_file_info = commandArgs(TRUE)
  if (length(data_file_info) < 3){
    stop("Provide the file path to functions,
         data file dir, and filename, eg. dev/ centriole_data_dir/ centriole_data.ext")
  } else {
    return(data_file_info)
  }
}

viability_and_brood_file_info = get_file_info()
viability_and_brood_datafile_directory = viability_and_brood_file_info[1]
viability_data_filename = viability_and_brood_file_info[2]
brood_data_filename = viability_and_brood_file_info[3]
## ---------------------------------------------------------------------------------

viability_dataframe = read.csv(paste0(viability_and_brood_datafile_directory,viability_data_filename))
brood_dataframe = read.csv(paste0(viability_and_brood_datafile_directory,brood_data_filename))

## ---------------------------------------------------------------------------------
by_species = viability_dataframe %>% group_by(strain)
brood_by_species = brood_dataframe %>% group_by(strain)

by_species_no_null_values = drop_na(by_species)
brood_by_species_no_null_values = drop_na(brood_by_species)

by_species_no_null_values = by_species_no_null_values[, -c(2)]
brood_by_species_no_null_values = brood_by_species_no_null_values[, -c(2)]



## ---------------------------------------------------------------------------------
avg_strain_viability = by_species_no_null_values %>% summarise(
                          day_1 = mean(day_1),
                          day_2 = mean(day_2),
                          day_3 = mean(day_3))
std_strain_viability = by_species_no_null_values %>% summarise(
                          day_1_sd = std.error(day_1),
                          day_2_sd = std.error(day_2),
                          day_3_sd = std.error(day_3))
avg_strain_brood = brood_by_species_no_null_values %>% summarise(
                          day_1 = mean(day_1),
                          day_2 = mean(day_2),
                          day_3 = mean(day_3))
std_strain_brood = brood_by_species_no_null_values %>% summarise(
                          day_1_sd = std.error(day_1),
                          day_2_sd = std.error(day_2),
                          day_3_sd = std.error(day_3))


## ---------------------------------------------------------------------------------
# reshpae avg viability and brood size to long
strain_viability_long = avg_strain_viability %>% pivot_longer(cols = c( 'day_1','day_2', 'day_3'), names_to = 'day_measured', values_to = 'viability')
strain_brood_long = avg_strain_brood %>% pivot_longer(cols = c('day_1','day_2', 'day_3'), names_to = 'day_measured', values_to = 'brood_size')

# reshape std viability and brood size to long
strain_viability_sd_long = std_strain_viability %>% pivot_longer(cols = c( 'day_1_sd','day_2_sd', 'day_3_sd'), names_to = 'day_measured', values_to = 'viability_sd')
strain_brood_sd_long = std_strain_brood %>% pivot_longer(cols = c( 'day_1_sd','day_2_sd', 'day_3_sd'), names_to = 'day_measured', values_to = 'brood_size_sd')


## ---------------------------------------------------------------------------------
strain_viability_long["date_time"] = rep(ymd_hms("2022-9-14-13-00-00"),length(strain_viability_long["day_measured"]))

strain_viability_long["date_time"][strain_viability_long["day_measured"] == 'day_1'] = ymd_hms("2023-11-15-17-30-00")
strain_viability_long["date_time"][strain_viability_long["day_measured"] == 'day_2'] = ymd_hms("2023-11-16-14-15-00")
strain_viability_long["date_time"][strain_viability_long["day_measured"] == 'day_3'] = ymd_hms("2023-11-17-11-30-00")


## ---------------------------------------------------------------------------------
strain_brood_long["date_time"] = rep(ymd_hms("2022-9-14-13-00-00"),length(strain_brood_long["day_measured"]))

strain_brood_long["date_time"][strain_brood_long["day_measured"] == 'day_1'] = ymd_hms("2023-11-15-17-30-00")
strain_brood_long["date_time"][strain_brood_long["day_measured"] == 'day_2'] = ymd_hms("2023-11-16-14-15-00")
strain_brood_long["date_time"][strain_brood_long["day_measured"] == 'day_3'] = ymd_hms("2023-11-17-11-30-00")


## ---------------------------------------------------------------------------------
strain_viability_long$strain = factor(strain_viability_long$strain, levels = unique(strain_viability_long$strain))

strain_brood_long$strain = factor(strain_brood_long$strain, levels = unique(strain_brood_long$strain))


## ---------------------------------------------------------------------------------
viability_bar_plot = ggplot(strain_viability_long, mapping = aes(date_time, viability, fill = strain))

viability_bar_plot + geom_bar(stat = "identity", position=position_dodge()) + theme_cowplot(16) + geom_errorbar(position = position_dodge(), color = "black", aes(ymin = viability - strain_viability_sd_long$viability_sd, ymax = viability + strain_viability_sd_long$viability_sd))


## ---------------------------------------------------------------------------------
ggsave(paste0(viability_and_brood_datafile_directory,substring(viability_data_filename,1,7),"_viability_barplot.pdf"), width = 20, height = 20, units = "cm", device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
brood_size_bar_plot = ggplot(strain_brood_long, mapping = aes(date_time, brood_size, fill = strain))

brood_size_bar_plot + geom_bar(stat = "identity", position=position_dodge()) + theme_cowplot(16) + geom_errorbar(position = position_dodge(), color = "black", aes(ymin = brood_size - strain_brood_sd_long$brood_size_sd, ymax = brood_size + strain_brood_sd_long$brood_size_sd))

## ---------------------------------------------------------------------------------
ggsave(paste0(viability_and_brood_datafile_directory,substring(brood_data_filename,1,7),"_brood_size_barplot.pdf"), width = 20, height = 20, units = "cm", device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
strain_viability_plot = ggplot(strain_viability_long, mapping = aes(date_time, viability, color = strain))

strain_brood_plot = ggplot(strain_brood_long,mapping = aes(date_time, brood_size, color = strain))


## ---------------------------------------------------------------------------------
viability_plot = (strain_viability_plot + geom_line(size = 2) + geom_point(aes(shape = strain), size = 3) + geom_errorbar(aes(ymin = viability - strain_viability_sd_long$viability_sd, ymax = viability + strain_viability_sd_long$viability_sd), width = 2) + scale_color_brewer(palette = "RdPu") + theme_cowplot(14))
viability_plot

# to change tic marks: scale_x_continuous(breaks = round(seq(min(dat$x), max(dat$x), by = 0.5),1))


## ---------------------------------------------------------------------------------
ggsave(paste0(viability_and_brood_datafile_directory,substring(viability_data_filename,1,7),"_viability_lineplot.pdf"), plot = viability_plot, width = 20, height = 20, units = "cm", device = "pdf", dpi = 300)


## ---------------------------------------------------------------------------------
brood_plot = (strain_brood_plot + geom_line(size = 2) + geom_point(aes(shape = strain), size = 3) + geom_errorbar(aes(ymin = brood_size - strain_brood_sd_long$brood_size_sd, ymax = brood_size + strain_brood_sd_long$brood_size_sd), width = 3) + scale_color_brewer(palette = "RdPu") + theme_cowplot(14))
brood_plot + scale_y_continuous(limits = c(0,150))


## ---------------------------------------------------------------------------------
ggsave(paste0(viability_and_brood_datafile_directory,substring(brood_data_filename,1,7),"_brood_size_lineplot.pdf"), width = 20, height = 20, units = "cm", device = "pdf", dpi = 300)


## cumulative brood size plot---------------------------------------------------------------------------------
brood_dataframe = drop_na(brood_dataframe)
cumulative_brood_size = brood_dataframe %>%
  mutate(cumulative_brood_size = rowSums(across(c(day_0, day_1, day_2, day_3))))


## ---------------------------------------------------------------------------------
cumulative_species_brood_size = cumulative_brood_size %>% group_by(strain) %>% summarise(
                          Sum = sum(cumulative_brood_size))
avg_strain_brood = cumulative_brood_size %>% group_by(strain) %>% summarise(
                          Mean = mean(cumulative_brood_size), std = std.error(cumulative_brood_size),
                          N = n())
alpha = 0.05
avg_strain_brood$degrees_freedom = avg_strain_brood$N - 1
avg_strain_brood$t_score = qt(p=alpha/2, df=avg_strain_brood$degrees_freedom,lower.tail=F)
avg_strain_brood$margin_error <- avg_strain_brood$t_score * avg_strain_brood$std
avg_strain_brood$lower_bound <- avg_strain_brood$Mean - avg_strain_brood$margin_error
avg_strain_brood$upper_bound <- avg_strain_brood$Mean + avg_strain_brood$margin_error


## ---------------------------------------------------------------------------------
viability_dataframe = drop_na(viability_dataframe)
cumulative_viability = viability_dataframe %>%
  mutate(cumulative_viability = rowSums(across(c(day_0, day_1, day_2, day_3))))
cumulative_viability$avg_viability = cumulative_viability$cumulative_viability/3
cumulative_viability$embryo_number = cumulative_brood_size$cumulative_brood_size


## ---------------------------------------------------------------------------------
cumulative_species_viability = cumulative_viability %>% group_by(strain) %>% summarise(
                          Sum = sum(cumulative_viability))
avg_strain_viability = cumulative_viability %>% group_by(strain) %>% summarise(
                          Mean = mean(avg_viability), std = std.error(cumulative_viability),
                          N = sum(embryo_number))
alpha = 0.05
avg_strain_viability$degrees_freedom = avg_strain_viability$N - 1
avg_strain_viability$t_score = qt(p=alpha/2, df=avg_strain_viability$degrees_freedom,lower.tail=F)
avg_strain_viability$margin_error <- avg_strain_viability$t_score * avg_strain_viability$std
avg_strain_viability$lower_bound <- avg_strain_viability$Mean - avg_strain_viability$margin_error
avg_strain_viability$upper_bound <- avg_strain_viability$Mean + avg_strain_viability$margin_error


## ---------------------------------------------------------------------------------
cumulative_brood_size_vis = ggplot(avg_strain_brood, aes(x = strain, y = Mean))
cumulative_brood_size_vis + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=Mean - std, ymax=Mean + std), width=.2) + theme_cowplot(14)

## ---------------------------------------------------------------------------------
ggsave(paste0(viability_and_brood_datafile_directory,substring(brood_data_filename,1,7),"_cumulative_brood_size_plot.pdf"), width = 20, height = 20, units = "cm", device = "pdf", dpi = 300)


