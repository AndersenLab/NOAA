library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(cowplot)
library(stationaRy)
library(ggplot2)

setwd("/Users/katieevans/Dropbox/AndersenLab/LabFolders/Katie/projects/noaa/data/")

#Read processed_mapping_dataframes
load("20161103.noaa.pmd.all.RData")
pmd3mo <- processed_mapping_df_3mo %>% mutate(trait = paste0(trait, "_3mo"))
pmd1yr <- processed_mapping_df_1yr %>% mutate(trait = paste0(trait, "_1yr"))
pmd3yr <- processed_mapping_df_3yr %>% mutate(trait = paste0(trait, "_3yr"))
pmdloc <- readRDS("/Users/katieevans/Dropbox/AndersenLab/LabFolders/Katie/projects/noaa/data/now/20160715.pmd.152.location.rds") %>%
  mutate(time_period = c(NA))
processed_mapping_df <- rbind(pmdloc, pmd3mo, pmd1yr, pmd3yr)

#Read wild isolate table
wi <- read.csv('/Users/katieevans/Dropbox/AndersenLab/LabFolders/Katie/projects/noaa/data/now/20160715.wild.isolates.adjusted.csv', header=T) %>% 
  dplyr::select(-X)

wi_stations <- wi_weather

#Get all NOAA weather stations
all_stations <- get_isd_stations() %>%
  mutate(id = paste(usaf, wban, sep = "-"))

#######################################FIGURE 1###########################################################

#Map of worms and all stations
mapWormStation <- ggplot() + 
  borders("world", ylim = c(-60, 90), fill="#F2EFE8") + 
  geom_point(aes(x=all_stations$lon, y=all_stations$lat), color="dodgerblue3", size=0.05) +
  geom_point(aes(x=wi$longitude, y = wi$latitude), color="red", size=.8, stroke = 1) +
  theme(panel.background = element_rect(fill = "#B7D1D1", colour = "#D9B8D8"),
        panel.grid.major = element_line(colour = "grey90"),
        axis.text=element_text(size=12), 
        axis.title=element_text(size=12, face = "bold"),
        axis.line.y = element_line(size = 0),
        axis.line.x = element_line(size = 0))+
  xlab("Longitude") + 
  ylab("Latitude") +
  ylim(-60, 90)

#Make histogram of station to wi for 3 year
wi_stations_3yr <- wi_stations %>% 
  filter(time_period == "3yr")

station_hist <- ggplot(data=wi_stations_3yr, aes(wi_stations_3yr$station_distance)) + 
  geom_histogram(binwidth = 15, fill = "dodgerblue3") + 
  xlab("Distance (km)") + 
  ylab("Frequency") +
  theme_bw() + 
  theme(axis.text=element_text(size=12, color = "black"), 
        axis.title=element_text(size=12, color = "black", face = "bold")) 

plot(plot_grid(mapWormStation, station_hist, rel_widths = c(2.75,1), labels = c("A", "B")))
cowplot::ggsave(paste0(getwd(), "/", Sys.Date(), ".mapWormStation.eps"), width = 7.5, height = 3.5)

pmd_full <- processed_mapping_df

####################################### FIGURE 2 & 4 ###########################################################
hm<-pmd_full
hm<-distinct(pmd_full, trait,strain,peak_id, .keep_all = T)
hm$allele <- factor(hm$allele,
                    levels = c(-1,1),
                    labels = c("REF", "ALT"))

gwasPxG <- function(trt,specific_peak,ylab){
  hm %>%
    filter(trait==trt,peak_id==specific_peak,!is.na(allele))%>%
    ggplot(.)+
    aes(x=allele,y = value,fill=as.factor(allele))+
    geom_boxplot(outlier.shape=NA,size =.5,color="gray52")+
    geom_point(size = .8,position=position_jitter(w=.4,  h=.025),na.rm=TRUE)+
    theme_bw()+
    theme(axis.text.x = element_text(size=9, color="black"),
          axis.text.y = element_text(size=9, color="black"),
          axis.title = element_text(size=9, color="black", vjust = 1, face = "bold"),
          strip.text.x = element_text(size=9, color="black"),
          strip.text.y = element_text(size=9, color="black"),
          plot.title = element_text(size=9, vjust=1),
          legend.title = element_text(size=9),
          panel.border = element_rect(size=1, colour = "black"),
          plot.margin = unit(c(1,.05,.05,.05), "cm"),
          legend.position = "none")+
    scale_y_continuous(breaks= pretty_breaks())+
    labs( x = "Genotype",y=ylab)+
    scale_fill_manual( values = c("darkgray", "burlywood2", "darkolivegreen","black"))
}

for(i in c("elevation", "temp_avg_3yr")) {
  pmd <- pmd_full %>% filter(trait == i)
  selection <- pmd %>% filter(trait == i)
  max <- max(selection$log10p)
    assign(paste0("plot_", i), pmd %>% 
             dplyr::filter(trait == i) %>% 
             dplyr::distinct(marker, .keep_all = T) %>% 
             ggplot2::ggplot(.) + 
             ggplot2::aes(x = POS/1e+06, y = log10p, fill = BF) + 
             ggplot2::geom_rect(mapping=aes(xmin=startPOS/1e6, xmax=endPOS/1e6, ymin=0, ymax= Inf),fill="thistle1", alpha=1) +
             ggplot2::geom_point(aes(color = ifelse(log10p > BF, 'red', 'black')), size = 0.8) +
             ggplot2::facet_grid(. ~ CHROM, scales = "free_x", space = "free_x") + 
             ggplot2::geom_hline(ggplot2::aes(yintercept = BF), color = "grey60", linetype = "dashed") + 
             ggplot2::theme(strip.background = element_blank(),
                            strip.text.x = element_text(size = 9, face = "bold"),
                            panel.background = element_rect(fill = "white"),
                            panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
                            panel.margin = unit(.6, "lines"),
                            panel.background = element_blank(),
                            axis.ticks = element_line(colour = "black"),
                            axis.text.y = element_text(colour = "black", size = 9),
                            axis.text.x = element_blank(),
                            axis.title = element_text(colour = "black", size = 9, face = "bold"),
                            axis.line.x = element_line(size = 0),
                            axis.line.y = element_line(size = 0),
                            legend.title = element_text(size = 10, colour = "blue", angle = 270),
                            legend.text = element_blank(),
                            legend.key.size=unit(0,"cm"),
                            legend.key = element_rect(colour = "pink"),
                            legend.position=('none'))+
             ggplot2::labs(x = "Genomic Position (Mb)", 
                           y = expression(bold(-log["10"](p))), colour = "red") +
             ggplot2::scale_color_identity()+
             scale_y_continuous(breaks= pretty_breaks(),expand=c(0,0),limits=c(0,max+.075*max),labels = function(x) format(x,width = 4))
    )
      
  if(i == "elevation") {
    label <- "Elevation (m)"
  } else if(i == "temp_avg_3yr") {
    label <- "Temperature (°C)"
  }
  setwd("~/Dropbox/AndersenLab/LabFolders/Katie/projects/noaa/figures")
  plot(plot_grid(get(paste0("plot_", i)), gwasPxG(i,1,label), rel_widths = c(2, 1), labels = c("A", "B")))
  ggsave(paste0(getwd(), "/", Sys.Date(), "_", i, ".manpxg.eps"), width = 7.5, height = 2.65, units = "in")
}

####################################### FIGURE 3 ###########################################################

#Filter to plot only the good traits for the main figure
processed_mapping_df <- processed_mapping_df %>% filter(trait == "rh_avg_3mo" | trait == "rh_avg_1yr" | trait == "rh_avg_3yr")

#Make manhattan plots (adapted from cegwas::manplot)
count <- 1
for(i in unique(processed_mapping_df$trait)) {
  selection <- processed_mapping_df %>% filter(trait == i)
  max <- max(selection$log10p)
  if(count == 1) {
    assign(paste0("plot_", i), processed_mapping_df %>% 
             dplyr::filter(trait == i) %>% 
             dplyr::distinct(marker, .keep_all = T) %>% 
             ggplot2::ggplot(.) + 
             ggplot2::aes(x = POS/1e+06, y = log10p, fill = BF) + 
             ggplot2::geom_rect(mapping=aes(xmin=startPOS/1e6, xmax=endPOS/1e6, ymin=0, ymax= Inf),fill="thistle1", alpha=1) +
             ggplot2::geom_point(aes(color = ifelse(log10p > BF, 'red', 'black')), size = 0.8) +
             ggplot2::facet_grid(. ~ CHROM, scales = "free_x", space = "free_x") + 
             ggplot2::geom_hline(ggplot2::aes(yintercept = BF), color = "grey60", linetype = "dashed") + 
             ggplot2::theme(strip.background = element_blank(),
                            strip.text.x = element_text(size = 12, face = "bold"),
                            panel.background = element_rect(fill = "white"),
                            panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
                            panel.margin = unit(.6, "lines"),
                            panel.background = element_blank(),
                            axis.ticks = element_line(colour = "black"),
                            axis.text.x = element_blank(),
                            axis.text.y = element_text(colour = "black", size = 9),
                            axis.title.y = element_blank(),
                            axis.line.x = element_line(size = 0),
                            axis.line.y = element_line(size = 0),
                            axis.title = element_text(size = 9),
                            plot.margin=unit(c(.05,.30,-.45,.50), "cm"),
                            legend.title = element_text(size = 10, colour = "blue", angle = 270),
                            legend.text = element_blank(),
                            legend.key.size=unit(0,"cm"),
                            legend.key = element_rect(colour = "pink"),
                            legend.position=('none'))+
            ggplot2::labs(x = "", y = "", colour = "red") +
            ggplot2::scale_color_identity()+
            scale_y_continuous(breaks= pretty_breaks(),expand=c(0,0),limits=c(0,max+.075*max),labels = function(x) format(x,width = 4))
    )
  } else {
    assign(paste0("plot_", i), processed_mapping_df %>% 
             dplyr::filter(trait == i) %>% 
             dplyr::distinct(marker, .keep_all = T) %>% 
             ggplot2::ggplot(.) + 
             ggplot2::aes(x = POS/1e+06, y = log10p, fill = BF) + 
             ggplot2::geom_rect(mapping=aes(xmin=startPOS/1e6, xmax=endPOS/1e6, ymin=0, ymax= Inf),fill="thistle1", alpha=1) +
             ggplot2::geom_point(aes(color = ifelse(log10p > BF, 'red', 'black')), size = 0.8) +
             ggplot2::facet_grid(. ~ CHROM, scales = "free_x", space = "free_x") + 
             ggplot2::geom_hline(ggplot2::aes(yintercept = BF), color = "grey60", linetype = "dashed") + 
             ggplot2::theme(strip.background = element_blank(),
                            strip.text.x = element_blank(),
                            panel.background = element_rect(fill = "white"),
                            panel.border = element_rect(color="black", size=0.5, linetype="solid", fill=NA),
                            panel.margin = unit(.6, "lines"),
                            panel.background = element_blank(),
                            axis.ticks = element_line(colour = "black"),
                            axis.text.x = element_blank(),
                            axis.text.y = element_text(colour = "black", size = 9),
                            axis.title.y = element_blank(),
                            axis.line.x = element_line(size = 0),
                            axis.line.y = element_line(size = 0),
                            axis.title = element_text(size = 9),
                            plot.margin=unit(c(.05,.30,-.45,.50), "cm"),
                            legend.title = element_text(size = 10, colour = "blue", angle = 270),
                            legend.text = element_blank(),
                            legend.key.size=unit(0,"cm"),
                            legend.key = element_rect(colour = "pink"),
                            legend.position=('none'))+
             ggplot2::labs(x = "", y = "", colour = "red") +
             ggplot2::scale_color_identity()+
             scale_y_continuous(breaks= pretty_breaks(),expand=c(0,0),limits=c(0,max+.075*max),labels = function(x) format(x,width = 4))
    )
  }
  count <- count + 1
}

#Plot manhattan plots with labels
all_traits <- plot_grid(plot_rh_avg_3mo, plot_rh_avg_1yr, plot_rh_avg_3yr, align = "v", ncol = 1, labels = c("A", "B", "C"))
df <- data.frame(1,2)
blank_plot<-ggplot(df,aes(x=1,y=1)) + geom_point(color="white") + theme(axis.line=element_blank(),axis.text =element_blank(),axis.ticks =element_blank(),axis.title =element_blank(),panel.background = element_blank(),panel.grid = element_blank())
all_traits<-plot_grid(all_traits,blank_plot, ncol=1,rel_heights = c(1, .03))
all_traits <- all_traits + draw_label(expression(bold(-log["10"](p))), x = .02, y = 0.5, hjust = .5, vjust = .5,
                                      fontfamily = "", fontface = "bold", colour = "black", size = 12,
                                      angle = 90, lineheight = 0.9, alpha = 1)
all_traits<- all_traits + draw_label("Genomic Position (Mb)", x = .5, y = 0.020, hjust = .5, vjust = .5,
                           fontfamily = "", fontface = "bold", colour = "black", size = 12,
                           angle = 0, lineheight = 0.9, alpha = 1)
plot(all_traits)
ggsave(paste0(getwd(), "/", Sys.Date(), "_", i, ".manpxg.eps"), width = 7.5, height = 6.5, units = "in")

####################################### FIGURE 5 ###########################################################

#Competition Assay
#Read in data and remove bad replicate (based on McGrath lab)
data25 <- read.csv("~/Dropbox/AndersenLab/LabFolders/Katie/projects/noaa/data/competition/ddpcr_results_25.csv")
data25 <- data25 %>% mutate(temp = 25)
data25[11,3] <- NA 
data25[11,5] <- NA
data25$Name <- as.character(data25$Name)
g0_25 <- c("g0", 100, 100, 50, 50, 25)
data25 <- rbind(g0_25, data25)
data15 <- read.csv("~/Dropbox/AndersenLab/LabFolders/Katie/projects/noaa/data/competition/ddpcr_results_15.csv")
data15 <- data15 %>% 
  mutate(temp = 15)
colnames(data15)[1] <- "Name"
data15$Name <- as.character(data15$Name)
g0_15 <- c("g0", 100, 100, 50, 50, 15)
data15 <- rbind(g0_15, data15)

#Find averages and proportion of JU to CX strains
data <- rbind(data15, data25)
data$Ch1.1 <- as.numeric(data$Ch1.1)
data$Ch1.2 <- as.numeric(data$Ch1.2)
data$Ch2.1 <- as.numeric(data$Ch2.1)
data$Ch2.2 <- as.numeric(data$Ch2.2)
data <- data %>%
  mutate(Ch1_avg = (Ch1.1 + Ch1.2)/2, 
         Ch2_avg = (Ch2.1 + Ch2.2)/2,
         JU_CX_prop = Ch2_avg / Ch1_avg,
         gen = extract_numeric(Name)) %>%
  dplyr::select(Name, gen, temp, Ch1.1, Ch1.2, Ch1_avg, Ch2.1, Ch2.2, Ch2_avg, JU_CX_prop) %>%
  mutate(gen_temp = paste0(gen, "_", temp))

#Plot line graph
avg <- data %>%
  group_by(gen_temp) %>%
  summarise_each(funs(mean(., na.rm=TRUE), sd(., na.rm=T))) %>%
  dplyr::select(gen_temp, JU_CX_prop_mean, JU_CX_prop_sd) %>%
  separate(gen_temp, c("gen", "temp"))

data$gen <- as.factor(data$gen)

plotA <- ggplot(data=avg, 
                aes(x=gen, y=JU_CX_prop_mean, color=temp, group=temp)) + 
  scale_colour_manual(values = c("#377DB8", "#E51A1D")) +
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymax=JU_CX_prop_mean+JU_CX_prop_sd, 
                    ymin=JU_CX_prop_mean-JU_CX_prop_sd), 
                width=0.15) +
  theme_bw() +
  theme(axis.text=element_text(size=9), 
        axis.title=element_text(size=9, face = "bold"),
        legend.position='none') +
  ylim(0.0, 0.6) +
  labs(x="Generation", y="JU847 Allele Frequency")

plot(plotA)
ggsave(paste0(getwd(), "/", Sys.Date(), ".competition.eps"), width = 6, height = 4.5, units = "in")

