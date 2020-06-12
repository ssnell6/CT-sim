
# Set options and load libraries
options(stringsAsFactors=F)
library(CTSim)
library(tidyverse)
library(cowplot)
library(abind)
library(reshape2)
library(raster)

# Define working directory
sum_dir = 'Summaries/EXP4'
fig_dir = 'Results/Plots/EXP4'
sim_dir = 'Code'
data_dir = 'Z:/Snell/core-transient-simulation-LL19/Results/EXP4/d-g2_hp-0.5'

xclass.out = c()
land.out = c()
xclass.ct = c()
scale = 3
data_dirs =c("Results/EXP4/d-g2_hp-0.5","Results/EXP4/d-g2_hp-0.6", "Results/EXP4/d-g2_hp-0.7", "Results/EXP4/d-g2_hp-0.8", "Results/EXP4/d-g2_hp-0.9")
for(d in data_dirs){
  for (r in 1:50) {
    data_dir = d
    d_hp = substring(d, 14)
    temp = pixelXclass(data_dir, d_hp, run=r, scale = 3, t_window = 186:200,
                       ct_threshold = 1/3, return = 'percent')
    tmp = pixelXclass(data_dir, d_hp, run=r, scale = 3, t_window = 186:200,
                      ct_threshold = 1/3, return = 'count')
    
    # Eliminate pixels within a distance 'scale' from the grid edge
    temp2 = temp$xclass[temp$landSim$x > scale & temp$landSim$x < max(temp$landSim) - scale &
                          temp$landSim$y > scale & temp$landSim$y < max(temp$landSim) - scale ,,]
    xclass.out = abind(xclass.out, temp2, along = 1)
    
    tmp2 = tmp$xclass[tmp$landSim$x > scale & tmp$landSim$x < max(tmp$landSim) - scale &
                        tmp$landSim$y > scale & tmp$landSim$y < max(tmp$landSim) - scale ,,]
    xclass.ct = abind(xclass.ct, tmp2, along = 1)
    
    
    land2 = temp$landSim[temp$landSim$x > scale & temp$landSim$x < max(temp$landSim) - scale &
                           temp$landSim$y > scale & temp$landSim$y < max(temp$landSim) - scale ,,]
    land2$d = ""
    land2$hp = ""
    land.out = rbind(land.out, land2) %>%
      mutate(hp = substring(d_hp, 9),
             d = substr(d_hp, 3, 4))
    
    print(paste(r, Sys.time()))
  }
}

save(xclass.out, xclass.ct, land.out, file = 'Results/Summary/EXP4/g2/pixel_xclass_summary.Rdata')

# Conduct cross-classification by species and pixel
#   (generates ~7M rows, may have memory issues)
# selected 496 bc it is midpoints between 480 and 512, the middle two rows of landscape
pix_all = c(448:544)
xclass.sp = c()
data_dirs =c("Results/EXP4/d-g2_hp-0.5", "Results/EXP4/d-g2_hp-0.6", "Results/EXP4/d-g2_hp-0.7", "Results/EXP4/d-g2_hp-0.8", "Results/EXP4/d-g2_hp-0.9") # 
for(d in data_dirs){
  for (r in 21:50) { # changed this to fewer runs
    for(pixel in c(500)){
      data_dir = d
      d_hp = substring(d, 14)
      tmp = pixelXclassBySpecies(data_dir, d_hp, run=r, scale = 3, t_window = 186:200, ct_threshold = 1/3) %>%
        filter(pix  == pixel) %>%
        mutate(hp = substring(d_hp, 9),
               d = substr(d_hp, 3, 4),
               run = r,
               landscape_sim = sim) %>% # added by SJST to filter to a single or set of focal pixed
        group_by(p) %>% # group by detection to generate length
        mutate(species_total = length(p)) %>% 
        group_by(pix, hp, p, species_total, run, landscape_sim, xc) %>%
        dplyr::summarise(count_xclass = n()) %>% # calculate each type of xclass within the group
        mutate(count_pct = count_xclass/species_total)
      
      xclass.sp = rbind(xclass.sp, tmp)
      rm(tmp)
      print(paste(r, d, Sys.time()))
    }
  }
}

save(xclass.sp, file = 'Results/Summary/EXP4/g2/pixel_xclass_summary_bysp_w_runs_1_50.Rdata')




##### pixel summary ######
load('Results/Summary/EXP4/g2/pixel_xclass_summary_bysp_w_runs_1_50.Rdata')
load('Results/Summary/EXP4/g2/pixel_xclass_summary.Rdata')

##### plotting ####
library(tidyverse)
library(cowplot)

landscape_bins <- xclass.sp %>%
  mutate(landscape_bin = case_when(landscape_sim > 0 & landscape_sim <= 0.1 ~ 0.1,
                                   landscape_sim > 0.2 & landscape_sim <= 0.2 ~ 0.2,
                                   landscape_sim > 0.2 & landscape_sim <= 0.3 ~ 0.3,
                                   landscape_sim > 0.3 & landscape_sim <= 0.4 ~ 0.4,
                                   landscape_sim > 0.4 & landscape_sim <= 0.5 ~ 0.5,
                                   landscape_sim > 0.5 & landscape_sim <= 0.6 ~ 0.6,
                                   landscape_sim > 0.6 & landscape_sim <= 0.7 ~ 0.7,
                                   landscape_sim > 0.7 & landscape_sim <= 0.8 ~ 0.8,
                                   landscape_sim > 0.8 & landscape_sim <= 0.9 ~ 0.9,
                                   landscape_sim > 0.9 & landscape_sim <= 1.0 ~ 1.0)) %>%
  dplyr::select(run, p, hp, landscape_bin) %>%
  distinct()

palette = colorRampPalette(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac'))

# raw plot #
raw_core <- xclass.sp %>%
  filter(xc %in% c("c-c")) %>%
  group_by(p, run, hp) %>%
  left_join(landscape_bins[,c("p","run", "hp", "landscape_bin")], by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  summarise(total_bio_core = mean(count_pct)) %>%
  ggplot(aes(p, landscape_bin)) +  geom_tile(aes(fill = total_bio_core)) +
  theme_classic() + 
  scale_fill_gradient2(midpoint = 0.5, 
                       name = "Proportion\nof core species",low = "#2166ac", mid = "#fddbc7", high ="#b2182b" ) +
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 8)) +
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + xlab("Detection probability") + ylab("Landscape similarity") +theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25)) 
ggsave("Results/Plots/EXP4/FigS1rawcore.pdf", height = 14, width = 16)

#### Figure 4 ####
# pull out df from tcore to generate mod. Need to fit BEFORE take the mean
tcoremod = xclass.sp %>%
  filter(xc %in% c("c-c", "c-t", "c-int")) %>%
  group_by(p, run, hp) %>%
  summarise(total_bio_core = sum(count_xclass)) %>%
  left_join(xclass.sp[,c("run", "hp", "p", "landscape_sim")], by = c("run", "hp", "p")) %>% # get landscape sim back and use this
  distinct()
tcore_mod <- lm(total_bio_core ~ p+landscape_sim, data = tcoremod)

# tcore var par
tcore_p <- lm(tcoremod$total_bio_core ~  tcoremod$p)  
# z scores separated out for env effects 
tcore_ls = lm(tcoremod$total_bio_core ~  tcoremod$landscape_sim)
# z scores separated out for env effects
tcore_both = lm(tcoremod$total_bio_core ~ tcoremod$p + tcoremod$landscape_sim)

p = summary(tcore_both)$r.squared - summary(tcore_ls)$r.squared #p only
ls = summary(tcore_both)$r.squared - summary(tcore_p)$r.squared #ls only
shared = summary(tcore_p)$r.squared - ls #shared variance
none = 1 - summary(tcore_both)$r.squared # neither

ttrans = xclass.sp %>%
  filter(xc %in% c("t-t", "t-c", "t-int")) %>%
  group_by(p, run, hp) %>%
  summarise(total_bio_trans = sum(count_xclass)) %>%
  left_join(xclass.sp[,c("run", "hp", "p", "landscape_sim")], by = c("run", "hp", "p")) %>% 
  distinct()

ttransmod <- lm(total_bio_trans ~ p+landscape_sim, data = ttrans)

# ttrans var par
ttrans_p <- lm(ttrans$total_bio_trans ~  ttrans$p)  
# z scores separated out for env effects 
ttrans_ls = lm(ttrans$total_bio_trans ~  ttrans$landscape_sim)
# z scores separated out for env effects
ttrans_both = lm(ttrans$total_bio_trans ~ ttrans$p + ttrans$landscape_sim)

p = summary(ttrans_both)$r.squared - summary(ttrans_ls)$r.squared #p only
ls = summary(ttrans_both)$r.squared - summary(ttrans_p)$r.squared #ls only
shared = summary(ttrans_p)$r.squared - ls #shared variance
none = 1 - summary(ttrans_both)$r.squared # neither


tcore = xclass.sp %>%
  filter(xc %in% c("c-c", "c-t", "c-int")) %>%
  group_by(p, run, hp) %>%
  summarise(total_bio_core = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  summarise(mean_bio_core = mean(total_bio_core),
            sd = sd(total_bio_core), 
            cv = sd/mean_bio_core) %>%
  ggplot(aes(p, landscape_bin)) +  geom_tile(aes(fill = mean_bio_core)) +
  theme_classic() + 
  scale_fill_gradient2(midpoint = 12, 
                       name = "Core species",low = "#2166ac", mid = "#fddbc7", high ="#b2182b", breaks = c(8,12,16)) +
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 8)) +
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + xlab("Detection probability") + ylab("Landscape similarity") +theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25)) 

tcore_sd = xclass.sp %>%
  filter(xc %in% c("c-c", "c-t", "c-int")) %>%
  group_by(p, run, hp) %>%
  summarise(total_bio_core = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  summarise(mean_bio_core = mean(total_bio_core),
            sd = sd(total_bio_core),
            cv = sd/mean_bio_core) %>%
  ggplot(aes(p, landscape_bin)) +  geom_tile(aes(fill = cv)) +
  theme_classic() + 
  scale_fill_gradient2(midpoint = 0.3, 
                       name = "Core species",low = "#2166ac", mid = "#fddbc7", high ="#b2182b") +
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 8)) +
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + xlab("Detection probability") + ylab("Landscape similarity") +theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25)) 

ttrans = xclass.sp %>%
  filter(xc %in% c("t-c", "t-t", "t-int")) %>%
  group_by(p, run, hp) %>%
  summarise(total_bio_trans = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  summarise(mean_bio_trans = mean(total_bio_trans),
            sd = sd(total_bio_trans)) %>%
  ggplot(aes(p, landscape_bin)) +  geom_tile(aes(fill = mean_bio_trans)) +
  theme_classic() + 
  scale_fill_gradient2(midpoint = 10, 
                       name = "Transient \nspecies",low = "#2166ac", mid = "#fddbc7", high ="#b2182b") +
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 8)) +
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + xlab("Detection probability") + ylab("Landscape similarity") +theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25))

ttrans_sd = xclass.sp %>%
  filter(xc %in% c("t-c", "t-t", "t-int")) %>%
  group_by(p, run, hp) %>%
  summarise(total_bio_trans = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  summarise(mean_bio_trans = mean(total_bio_trans),
            sd = sd(total_bio_trans),
            cv = sd/mean_bio_trans) %>%
  ggplot(aes(p, landscape_bin)) +  geom_tile(aes(fill = cv)) +
  theme_classic() + 
  scale_fill_gradient2(midpoint = 0.4, 
                       name = "Transient \nspecies",low = "#2166ac", mid = "#fddbc7", high ="#b2182b") +
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 8)) +
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + xlab("Detection probability") + ylab("Landscape similarity") +theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25))

theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
fig4 <- plot_grid(tcore + theme(legend.position="right"),
                  ttrans + theme(legend.position="right"),
                  align = 'hv',
                  labels = c("A","B"),
                  label_size = 30,
                  nrow = 2) 
ggsave("Results/Plots/EXP4/g2/Fig4_t_core_trans.pdf", height = 17, width = 14)

#### line plots fig 4 ######
tcore_det_single = xclass.sp %>%
  filter(xc %in% c("c-c", "c-t", "c-int")) %>%
  group_by(p, run, hp) %>%
  summarise(total_bio_core = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  summarise(mean_bio_core = mean(total_bio_core),
            sd = sd(total_bio_core),
            lower = mean_bio_core - 1.96*sd,
            upper = mean_bio_core + 1.96*sd) %>%
  filter(landscape_bin %in% c(0.3, 0.8)) %>%
  mutate(`Landscape \nsimilarity` = as.factor(landscape_bin)) %>%
  ggplot(aes(p, mean_bio_core)) + 
  geom_point(size =5, col = "goldenrod2", aes(shape = `Landscape \nsimilarity`)) + 
  geom_line(lwd = 1.5, col = "goldenrod2", aes(lty = `Landscape \nsimilarity`)) +
  scale_y_continuous(limits = c(0, 20)) +
  # geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) +
  scale_shape_manual(values = c(15, 16))+
  xlab("Detection probability") + ylab("Mean count of core species") +
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), axis.text.x=element_text(size=25),axis.text.y=element_text(size=25), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25)) 

ttrans_det_single = xclass.sp %>%
  filter(xc %in% c("t-c", "t-t", "t-int")) %>%
  group_by(p, run, hp) %>%
  summarise(total_bio_trans = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  summarise(mean_bio_trans = mean(total_bio_trans),
            sd = sd(total_bio_trans),
            lower = mean_bio_trans - 1.96*sd,
            upper = mean_bio_trans + 1.96*sd) %>% 
  filter(landscape_bin %in% c(0.3, 0.8)) %>%
  mutate(`Landscape \nsimilarity` = as.factor(landscape_bin)) %>%
  ggplot(aes(p, mean_bio_trans)) + 
  geom_point(size =5, col = "goldenrod2", aes(shape = `Landscape \nsimilarity`)) + 
  geom_line(lwd = 1.5,col = "goldenrod2", aes(lty = `Landscape \nsimilarity`)) +
  # geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) +
  geom_line(lwd = 1.5, col = "goldenrod2", aes(lty = `Landscape \nsimilarity`)) +
  scale_y_continuous(limits = c(0, 20)) +
  scale_shape_manual(values = c(15, 16))+
  xlab("Detection probability") + ylab("Mean count of transient species") +
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), axis.text.x=element_text(size=25),axis.text.y=element_text(size=25), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25)) 

tcore_ls_single = xclass.sp %>%
  filter(xc %in% c("c-c", "c-t", "c-int")) %>%
  group_by(p, run, hp) %>%
  summarise(total_bio_core = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  summarise(mean_bio_core = mean(total_bio_core)) %>%
  filter(p %in% c(0.1, 0.9)) %>%
  mutate(`Detection` = as.factor(p)) %>%
  ggplot(aes(landscape_bin, mean_bio_core)) + 
  geom_point(size =5, aes(shape = `Detection`)) + 
  scale_y_continuous(limits = c(0, 20)) +
  geom_line(lwd = 1.5, aes(lty = `Detection`)) +
  scale_shape_manual(values = c(15, 16))+
  xlab("Landscape similarity") + ylab("Mean count of core species") +
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), axis.text.x=element_text(size=25),axis.text.y=element_text(size=25), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25)) 

ttrans_ls_single = xclass.sp %>%
  filter(xc %in% c("t-c", "t-t", "t-int")) %>%
  group_by(p, run, hp) %>%
  summarise(total_bio_trans = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  summarise(mean_bio_trans = mean(total_bio_trans)) %>% 
  filter(p %in% c(0.1, 0.9)) %>%
  mutate(`Detection` = as.factor(p)) %>%
  ggplot(aes(landscape_bin, mean_bio_trans)) + 
  geom_point(size =5, aes(shape = `Detection`)) + 
  geom_line(lwd = 1.5, aes(lty = `Detection`)) +
  scale_y_continuous(limits = c(0, 20)) +
  scale_shape_manual(values = c(15, 16)) +
  xlab("Landscape similarity") + ylab("Mean count of transient species") +
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), axis.text.x=element_text(size=25),axis.text.y=element_text(size=25), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25)) 

theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
grid <- plot_grid(tcore_det_single + theme(legend.position="right"),
                  tcore_ls_single + theme(legend.position="right"),
                  ttrans_det_single + theme(legend.position="right"),
                  ttrans_ls_single + theme(legend.position="right"),
                  align = 'hv',
                  # labels = c("A","B", "C", "D"),
                  # label_size = 30,
                  # hjust = -2,
                  nrow = 2) 
ggsave("Results/Plots/EXP4/g2/Fig4_lines.pdf", height = 14, width = 18)  

core <- plot_grid(tcore + theme(legend.position="right"),
                  tcore_det_single + theme(legend.position="right"),
                  tcore_ls_single + theme(legend.position="right"),
                  align = 'hv',
                  ncol = 3)

trans <- plot_grid(ttrans + theme(legend.position="right"),
                   ttrans_det_single + theme(legend.position="right"),
                   ttrans_ls_single + theme(legend.position="right"),
                   align = 'hv',
                   ncol = 3)

plot_grid(core, 
          trans,
          align = 'hv',
          nrow = 2)
ggsave("Results/Plots/EXP4/g2/Fig4_grid_lines.pdf", height = 14, width = 30)  

#### Figure 5 ####
# t.test(ncore$`c-t`, ncore$`c-c`)
# change color ramp of range based on midpoint, change line plots and fig 5b and all stats
ncore = xclass.sp %>%
  filter(xc %in% c("c-c", "c-t", "c-int")) %>%
  group_by(p, run, hp, xc) %>%
  summarise(total_bio_core = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  spread(xc, total_bio_core) %>%
  mutate(`c-c_0` = ifelse(is.na(`c-c`), 0, `c-c`),
         `c-t_0` = ifelse(is.na(`c-t`), 0, `c-t`),
         `c-int_0` = ifelse(is.na(`c-int`), 0, `c-int`)) %>%
  mutate(perc_ncore = `c-t_0`/(`c-c_0`+`c-t_0`+`c-int_0`)) %>% 
  group_by(p, landscape_bin) %>%
  summarise(mean_perc_core = mean(perc_ncore),
            sd = sd(perc_ncore)) %>%
  ggplot(aes(p, landscape_bin)) + geom_tile(aes(fill = mean_perc_core*100)) +
  theme_classic() + 
  scale_fill_gradient2(midpoint = 33, 
                       name = "Core species: \n% incorrect",low = "#2166ac", mid = "#fddbc7", high ="#b2182b") +
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 8)) +
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + xlab("Detection probability") + ylab("Landscape similarity") +theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25))

ncore_sd = xclass.sp %>%
  filter(xc %in% c("c-c", "c-t", "c-int")) %>%
  group_by(p, run, hp, xc) %>%
  summarise(total_bio_core = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  spread(xc, total_bio_core) %>%
  mutate(`c-c_0` = ifelse(is.na(`c-c`), 0, `c-c`),
         `c-t_0` = ifelse(is.na(`c-t`), 0, `c-t`),
         `c-int_0` = ifelse(is.na(`c-int`), 0, `c-int`)) %>%
  mutate(perc_ncore = `c-t_0`/(`c-c_0`+`c-t_0`+`c-int_0`)) %>% 
  group_by(p, landscape_bin) %>%
  summarise(mean_perc_core = mean(perc_ncore),
            sd = sd(perc_ncore),
            cv = sd/mean_perc_core) %>%
  ggplot(aes(p, landscape_bin)) +  geom_tile(aes(fill = cv)) +
  theme_classic() + 
  scale_fill_gradient2(midpoint = 0.9, 
                       name = "Core species: \n% incorrect",low = "#2166ac", mid = "#fddbc7", high ="#b2182b") +
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 8)) +
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + xlab("Detection probability") + ylab("Landscape similarity") +theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25))
# t.test(ntrans$`t-c`, ntrans$`t-t`)
ntrans = xclass.sp %>%
  filter(xc %in% c("t-c", "t-t", "t-int")) %>%
  group_by(p, run, hp, xc) %>%
  summarise(total_bio_trans = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  spread(xc, total_bio_trans) %>%
  mutate(`t-c_0` = ifelse(is.na(`t-c`), 0, `t-c`),
         `t-t_0` = ifelse(is.na(`t-t`), 0, `t-t`),
         `t-int_0` = ifelse(is.na(`t-int`), 0, `t-int`)) %>%
  mutate(perc_ntrans = `t-c_0`/(`t-c_0`+`t-t_0`+`t-int_0`)) %>% 
  group_by(p, landscape_bin) %>%
  summarise(mean_perc_trans = mean(perc_ntrans),
            sd = sd(perc_ntrans)) %>%
  ggplot(aes(p, landscape_bin)) +  geom_tile(aes(fill = mean_perc_trans*100)) +
  theme_classic() + 
  scale_fill_gradient2(midpoint = 22, 
                       name = "Transient species: \n% incorrect",low = "#2166ac", mid = "#fddbc7", high ="#b2182b", limits = c(0, 50)) +
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 8)) +
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + xlab("Detection probability") + ylab("Landscape similarity") +theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25))

ntrans_sd =  xclass.sp %>%
  filter(xc %in% c("t-c", "t-t", "t-int")) %>%
  group_by(p, run, hp, xc) %>%
  summarise(total_bio_trans = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  spread(xc, total_bio_trans) %>%
  mutate(`t-c_0` = ifelse(is.na(`t-c`), 0, `t-c`),
         `t-t_0` = ifelse(is.na(`t-t`), 0, `t-t`),
         `t-int_0` = ifelse(is.na(`t-int`), 0, `t-int`)) %>%
  mutate(perc_ntrans = `t-c_0`/(`t-c_0`+`t-t_0`+`t-int_0`)) %>% 
  group_by(p, landscape_bin) %>%
  summarise(mean_perc_trans = mean(perc_ntrans),
            sd = sd(perc_ntrans),
            cv = sd/mean_perc_trans) %>%
  ggplot(aes(p, landscape_bin)) +  geom_tile(aes(fill = cv)) +
  theme_classic() + 
  scale_fill_gradient2(midpoint = 3, 
                       name = "Transient species: \n% incorrect",low = "#2166ac", mid = "#fddbc7", high ="#b2182b") +
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 8)) +
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + xlab("Detection probability") + ylab("Landscape similarity") +theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25))

theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
grid <- plot_grid(ncore + theme(legend.position="right"),
                  ntrans + theme(legend.position="right"),
                  align = 'hv',
                  labels = c("A","B"),
                  label_size = 30,
                  nrow = 2) 
ggsave("Results/Plots/EXP4/g2/Fig5_n_core_trans.pdf", height = 18, width = 14)  


#### line plots fig 5 ######
tcore_det_single = xclass.sp %>%
  filter(xc %in% c("c-c", "c-t", "c-int")) %>%
  group_by(p, run, hp, xc) %>%
  summarise(total_bio_core = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  spread(xc, total_bio_core) %>%
  mutate(`c-c_0` = ifelse(is.na(`c-c`), 0, `c-c`),
         `c-t_0` = ifelse(is.na(`c-t`), 0, `c-t`),
         `c-int_0` = ifelse(is.na(`c-int`), 0, `c-int`)) %>%
  mutate(perc_ncore = `c-t_0`/(`c-c_0`+`c-t_0`+`c-int_0`)) %>% 
  filter(landscape_bin %in% c(0.3, 0.8)) %>%
  # summarise to get mean of points
  group_by(p, landscape_bin) %>%
  summarise(mean_perc_ncore = mean(perc_ncore),
            sd = sd(perc_ncore)) %>%
  mutate(`Landscape \nsimilarity` = as.factor(landscape_bin)) %>%
  ggplot(aes(p, mean_perc_ncore)) + 
  geom_point(size =5, col = "goldenrod2", aes(shape = `Landscape \nsimilarity`)) + 
  geom_line(lwd = 1.5,col = "goldenrod2", aes(lty = `Landscape \nsimilarity`)) +
  scale_y_continuous(limits = c(0, 0.6)) +
  scale_shape_manual(values = c(15, 16)) +
  xlab("Detection probability") + ylab("Core species: % incorrect") +
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), axis.text.x=element_text(size=25),axis.text.y=element_text(size=25),legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 20)) 

ttrans_det_single = xclass.sp %>%
  filter(xc %in% c("t-c", "t-t", "t-int")) %>%
  group_by(p, run, hp, xc) %>%
  summarise(total_bio_trans = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  spread(xc, total_bio_trans) %>%
  mutate(`t-c_0` = ifelse(is.na(`t-c`), 0, `t-c`),
         `t-t_0` = ifelse(is.na(`t-t`), 0, `t-t`),
         `t-int_0` = ifelse(is.na(`t-int`), 0, `t-int`)) %>% 
  mutate(perc_ntrans = `t-c_0`/(`t-c_0`+`t-t_0`+`t-int_0`)) %>% 
  filter(landscape_bin %in% c(0.3, 0.8)) %>%
  group_by(p, landscape_bin) %>%
  summarise(mean_perc_ntrans = mean(perc_ntrans),
            sd = sd(perc_ntrans)) %>%
  mutate(`Landscape \nsimilarity` = as.factor(landscape_bin)) %>%
  ggplot(aes(p, mean_perc_ntrans)) + 
  geom_point(size =5, col = "goldenrod2", aes(shape = `Landscape \nsimilarity`)) + 
  geom_line(lwd = 1.5,col = "goldenrod2", aes(lty = `Landscape \nsimilarity`)) +
  scale_y_continuous(limits = c(0, 0.6)) +
  scale_shape_manual(values = c(15, 16)) +
  xlab("Detection probability") + ylab("Transient species: % incorrect") +
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), axis.text.x=element_text(size=25),axis.text.y=element_text(size=25),legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 20)) 

tcore_ls_single = xclass.sp %>%
  filter(xc %in% c("c-c", "c-t", "c-int")) %>%
  group_by(p, run, hp, xc) %>%
  summarise(total_bio_core = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  spread(xc, total_bio_core) %>%
  mutate(`c-c_0` = ifelse(is.na(`c-c`), 0, `c-c`),
         `c-t_0` = ifelse(is.na(`c-t`), 0, `c-t`),
         `c-int_0` = ifelse(is.na(`c-int`), 0, `c-int`)) %>%
  mutate(perc_ncore = `c-t_0`/(`c-c_0`+`c-t_0`+`c-int_0`)) %>% 
  filter(p %in% c(0.1, 0.9)) %>%
  # summarise to get mean of points
  group_by(p, landscape_bin) %>%
  summarise(mean_perc_ncore = mean(perc_ncore),
            sd = sd(perc_ncore)) %>%
  mutate(`Detection` = as.factor(p)) %>%
  ggplot(aes(landscape_bin, mean_perc_ncore)) + 
  geom_point(size =5, aes(shape = `Detection`)) + 
  geom_line(lwd = 1.5,aes(lty = `Detection`)) +
  scale_y_continuous(limits = c(0, 0.6)) +
  scale_shape_manual(values = c(15, 16)) +
  xlab("Landscape similarity") + ylab("Core species: % incorrect") +
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), axis.text.x=element_text(size=25),axis.text.y=element_text(size=25),legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 20)) 

ttrans_ls_single = xclass.sp %>%
  filter(xc %in% c("t-c", "t-t", "t-int")) %>%
  group_by(p, run, hp, xc) %>%
  summarise(total_bio_trans = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  spread(xc, total_bio_trans) %>%
  mutate(`t-c_0` = ifelse(is.na(`t-c`), 0, `t-c`),
         `t-t_0` = ifelse(is.na(`t-t`), 0, `t-t`),
         `t-int_0` = ifelse(is.na(`t-int`), 0, `t-int`)) %>% 
  mutate(perc_ntrans = `t-c_0`/(`t-c_0`+`t-t_0`+`t-int_0`)) %>% 
  filter(p %in% c(0.1, 0.9)) %>%
  group_by(p, landscape_bin) %>%
  summarise(mean_perc_ntrans = mean(perc_ntrans),
            sd = sd(perc_ntrans)) %>%
  mutate(`Detection` = as.factor(p)) %>%
  ggplot(aes(landscape_bin, mean_perc_ntrans)) + 
  geom_point(size =5, aes(shape = `Detection`)) + 
  geom_line(lwd = 1.5,aes(lty = `Detection`)) +
  scale_y_continuous(limits = c(0, 0.6)) +
  scale_shape_manual(values = c(15, 16)) +
  xlab("Landscape similarity") + ylab("Transient species: % incorrect") +
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), axis.text.x=element_text(size=25),axis.text.y=element_text(size=25),legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 20)) 


theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
grid <- plot_grid(tcore_det_single + theme(legend.position="right"),
                  tcore_ls_single + theme(legend.position="right"),
                  ttrans_det_single + theme(legend.position="right"),
                  ttrans_ls_single + theme(legend.position="right"),
                  align = 'hv',
                  # labels = c("A","B", "C", "D"),
                  # label_size = 30,
                  # hjust = -2,
                  nrow = 2) 
ggsave("Results/Plots/EXP4/g2/Fig5_lines.pdf", height = 14, width = 18)  

core <- plot_grid(ncore + theme(legend.position="right"),
                  tcore_det_single + theme(legend.position="right"),
                  tcore_ls_single + theme(legend.position="right"),
                  align = 'hv',
                  ncol = 3)

trans <- plot_grid(ntrans + theme(legend.position="right"),
                   ttrans_det_single + theme(legend.position="right"),
                   ttrans_ls_single + theme(legend.position="right"),
                   align = 'hv',
                   ncol = 3)

plot_grid(core, 
          trans,
          align = 'hv',
          nrow = 2)
ggsave("Results/Plots/EXP4/g2/Fig5_grid_lines.pdf", height = 14, width = 30)  

#### tcorex ######
# pull out df from tcore to generate mod. Need to fit BEFORE take the mean
tcorex = xclass.sp %>%
  filter(xc %in% c("c-c", "c-t", "c-int")) %>%
  group_by(p, run, hp, xc) %>%
  summarise(total_bio_core = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  spread(xc, total_bio_core) %>%
  mutate(`c-c_0` = ifelse(is.na(`c-c`), 0, `c-c`),
         `c-t_0` = ifelse(is.na(`c-t`), 0, `c-t`),
         `c-int_0` = ifelse(is.na(`c-int`), 0, `c-int`)) %>%
  mutate(perc_ncore = `c-t_0`/(`c-c_0`+`c-t_0`+`c-int_0`)) %>%
  distinct()
tcorex_mod <- lm(perc_ncore ~ p+landscape_sim, data = tcorex)

# tcore var par
tcore_p <- lm(tcorex$perc_ncore ~  tcorex$p)  
# z scores separated out for env effects 
tcore_ls = lm(tcorex$perc_ncore ~  tcorex$landscape_sim)
# z scores separated out for env effects
tcore_both = lm(tcorex$perc_ncore ~ tcorex$p + tcorex$landscape_sim)

p = summary(tcore_both)$r.squared - summary(tcore_ls)$r.squared #p only
ls = summary(tcore_both)$r.squared - summary(tcore_p)$r.squared #ls only
shared = summary(tcore_p)$r.squared - ls #shared variance
none = 1 - summary(tcore_both)$r.squared # neither


ttransx = xclass.sp %>%
  filter(xc %in% c("t-c", "t-t", "t-int")) %>%
  group_by(p, run, hp, xc) %>%
  summarise(total_bio_trans = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  spread(xc, total_bio_trans) %>%
  mutate(`t-c_0` = ifelse(is.na(`t-c`), 0, `t-c`),
         `t-t_0` = ifelse(is.na(`t-t`), 0, `t-t`),
         `t-int_0` = ifelse(is.na(`t-int`), 0, `t-int`)) %>%
  mutate(perc_ntrans = `t-c_0`/(`t-c_0`+`t-t_0`+`t-int_0`)) %>%
  distinct()
ttransx_mod <- lm(perc_ntrans~ p+landscape_bin, data = ttransx)

# tcore var par
ttrans_p <- lm(ttransx$perc_ntrans ~  ttransx$p)  
# z scores separated out for env effects 
ttrans_ls = lm(ttransx$perc_ntrans ~  ttransx$landscape_sim)
# z scores separated out for env effects
ttrans_both = lm(ttransx$perc_ntrans ~ ttransx$p + ttransx$landscape_sim)

p = summary(ttrans_both)$r.squared - summary(ttrans_ls)$r.squared #p only
ls = summary(ttrans_both)$r.squared - summary(ttrans_p)$r.squared #ls only
shared = summary(ttrans_p)$r.squared - ls #shared variance
none = 1 - summary(ttrans_both)$r.squared # neither

##### Figure 6 ####
pixelSummary("Results/EXP4/d-g2_hp-0.5", "d-g2_hp-0.5", run=1:20, "Results/Plots/EXP4", plot.pdf = TRUE)

xclass.ngrid = c()
for(r in 1:50){
  tmp = pixelXclassBySpecies("Results/EXP4/d-g2_hp-0.5", "d-g2_hp-0.5", run=r, scale = 3, t_window = 186:200, ct_threshold = 1/3)
  xclass.ngrid = rbind(xclass.ngrid, tmp)
}
xclass.ngrid.filter <- filter(xclass.ngrid, pix == 500 & p == 0.5)

xc.sp.p1 = xclass.sp[xclass.sp$sim < 2/3 & xclass.sp$p == 1,]
xc.sp.p.5 = xclass.sp[xclass.sp$sim > 2/3 & xclass.sp$p == 0.5,]

scratch <- xclass.sp %>%
  filter(p == 0.5) %>%
  filter(xc %in% c("c-c", "c-t")) %>%
  filter(pix == 500) %>%
  pivot_wider(names_from = xc, values_from = Ngrid)  
t.test(scratch$`c-c`, scratch$`c-t`)

vplot <- xclass.ngrid.filter %>%
  filter(p == 0.5) %>%
  filter(xc %in% c("c-c", "c-t")) %>%
  filter(pix == 500) %>%
  mutate(class = case_when(xc == "c-c" ~ "Correct",
                           xc == "c-t" ~ "Incorrect")) %>%
  mutate(correct = case_when(xc == "c-c" ~ 1,
                             xc == "c-t" ~ 0))
# vplot <- read.csv("Results/EXP4/d-g2_hp-0.5/vplot.csv", header = TRUE)

ggplot(vplot, aes(as.factor(class), Ngrid)) + geom_violin(linetype = "blank", scale ="count", aes(fill = factor(vplot$xc))) +
  theme_classic() + 
  scale_fill_manual(values=c("black", "grey"), labels = c("Correct", "Incorrect")) +
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + 
  xlab("Core species classification") + ylab("Abundance") + 
  theme(axis.title.x=element_text(size=30, vjust = -.2),axis.title.y=element_text(size=30, vjust = 1), legend.text=element_blank(), legend.title=element_blank())
# ggsave(violin, "Results/Plots/EXP4/Fig6.pdf")

#### abun mods ####
logitMod <- glm(correct ~ log10(Ngrid), data=vplot, family=binomial(link="logit"))
summary(logitMod)
# plot(log10(vplot$Ngrid), vplot$correct)
# add fit of modelled curve (s), telling us at what abun level transition is occurring
# figur eout how to calc from parameters (4-4.5 log)
pred.val <- predict(logitMod, type="response")
p <- 0.5
x <- (log(p/(1-p)) - coef(logitMod)[1]) / coef(logitMod)[2]

ggplot(aes((log10(Ngrid/max(Ngrid))), correct), data = vplot) +
  geom_jitter(pch = 16, size = 2, height = 0.025, alpha = 0.6) +
  geom_smooth(data = vplot, aes((log10(Ngrid/max(Ngrid))), correct), method = "glm", method.args = list(family = "binomial")) +
  geom_vline(xintercept = log10((10^x)/max(vplot$Ngrid)), lwd = 1.25, linetype = "dashed") +
  theme_classic() + 
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + 
  scale_x_continuous(breaks = log10(c(0.01, 0.03, 0.1, 0.3, 1)), label = c("1%", "3%", "10%", "30%", "100%")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), label = c("0%", "25%", "50%", "75%", "100%")) +
  xlab("Relative abundance") + ylab("Probability of correct classification") + 
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_blank(), legend.title=element_blank())
ggsave("Results/Plots/EXP4/g2/Fig6.pdf", height = 10, width = 12)
#### Figure S2 ####
ncore = xclass.sp %>%
  filter(xc %in% c("c-c", "t-c", "c-int")) %>%
  group_by(p, run, hp, xc) %>%
  summarise(total_id_core = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  spread(xc, total_id_core) %>%
  mutate(`t-c_0` = ifelse(is.na(`t-c`), 0, `t-c`),
         `c-c_0` = ifelse(is.na(`c-c`), 0, `c-c`),
         `c-int_0` = ifelse(is.na(`c-int`), 0, `c-int`)) %>%
  mutate(perc_ncore1 = `t-c_0`/(`c-c_0`+`t-c_0`),
         perc_ncore = if_else(is.nan(perc_ncore1), 0, perc_ncore1)) %>% 
  ggplot(aes(p, landscape_bin)) +  geom_tile(aes(fill = perc_ncore*100)) +
  theme_classic() + 
  scale_fill_gradient2(midpoint = 50, 
                       name = "Core species: \n% incorrect",low = "#2166ac", mid = "#fddbc7", high ="#b2182b") +
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 8)) +
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + xlab("Detection probability") + ylab("Landscape similarity") +theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25))

ntrans = xclass.sp %>%
  filter(xc %in% c("t-t", "c-t", "t-int")) %>%
  group_by(p, run, hp, xc) %>%
  summarise(total_id_trans = sum(count_xclass)) %>%
  left_join(landscape_bins, by = c("run", "hp", "p")) %>%
  group_by(p, landscape_bin) %>%
  spread(xc, total_id_trans) %>%
  mutate(`c-t_0` = ifelse(is.na(`c-t`), 0, `c-t`),
         `t-t_0` = ifelse(is.na(`t-t`), 0, `t-t`),
         `t-int_0` = ifelse(is.na(`t-int`), 0, `t-int`)) %>%
  mutate(perc_ntrans = `c-t_0`/(`t-t_0`+`c-t_0`)) %>% 
  ggplot(aes(p, landscape_bin)) +  geom_tile(aes(fill = perc_ntrans*100)) +
  theme_classic() + 
  scale_fill_gradient2(midpoint = 50, 
                       name = "Transient species: \n% incorrect",low = "#2166ac", mid = "#fddbc7", high ="#b2182b") +
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 8)) +
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + xlab("Detection probability") + ylab("Landscape similarity") +theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25))


theme_set(theme_cowplot(font_size=20,font_family = "URWHelvetica"))
grid <- plot_grid(ncore + theme(legend.position="right"),
                  ntrans + theme(legend.position="right"),
                  align = 'hv',
                  labels = c("A","B"),
                  label_size = 30,
                  nrow = 2) 
ggsave("Results/Plots/EXP4/Fig_S2_ns_core_trans.pdf", height = 17, width = 14) 


