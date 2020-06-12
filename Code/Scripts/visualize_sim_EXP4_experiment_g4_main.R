## This script is used to visualize simulation summaries from EXP 3
#   (specifically where dispersal is set as a Gaussian (99% within 2 pixels),
#   and immigration to any pixel from the metacommunity is 0.001.)

# Set options and load libraries
options(stringsAsFactors=F)
library(CTSim)
library(dplyr)
library(tidyr)
library(abind)
library(reshape2)
library(raster)

# Define working directory
sum_dir = 'Summaries/EXP4'
fig_dir = 'Results/Plots/EXP4'
sim_dir = 'Code'
data_dir = 'Z:/Snell/core-transient-simulation-LL19/Results/EXP4/d-g4_hp-0.5'

# Load simulation functions
#source(file.path(sim_dir, 'simulation_functions.R'))

setwd("Z:/Snell/core-transient-simulation-LL19")
##########################
## Functions


# A function that extracts parameter values from a runID
#	runID = a string in the filename that gives the parameters
#	assoc_str = character pairing a parameter name with its value
#	sep_str = character separating different parameter-value pairs
get_parms = function(runID, assoc_str='-', sep_str='_'){
  parm_pairs = sapply(strsplit(runID, sep_str), function(x) strsplit(x, assoc_str))
  
  vals = sapply(parm_pairs, function(x) x[2])
  vals = as.data.frame(t(vals))
  names(vals) = sapply(parm_pairs, function(x) x[1])
  
  vals
}


make_plot = function(xlim, ylim, xlab=NULL, ylab=NULL, cex=1){	
  plot.new()
  plot.window(xlim=xlim, ylim=ylim)
  axis(1)
  abline(h=par('usr')[3], lwd=3)
  axis(2, las=1)
  abline(v=par('usr')[1], lwd=3)
  if(!is.null(xlab)) mtext(xlab, 1, 2.5, cex=cex)
  if(!is.null(ylab)) mtext(ylab, 2, 3, cex=cex)
  
}


# Function to plot core/transient dynamics for a pixel
pixDyn = function(results, this_land, row, col, lab = NULL, timewindow = NULL, scale = 3) {
  
  if(this_land[row, col] == 1) { hab = 'A' } else { hab = 'B' }
  
  if(!is.null(lab)) { lab = paste(lab, "; ", sep = '') }
  
  # Fraction of identical landscape over (2*scale+1)x(2*scale+1) region
  het = sum(this_land[max(row-scale, 1):min(row+scale, 32), max(col-scale, 1):min(col+scale, 32)] == this_land[row, col])/
    length(this_land[max(row-scale, 1):min(row+scale, 32), max(col-scale, 1):min(col+scale, 32)])
  
  # Species #1-20 are by definition core in Habitat B; species 21-40 in Habitat A
  if (hab == 'B') {
    core = unlist(lapply(results[row, col, ], function(x) sum(unique(x) <= 20)))
    tran = unlist(lapply(results[row, col, ], function(x) sum(unique(x) > 20)))
  } else {
    core = unlist(lapply(results[row, col, ], function(x) sum(unique(x) > 20)))
    tran = unlist(lapply(results[row, col, ], function(x) sum(unique(x) <= 20)))
  }

  plot(core, type = 'l', xlab = 'Time', ylab = 'Number of species', col = 'red',
       frame.plot=FALSE,
       main = paste(lab, "Landscape similarity ", round(het,2), ";\ncore (blue), transient (red)", sep = ''), 
       lwd = 2, ylim = c(0, max(c(core, tran))))
  points(tran, type = 'l', col = 'dark gray', lwd = 2)
  # vec = c("Core", "Transient")
  # legend('center', legend = vec, lty=1,lwd=6,col = c("red", "dark gray"), cex = 1.5, bty = "n")
  
  # Optionally plot % transient species aggregated over timewindow
  # if(!is.null(timewindow)) {
  #   pct.trans = c()
  #   times = timewindow*1:floor(200/timewindow) - timewindow/2
  #   for (t in 1:floor(200/timewindow)) {
  #     uniqsp = unique(unlist(results[row, col, ((t-1)*timewindow + 2):(t*timewindow+1)]))
  #     if (hab == 'A') {
  #       pct.trans = c(pct.trans, sum(uniqsp <= 20)/length(uniqsp))
  #     } else {
  #       pct.trans = c(pct.trans, sum(uniqsp > 20)/length(uniqsp))
  #     }
  #   }
  #   par(new=T)
  #   plot(times, pct.trans, xlim = c(0,200), ylim = c(0,1), xlab = '', ylab = '', yaxt = 'n',
  #        xaxt = 'n', type = 'l', frame.plot=FALSE)
  #   axis(4, at = seq(0,1, by = 0.25), tcl = .3, labels = F)
  #   mtext(c(1.0, 0.5, 0), 4, at = c(1, .5, 0), cex = .75)
  # }
}

# Plot occupancy histogram for the final time window 162:201 (time 161:200)
pixOccHist = function(results, row, col, lab, timewindow = 15, binwidth = 4) {
  
  tmp = results[row, col, (202 - timewindow):201]
  unq = lapply(tmp, function(x) unique(x))
  occs = table(unlist(unq))
  
  # Define histogram breaks
  xrange = seq(0, timewindow, binwidth)
  
  # Split into two groups (species ID <= or > 20)
  ocore = occs[as.numeric(names(occs)) <= 20]
  otran = occs[as.numeric(names(occs)) > 20]
  
  # compute the counts per interval
  hv1 = hist(ocore,breaks=xrange,plot=F)$counts
  hv2 = hist(otran,breaks=xrange,plot=F)$counts
  
  # Generate a a stacked histogram
  barplot(rbind(hv1,hv2), col=c("#FF0000FF", "#FFFFBFFF"), names.arg = xrange[-1]/timewindow,
          space = 0, las = 1, ylab = "Number of species", xlab = "Occupancy",
          main = paste("Pixel", lab))
}


# Generic (fixed pixel location) dynamics for any simulation run
# data_dir : directory where raw simulation results are stored
# sim      : sim name specifying parameter combinations, e.g. hp-0.9
# run      : simulation run #
# plot_dir : directory for saving output plot
pixelSummary = function(data_dir, sim, run=1:20, plot_dir, plot.pdf = TRUE) {
  
  if (plot.pdf) {
    pdf(paste(plot_dir, '/', sim, '_run', run[1],'-',run[length(run)], '_dynamics.pdf', sep = ''), 
        height = 10, width = 8)
  }
  
  par(mfrow = c(4,3))
  
  for (r in run) {
    suppressWarnings(rm(list = c('results', 'res', 'this_land', 'this_metacomm', 'this_species')))
    load(paste(data_dir, '/', sim, '_run', r, '.Rdata', sep = ''))
    
    # Results from early sims (pre-turnover) are in slightly different structure
    if (class(results) == 'array') {
      res = results
    } else {
      res = results$sim
    }
    
    image(this_land, main = paste(sim, '_run', r, sep = ''))
    
    # Gridded pixels for investigation, remove if d/n want labels
     sites = data.frame(id = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'),
                       row = c(5, 5, 5, 16, 16, 16, 16, 28, 28, 28, 28),
                       col = c(5, 16, 27, 4, 12, 20, 28, 4, 12, 20, 28))

    text(sites$col, 33-sites$row, sites$id, cex = .5)
    
    sapply(1:nrow(sites), function(x) pixDyn(res, this_land, sites$row[x], sites$col[x], sites$id[x], timewindow = 40))
    
    image(this_land, main = paste(sim, '_run', r, sep = ''))
    # text(sites$col, 33-sites$row, sites$id, cex = .5)
    sapply(1:nrow(sites), function(x) 
      pixOccHist(res, sites$row[x], sites$col[x], sites$id[x], timewindow = 40, binwidth = 4))
    
  }
  
  if (plot.pdf) {
    dev.off()
  }
}

# Calculate similarity of landscape to focal pixel
# (Fraction of identical landscape over (2*scale+1)x(2*scale+1) region)

#   loc.xy is a matrix of x and y coordinates on the landscape grid
#   land is the landscape grid raster
#   scale is the (square) radius over which similarity is calculated
land_similarity = function(loc.xy, land, scale) {
  het = apply(loc.xy, 1, function(i) {
    sum(land[max(i[1]-scale, 1):min(i[1]+scale, 32), max(i[2]-scale, 1):min(i[2]+scale, 32)] == land[i[1], i[2]])/
      length(land[max(i[1]-scale, 1):min(i[1]+scale, 32), max(i[2]-scale, 1):min(i[2]+scale, 32)])
  })
  return(het)
}


# Generic pixel cross-classification analysis for all (non-edge) pixels 
# data_dir : directory where raw simulation results are stored
# sim      : sim name specifying parameter combinations, e.g. hp-0.9
# run      : simulation run #
# scale    : radius in pixels over which landscape heterogeneity is calculated
pixelXclass = function(data_dir, sim, run=1:50, scale = 3, t_window = 186:200, 
                       ct_threshold = 1/3, return = 'percent') {
  
  timewindow = length(t_window)
  
  suppressWarnings(rm(list = c('results', 'res', 'this_land', 'this_metacomm', 'this_species')))
  load(paste(data_dir, '/', sim, '_run', run, '.Rdata', sep = ''))
  
  # Results from early sims (pre-turnover) are in slightly different structure
  if (class(results) == 'array') {
    res = results
  } else {
    res = results$sim
  }
  
  # Check that required objects exist
  if(!exists('this_species')|!exists('this_land')|!exists('this_gsad')){
    stop(paste(sim, 'may not contain this_species, this_land, or this_gsad. Summary NOT run.'))		
  } else {
    # Assign objects
    species = this_species
    land = this_land
    gsad = this_gsad
  }
  
  # Number of species
  N_S = dim(species)[1]
  
  # Define different scales of spatial aggregation (in a partition of the grid)
  # Must supply grid dimensions
  scale_locs = sapply(2^c(0:4), function(fact) aggregate_cells(X=c(32,32), dX=fact, dY=fact, form='partition'))
  
  # Locations to be aggregated and evaluated
  locs = scale_locs[[1]]
  
  # Calculate species abundance profiles at the spatial and temporal resolution given by locs and t_window
  abuns_act = calc_abun_profile(locs, t_window, res, N_S)
  
  # Range of detection probabilities to examine
  P_obs = seq(0.1, 1, by = .1)
  
  # Apply observation bias
  abuns_obs = sapply(P_obs, function(p){
    sample_sim(abuns_act, probs = p, return='abundance')
  }, simplify='array')
  dimnames(abuns_obs)[[2]] = 1:N_S # Name columns with species names
  # dims are now [timepoint, species, spatial unit, P]
  
  # Calculate species occupancy
  occ = sapply(1:length(P_obs), function(i){
    use_abun = abuns_obs[,,,i, drop=FALSE]
    use_abun = abind::adrop(use_abun, 4) 
    calc_occupancy(abuns=use_abun, agg_times=NULL, do_freq=F)
  }, simplify='array')
  # dims are now: [spatial unit, species, P]
  
  breaks = seq(ct_threshold, 0.9999, by = ct_threshold)
  
  # Get species birth rates
  b_rates = species[,,'b']
  
  # Get habitat types for spatial units
  habitats = sapply(locs, function(x) average_habitat(x, land))
  
  # Calculate classification of species based on birth rates
  cores = t(sapply(habitats, function(h) b_rates[,h]>0))
  classification = apply(cores, 1:2, function(x) ifelse(x, 'core', 'trans'))	
  
  # Calculate proportion mis-classified
  xclass = sapply(1:length(P_obs), function(i){
    use_occ = occ[,,i, drop=FALSE]
    use_occ = abind::adrop(use_occ, 3)
    tabs = cross_classify(use_occ, breaks, classification=classification, do_each=T, return='counts')
    reshape2::acast(reshape2::melt(tabs, varnames=c('bio','occ','sp_unit')), sp_unit ~ bio+occ)
  }, simplify='array')
  # dims are now [spatial unit, category, P]
  
  xclass.pct = array(dim = dim(xclass))
  for (d in 1:dim(xclass)[1]) {
    co = xclass[d,1:2,] #biologically core spp classified as either core or transient
    tr = xclass[d,3:4,] #biologically transient spp classified as either core or transient
    corepct = sapply(1:ncol(co), function(i) {
      if (colSums(co)[i] == 0) { 
        c(NA, NA) 
      } else { 
        co[,i]/matrix(colSums(co)[i], nrow=1)
      }
    })
    
    tranpct = sapply(1:ncol(tr), function(i) {
      if (colSums(tr)[i] == 0) { 
        c(NA, NA) 
      } else { 
        tr[,i]/matrix(colSums(tr)[i], nrow=1)
      }
    })
    
    xclass.pct[d,,] = rbind(corepct, tranpct)
  }
  
  loc.xy = matrix(unlist(locs), ncol = 2, byrow = TRUE)
  
  landsim = land_similarity(loc.xy, land, scale)
  
  landscapeS = data.frame(x = loc.xy[,1], y = loc.xy[,2], sim = landsim)
  
  if(return == 'count') {
    return(list(xclass=xclass, landSim = landscapeS))
  } else if (return == 'percent') {
    return(list(xclass = xclass.pct, landSim = landscapeS))
  }
}



# Generic pixel cross-classification analysis for all (non-edge) pixels 
# data_dir : directory where raw simulation results are stored
# sim      : sim name specifying parameter combinations, e.g. hp-0.9
# run      : simulation run #
# scale    : radius in pixels over which landscape heterogeneity is calculated
pixelXclassBySpecies = function(data_dir, sim, run=1, scale = 3, t_window = 186:200, 
                       ct_threshold = 1/3) {
  
  timewindow = length(t_window)
  
  suppressWarnings(rm(list = c('results', 'res', 'this_land', 'this_metacomm', 'this_species')))
  load(paste(data_dir, '/', sim, '_run', run, '.Rdata', sep = ''))
  
  # Results from early sims (pre-turnover) are in slightly different structure
  if (class(results) == 'array') {
    res = results
  } else {
    res = results$sim
  }
  
  # Check that required objects exist
  if(!exists('this_species')|!exists('this_land')|!exists('this_gsad')){
    stop(paste(sim, 'may not contain this_species, this_land, or this_gsad. Summary NOT run.'))		
  } else {
    # Assign objects
    species = this_species
    land = this_land
    gsad = this_gsad
  }
  
  # Number of species
  N_S = dim(species)[1]
  
  # Define different scales of spatial aggregation (in a partition of the grid)
  # Must supply grid dimensions
  scale_locs = sapply(2^c(0:4), function(fact) aggregate_cells(X=c(32,32), dX=fact, dY=fact, form='partition'))
  
  # Locations to be aggregated and evaluated
  locs = scale_locs[[1]]
  
  # Calculate species abundance profiles at the spatial and temporal resolution given by locs and t_window
  abuns_act = calc_abun_profile(locs, t_window, res, N_S)
  
  # Get total abundance of each species across the entire grid (remove 1st col describing empty cells)
  sad.grid = data.frame(sp = 1:40, Ngrid = apply(abuns_act[, -1, ], 2, sum))
  
  # Get abundance of each species in each pixel averaged over timewindow
  sad.cell = t(apply(abuns_act[, -1, ], 2:3, sum)) %>% data.frame() %>%
    mutate(pix = 1:1024) %>%
    gather("sp", "Ncell", 1:40) %>%
    mutate(sp = as.numeric(substr(sp, 2, nchar(sp)))) 
           
  
  # Range of detection probabilities to examine
  P_obs = seq(0.1, 1, by = .1)
  
  # Apply observation bias
  abuns_obs = sapply(P_obs, function(p){
    sample_sim(abuns_act, probs = p, return='abundance')
  }, simplify='array')
  dimnames(abuns_obs)[[2]] = 1:N_S # Name columns with species names
  # dims are now [timepoint, species, spatial unit, P]
  
  # Calculate species occupancy
  occ = sapply(1:length(P_obs), function(i){
    use_abun = abuns_obs[,,,i, drop=FALSE]
    use_abun = abind::adrop(use_abun, 4) 
    calc_occupancy(abuns=use_abun, agg_times=NULL, do_freq=F)
  }, simplify='array')
  # dims are now: [spatial unit, species, P]
  
  # Flatten occupancy data
  occ2 = apply(occ, 2, I) %>% data.frame() %>%
    mutate(pix = rep(1:1024, times = 10),
           p = rep(seq(0.1, 1, 0.1), each = 1024)) %>%
    gather("sp", "occ", 1:40) %>%
    mutate(sp = as.numeric(substr(sp, 2, nchar(sp)))) 
  
  breaks = seq(ct_threshold, 0.9999, by = ct_threshold)
  
  # Get species birth rates
  b_rates = species[,,'b']
  
  # Get habitat types for spatial units
  habitats = sapply(locs, function(x) average_habitat(x, land))
  
  # Calculate classification of species based on birth rates
  cores = t(sapply(habitats, function(h) b_rates[,h]>0))
  classification = apply(cores, 1:2, function(x) ifelse(x, 'core', 'trans'))	
  
  # compare classification to occupancy-based designation
  xclass.sp = sapply(1:length(P_obs), function(i){
    use_occ = occ[,,i, drop=FALSE]
    use_occ = abind::adrop(use_occ, 3)
    xclass = classification
    xclass[classification == 'core' & use_occ > (1 - ct_threshold)] = 'c-c'
    xclass[classification == 'core' & use_occ <= ct_threshold & use_occ > 0] = 'c-t'
    xclass[classification == 'trans' & use_occ > (1 - ct_threshold)] = 't-c'
    xclass[classification == 'trans' & use_occ <= ct_threshold & use_occ > 0] = 't-t'
    xclass[xclass == 'core'] = 'c-int'  # biologically core but intermediate occupancy
    xclass[xclass =='trans'] = 't-int'  # biologically transient but intermediate occupancy
    xclass[use_occ == 0] = NA           # absent over the time period 
    return(xclass)
  }, simplify='array')
  
  
  # Calculate landscape similarity to focal pixel
  loc.xy = matrix(unlist(locs), ncol = 2, byrow = TRUE)
  landsim = land_similarity(loc.xy, land, scale)
  landscapeS = data.frame(pix = 1:1024, x = loc.xy[,1], y = loc.xy[,2], sim = landsim)
  # xclass.sp flattened out into 2 dimensions 
  #   (rows increase in pixel.id first, then species, then P_obs)
  #   and then adding in landscape sim, abundances, etc
  xc2 = apply(xclass.sp, 2, I) %>% data.frame() %>% 
    mutate(pix = rep(1:1024, times = 10), 
           p = rep(P_obs, each = 1024)) %>% 
    gather("sp", "xc", 1:40) %>%
    mutate(sp = as.numeric(substr(sp, 2, nchar(sp)))) %>%
    full_join(landscapeS, by = "pix") %>%
    filter(x > scale & x < max(x) - scale & y > scale & y < max(y) - scale, 
           !is.na(xc)) %>%
    inner_join(sad.cell, by = c("pix", "sp")) %>%
    inner_join(sad.grid, by = "sp") %>%
    inner_join(occ2, by = c("pix", "p", "sp")) %>%
    arrange(pix, p, sp) %>%
    dplyr::select(pix, sim, p, sp, Ncell, Ngrid, occ, xc)
    
  return(xclass = xc2)
}




#---------------------------------------------------------------------------------
xclass.out = c()
land.out = c()
xclass.ct = c()
scale = 3
data_dirs =c("Results/EXP4/d-g4_hp-0.5","Results/EXP4/d-g4_hp-0.6", "Results/EXP4/d-g4_hp-0.7", "Results/EXP4/d-g4_hp-0.8", "Results/EXP4/d-g4_hp-0.9")
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

save(xclass.out, xclass.ct, land.out, file = 'Results/Summary/EXP4/g4/pixel_xclass_summary.Rdata')

# Conduct cross-classification by species and pixel
#   (generates ~7M rows, may have memory issues)
# selected 496 bc it is midpoints between 480 and 512, the middle two rows of landscape
pix_all = c(448:544)
xclass.sp = c()
data_dirs =c("Results/EXP4/d-g4_hp-0.5") # ,"Results/EXP4/d-g4_hp-0.6", "Results/EXP4/d-g4_hp-0.7", "Results/EXP4/d-g4_hp-0.8", "Results/EXP4/d-g4_hp-0.9"
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

save(xclass.sp, file = 'Results/Summary/EXP4/g4/pixel_xclass_summary_bysp_w_runs_1_50.Rdata')


####
xclass.ngrid = c()
for(r in 1:50){
tmp = pixelXclassBySpecies("Results/EXP4/g4/d-g4_hp-0.5", "d-g4_hp-0.5", run=r, scale = 3, t_window = 186:200, ct_threshold = 1/3)
xclass.ngrid = rbind(xclass.ngrid, tmp)
}
####

##### pixel summary ######
# CAN START FROM HERE ONCE THE ABOVE HAS BEEN RUN ONCE:
# Load
load('Results/Summary/EXP4/g4/pixel_xclass_summary_bysp_w_runs_1_50.Rdata')
load('Results/Summary/EXP4/g4/pixel_xclass_summary.Rdata')


pixelSummary("Results/EXP4/d-g4_hp-0.5", "d-g4_hp-0.5", run=1:20, "Results/Plots/EXP4", plot.pdf = TRUE)

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
  scale_fill_gradient2(midpoint = 10, 
       name = "Core species",low = "#2166ac", mid = "#fddbc7", high ="#b2182b", breaks = c(9,14,19)) +
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
ggsave("Results/Plots/EXP4/Fig4_t_core_trans.pdf", height = 17, width = 14)

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
  # geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) +
  scale_y_continuous(limits = c(0, 20)) +
  scale_shape_manual(values = c(15, 16)) +
  xlab("Detection probability") + ylab("Mean count of core species") +
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), axis.text.x=element_text(size=25),axis.text.y=element_text(size=25), legend.text=element_text(size=30, hjust = 1, vjust = 0.5), legend.title=element_text(size = 30)) 

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
  scale_y_continuous(limits = c(0, 20)) +
  # geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) +
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
  geom_line(lwd = 1.5, aes(lty = `Detection`)) +
  scale_y_continuous(limits = c(0, 20)) +
  scale_shape_manual(values = c(15, 16))+
  xlab("Landscape similarity") + ylab("Mean count of core species") +
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), axis.text.x=element_text(size=25),axis.text.y=element_text(size=25), legend.text=element_text(size=30, hjust = 1, vjust = 0.5), legend.title=element_text(size = 30)) 

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
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), axis.text.x=element_text(size=25),axis.text.y=element_text(size=25), legend.text=element_text(size=30, hjust = 1, vjust = 0.5), legend.title=element_text(size = 30)) 

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
  ggsave("Results/Plots/EXP4/Fig4_lines.pdf", height = 14, width = 18)  
  
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
ggsave("Results/Plots/EXP4/Fig4_grid_lines.pdf", height = 14, width = 30)  

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
  scale_fill_gradient2(midpoint = 45, 
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
  scale_fill_gradient2(midpoint = 50, 
       name = "Transient species: \n% incorrect",low = "#2166ac", mid = "#fddbc7", high ="#b2182b", limits = c(0, 100)) +
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
ggsave("Results/Plots/EXP4/Fig5_n_core_trans.pdf", height = 18, width = 14)  


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
  scale_y_continuous(limits = c(0, 0.8)) +
  scale_shape_manual(values = c(15, 16)) +
  xlab("Detection probability") + ylab("Core species: % incorrect") +
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30),axis.text.x=element_text(size=25),axis.text.y=element_text(size=25), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 20)) 

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
  scale_y_continuous(limits = c(0, 0.8)) +
  scale_shape_manual(values = c(15, 16)) +
  xlab("Detection probability") + ylab("Transient species: % incorrect") +
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30),axis.text.x=element_text(size=25),axis.text.y=element_text(size=25), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 20)) 

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
  scale_y_continuous(limits = c(0, 0.8)) +
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
  scale_y_continuous(limits = c(0, 0.8)) +
  geom_line(lwd = 1.5,aes(lty = `Detection`)) +
  scale_shape_manual(values = c(15, 16)) +
  xlab("Landscape similarity") + ylab("Transient species: % incorrect") +
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), axis.text.x=element_text(size=25),axis.text.y=element_text(size=25), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 20)) 


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
ggsave("Results/Plots/EXP4/Fig5_lines.pdf", height = 14, width = 18)  

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
ggsave("Results/Plots/EXP4/Fig5_grid_lines.pdf", height = 14, width = 30)  

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
load('Results/Summary/EXP4/d-g4_hp-0.5/pixel_xclass_summary_bysp.Rdata')

xc.sp.p1 = xclass.sp[xclass.sp$sim < 2/3 & xclass.sp$p == 1,]
xc.sp.p.5 = xclass.sp[xclass.sp$sim > 2/3 & xclass.sp$p == 0.5,]

scratch <- xclass.sp %>%
  filter(p == 0.5) %>%
  filter(xc %in% c("c-c", "c-t")) %>%
  filter(pix == 500) %>%
  pivot_wider(names_from = xc, values_from = Ngrid)  
t.test(scratch$`c-c`, scratch$`c-t`)

vplot <- xclass.sp %>%
  filter(p == 0.5) %>%
  filter(xc %in% c("c-c", "c-t")) %>%
  filter(pix == 500) %>%
  mutate(class = case_when(xc == "c-c" ~ "Correct",
                           xc == "c-t" ~ "Incorrect")) %>%
  mutate(correct = case_when(xc == "c-c" ~ 1,
                             xc == "c-t" ~ 0))


ggplot(aes((log10(Ngrid/max(Ngrid))), correct), data = vplot) +
  geom_jitter(pch = 16, size = 2, height = 0.025, alpha = 0.6) +
  geom_smooth(data = vplot, aes((log10(Ngrid/max(Ngrid))), correct), method = "glm", method.args = list(family = "binomial")) +
 geom_vline(xintercept = log10(18053.63/max(vplot$Ngrid)), lwd = 1.25, linetype = "dashed") +
  theme_classic() + 
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + 
  scale_x_continuous(breaks = log10(c(0.01, 0.03, 0.1, 0.3, 1)), label = c("1%", "3%", "10%", "30%", "100%")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), label = c("0%", "25%", "50%", "75%", "100%")) +
  xlab("Relative abundance") + ylab("Probability of correct classification") + 
  theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_blank(), legend.title=element_blank())


# ggsave(violin, "Results/Plots/EXP4/Fig6.pdf")

#### abun mods ####
logitMod <- glm(correct ~ log10(Ngrid), data=vplot, family=binomial(link="logit"))
summary(logitMod)
# plot(log10(vplot$Ngrid), vplot$correct)
# add fit of modelled curve (s), telling us at what abun level transition is occurring
# figure out how to calc from parameters (4-4.5 log)
pred.val <- predict(logitMod, type="response")
p <- 0.5
x <- (log(p/(1-p)) - coef(logitMod)[1]) / coef(logitMod)[2]

relabundances = seq(0.01, 1.0, by = .01)
# for all species with abundance greater than this threshold, what % of them were misclassified?

# t <- data.frame("rel_abun_bins" = i, "e_rate" = e_rate)
# for(i in relabundances){
#   v = filter(vplot, Ngrid/max(Ngrid) > i)
#   e_rate = length(v$class[v$class == "Incorrect"])/length(v$class[v$class == "Correct"])
#   t = rbind(t, c(i, e_rate))
# }
# t$rel_abun_bins <- as.factor(t$rel_abun_bins)

# t2 <- vplot %>%
#   mutate(rel_abun = Ngrid/max(Ngrid),
#          rel_abun_bins = as.factor(round(rel_abun, digits = 1))) %>%
#   left_join(t, by = "rel_abun_bins")

# ggplot(aes(rel_abun_bins, e_rate), data = t2) +
#   geom_jitter(pch = 16, size = 2, width = 0.4, height = 0.01, alpha = 0.6) +
#   theme_classic() +
#   theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) +
#   # scale_x_continuous(breaks = log10(c(0.01, 0.03, 0.1, 0.3, 1)), label = c("1%", "3%", "10%", "30%", "100%")) +
#   scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3), label =  c("0%", "10%",  "20%", "30%")) +
#   xlab("Relative abundance") + ylab("Error rate") +
#   theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_blank(), legend.title=element_blank())
# 
# n_correct <- vplot %>%
#   mutate(rel_abun = (Ngrid/max(Ngrid))*20,
#          rel_bins =round(rel_abun, digits = 0),
#          rel_abun_bins = rel_bins/20) %>%
#   group_by(rel_abun_bins) %>%
#   count(correct) %>%
#   filter(correct == 1)
# 
# t2 <- vplot %>%
#   mutate(rel_abun = (Ngrid/max(Ngrid))*20,
#          rel_bins =round(rel_abun, digits = 0),
#          rel_abun_bins = rel_bins/20) %>%
#   group_by(rel_abun_bins) %>%
#   count(correct) %>%
#   group_by(rel_abun_bins) %>%
#   summarise(sum = sum(n))
# 
# t3 <- full_join(n_correct, t2, by = "rel_abun_bins") %>%
#   mutate(perc_correct = n/sum)
# t3$rel_abun_bins_5 = if_else(t3$rel_abun_bins > 0.5, 0.5, t3$rel_abun_bins)
# t3$rel_abun_bins_5 = as.factor(t3$rel_abun_bins_5)
# 
# ggplot(aes(rel_abun_bins_5, perc_correct), data = t3) +
#   geom_point(size = 4) +
#   theme_classic() +
#   theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) +
#   xlab("Relative abundance") + ylab("Percent correct") +
#   theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_blank(), legend.title=element_blank())

# ( species misclassified whose abundances were > 12%) / (total  species whose abundances were > 12%)

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


###### bimodality #######
pixelXclass_occ = function(data_dir, sim, run=1:50, scale = 3, t_window = 186:200, 
                       ct_threshold = 1/3, return = 'percent') {
  
  timewindow = length(t_window)
  
  suppressWarnings(rm(list = c('results', 'res', 'this_land', 'this_metacomm', 'this_species')))
  load(paste(data_dir, '/', sim, '_run', run, '.Rdata', sep = ''))
  
  # Results from early sims (pre-turnover) are in slightly different structure
  if (class(results) == 'array') {
    res = results
  } else {
    res = results$sim
  }
  
  # Check that required objects exist
  if(!exists('this_species')|!exists('this_land')|!exists('this_gsad')){
    stop(paste(sim, 'may not contain this_species, this_land, or this_gsad. Summary NOT run.'))		
  } else {
    # Assign objects
    species = this_species
    land = this_land
    gsad = this_gsad
  }
  
  # Number of species
  N_S = dim(species)[1]
  
  # Define different scales of spatial aggregation (in a partition of the grid)
  # Must supply grid dimensions
  scale_locs = sapply(2^c(0:4), function(fact) aggregate_cells(X=c(32,32), dX=fact, dY=fact, form='partition'))
  
  # Locations to be aggregated and evaluated
  locs = scale_locs[[1]]
  
  # Calculate species abundance profiles at the spatial and temporal resolution given by locs and t_window
  abuns_act = calc_abun_profile(locs, t_window, res, N_S)
  
  # Range of detection probabilities to examine
  P_obs = seq(0.1, 1, by = .1)
  
  # Apply observation bias
  abuns_obs = sapply(P_obs, function(p){
    sample_sim(abuns_act, probs = p, return='abundance')
  }, simplify='array')
  dimnames(abuns_obs)[[2]] = 1:N_S # Name columns with species names
  # dims are now [timepoint, species, spatial unit, P]
  
  # Calculate species occupancy
  occ = sapply(1:length(P_obs), function(i){
    use_abun = abuns_obs[,,,i, drop=FALSE]
    use_abun = abind::adrop(use_abun, 4) 
    calc_occupancy(abuns=use_abun, agg_times=NULL, do_freq=F)
  }, simplify='array')
  # dims are now: [spatial unit, species, P]
  
  # Flatten occupancy data
  occ2 = apply(occ, 2, I) %>% data.frame() %>%
    mutate(pix = rep(1:1024, times = 10),
           p = rep(seq(0.1, 1, 0.1), each = 1024)) %>%
    gather("sp", "occ", 1:40) %>%
    mutate(sp = as.numeric(substr(sp, 2, nchar(sp)))) 
  
  # Get habitat types for spatial units
  habitats = sapply(locs, function(x) average_habitat(x, land))
  loc.xy = matrix(unlist(locs), ncol = 2, byrow = TRUE)
  landsim = land_similarity(loc.xy, land, scale)
  landscapeS = data.frame(pix = 1:1024, x = loc.xy[,1], y = loc.xy[,2], sim = landsim)
  
  # xclass.sp flattened out into 2 dimensions 
  #   (rows increase in pixel.id first, then species, then P_obs)
  #   and then adding in landscape sim, abundances, etc
  xc2 = full_join(occ2, landscapeS, by = "pix")
}

d4hp5 <- pixelXclass_occ("Results/EXP4/d-g4_hp-0.5", "d-g4_hp-0.5", run = 1, 3, 186:200, 1/3, "percent")
d4hp5$hp = 0.5
d4hp6 <- pixelXclass_occ("Results/EXP4/d-g4_hp-0.6", "d-g4_hp-0.6", run = 1, 3, 186:200, 1/3, "percent")
d4hp6$hp = 0.6
d4hp7 <- pixelXclass_occ("Results/EXP4/d-g4_hp-0.7", "d-g4_hp-0.7", run = 1, 3, 186:200, 1/3, "percent")
d4hp7$hp = 0.7
d4hp8 <- pixelXclass_occ("Results/EXP4/d-g4_hp-0.8", "d-g4_hp-0.8", run = 1, 3, 186:200, 1/3, "percent")
d4hp8$hp = 0.8
d4hp9 <- pixelXclass_occ("Results/EXP4/d-g4_hp-0.9", "d-g4_hp-0.9", run = 1, 3, 186:200, 1/3, "percent")
d4hp9$hp = 0.9

#------------------------------------------------------------------------------------------------------*
# ---- Function for calculating bimodality ----
#======================================================================================================*
# Note 1: Bimodality is the fraction of species occurring at either end of 
# occupancy distribution. We use a randomization approach to test whether the 
# distribution is significantly bimodal.
# Note 2: To run this function the number of time samples for the site (nt) needs
# to be specified. This is done so in the wrapper summary table function.

# In these functions, propOcc refers to a vector of occupancy values for the
# species at a single site, and nTime is the number of time samples (typically
# years) as an integer.


bimodalityFun = function(propOcc_or_RandomPropOcc, nTime){
  occs = propOcc_or_RandomPropOcc
  maxvar = var(c(rep(1/nTime,floor(length(occs)/2)),
                 rep(1,ceiling(length(occs)/2))))
  return(var(occs)/maxvar)
}

# Random sample of occurences for a given site (to be used in randomization, below):

randomOccsFun = function(propOcc, nTime){
  # Generate a table (data frame) of occProps and frequencies:
  occPropTable = data.frame(table(propOcc))
  # Create a data frame of possible occProps:
  occPropDummyTable = data.frame(propOcc = seq(1/nTime, 1, by = 1/nTime))
  # Merge the two data frames:
  combinedTable = merge(occPropDummyTable, occPropTable, all.x = T)
  combinedTable[is.na(combinedTable[,2]),2]<-0                            # Replace NA's with zeros
  # Reassign bin values randomly and add to frame:
  newFreq = sample(combinedTable$Freq, length(combinedTable[,1]))
  randomTable = data.frame(combinedTable[,1], newFreq)
  randomOccs=unlist(apply(randomTable, 1, function(x) rep(x[1], x[2])))
  return(as.vector(randomOccs))
}
# currently running this manually bc not sure how to rename each DF
bimod_input <- bind_rows(d4hp5, d4hp6, d4hp7, d4hp8, d4hp9) %>%
  filter(pix  == 500)

bimod_5 = data.frame(res = c(), hp = c(), bimod = c())
for(h in unique(bimod_input$hp)){
  for(pobs in unique(bimod_input$p)){
    for(i in unique(bimod_input$x)){
      for(j in unique(bimod_input$y)){
  propOcc = filter(bimod_input, x == i & y == j & p == pobs & hp == h) %>%
    dplyr::select(occ) %>%
    na.omit()
  propOcc= as.vector(propOcc$occ)
  bi_sub <- bimodalityFun(propOcc, 15)
  bimod_5 = rbind(bimod_5, c(i, j, pobs, h, bi_sub))
      }
    }
  }
}
names(bimod_5) <- c("x", "y", "p", "hp", "bimod")
bimod_5 <- data.frame(bimod_5)

# d = d4hp5
# d2 = "d4hp5"
# hp = substring(as.character(d2), 5, 5)
# for(i in colnames(d)){
#   propOcc = d[,i]
#   propOcc2 = propOcc[!is.na(propOcc)]
#   bi_sub <- bimodalityFun(propOcc2, 15)
#   bimod_5 = rbind(bimod_5, data.frame(i, hp, bi_sub))
# }


bimod_full <- bimod_5 %>%
  left_join(bimod_input[,c("p", "x", "y", "sim", "hp")], by = c("x", "y", "p", "hp")) %>%
  distinct() %>%
  mutate(landscape_bin = case_when(sim < 0.2 ~ 0.1,
                                   sim >= 0.2 & sim < 0.3 ~ 0.2,
                                   sim >= 0.3 & sim < 0.4 ~ 0.3,
                                   sim >= 0.4 & sim < 0.5 ~ 0.4,
                                   sim >= 0.5 & sim < 0.6 ~ 0.5,
                                   sim >= 0.6 & sim < 0.7 ~ 0.6,
                                   sim >= 0.7 & sim < 0.8 ~ 0.7,
                                   sim >= 0.8 & sim < 0.9 ~ 0.8,
                                   sim >= 0.9 & sim < 1.0 ~ 0.9,
                                   sim == 1 ~ 1)) %>%
  ggplot(aes(p, landscape_bin)) +  geom_tile(aes(fill = bimod)) +
  scale_fill_gradient(name = "Bimodality",low = "white", high = "#b2182b") +
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 8)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=25),axis.text.y=element_text(size=25)) + xlab("Detection probability") + ylab("Proportion landscape similarity") +theme(axis.title.x=element_text(size=30),axis.title.y=element_text(size=30), legend.text=element_text(size=25, hjust = 1, vjust = 0.5), legend.title=element_text(size = 25))
ggsave("Results/Plots/EXP4/bimod.pdf", height = 14, width = 16) 

