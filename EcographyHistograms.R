# Chunk 1 read in adult data and prep adult data for analyses

# set working directory
fp <- 'C:/Users/srecord/Dropbox/NASA/ForestScaling/'
setwd(fp)

# load relevant libraries
library(cowplot)
library(grid)
library(tidyverse)

# read in Harvard Forest FORESTGEO data for stems >1 cm dbh
HFforestgeo <- read.csv("HARVdata/FieldData/ForestGEO/hf253-05-stems-2019.csv", header=TRUE)
dim(HFforestgeo) # 123218

# select live stems from 2019 census
HFforestgeo <- subset(HFforestgeo, df.status=="alive")
HFforestgeo <- subset(HFforestgeo, census.id==2)
dim(HFforestgeo) # 61217

# Remove 15 rows with NA values for 
HFforestgeo <- HFforestgeo[complete.cases(HFforestgeo[,c("dbh")]),]
dim(HFforestgeo)

# Abundance by diameter plots
# Binning and plotting custom functions
logbin <- function(x, y = NULL, n) {
  logx <- log10(x)                                           # log transform x value (biomass or diameter)
  bin_edges <- seq(min(logx), max(logx), length.out = n + 1) # get edges of bins
  logxbin <- rep(NA, length(logx))                           # create data structure to assign trees to bins
  b <- bin_edges                                             # add a little to the biggest bin temporarily
  b[length(b)] <- b[length(b)] + 1                           # (so that the biggest single tree is put in a bin)
  for (i in 1:length(logx)) {
    logxbin[i] <- sum(logx[i] >= b)                          # assign each tree to a bin
  }
  bin_midpoints <- numeric(n)
  for (i in 1:n) {
    bin_midpoints[i] <- mean(10^(bin_edges[i:(i+1)]))        # backtransform bin edges to linear, and get midpoints
  }
  bin_widths <- diff(10^bin_edges)                           # get linear width of each bin
  bin_factor <- factor(logxbin, levels=1:n)                  # convert bin to factor (required to deal with zeroes if present)
  bin_counts <- table(bin_factor)                            # find number of trees in each bin
  if (!is.null(y)) {
    rawy <- tapply(y, bin_factor, sum)                       # sum y value (production) in each bin
    rawy[is.na(rawy)] <- 0                                   # add zeroes back in if present
    bin_values <- as.numeric(rawy/bin_widths)                # divide production by width for each bin 
  }
  else {
    bin_values <- as.numeric(bin_counts/bin_widths)          # 1-dimensional case.
  }
  
  return(data.frame(bin_midpoint = bin_midpoints,            # return result!
                    bin_value = bin_values,                  # also add bin min and max for bar plot purposes
                    bin_count = as.numeric(bin_counts),
                    bin_min = 10^bin_edges[1:n],
                    bin_max = 10^bin_edges[2:(n+1)]))
  
}

# Function for plotting.

plotlogbin_cutoff <- function(dat, xl, yl, plottitle, plotsubtitle=NULL, reg = FALSE, cutoff = NA, y_min=0.1, y_max=53, x_min=1.1, x_max=110, plotarea=35, y_values=-1:3, x_values = 0:2, barfill) {
  
  dat <- transform(dat, bin_value = bin_value/plotarea) # kg per hectare.
  
  p <- ggplot(dat, aes(xmin=bin_min, xmax=bin_max, ymin=0, ymax=bin_value)) + 
    geom_rect(alpha = 0.5, fill=barfill) +
    scale_x_log10(name = xl, expand = c(0,0),
                  breaks = 10^x_values,
                  labels = as.character(10^x_values), limits=c(x_min,x_max)) +
    scale_y_log10(name = yl, expand = c(0,0), limits = c(y_min, y_max),
                  breaks = 10^y_values,
                  labels = as.character(10^y_values)) +
    panel_border(colour = 'black') + 
    ggtitle(plottitle, plotsubtitle) +
    labs(ylab="Log10 Abundance", xlab="Log10 Diameter")
  if (reg) {
    p <- p +
      stat_smooth(method = 'lm', se = FALSE, color = 'forestgreen', size = 2,
                  aes(x = bin_midpoint, y = bin_value)) +
      geom_text(x = -Inf, y = -Inf, 
                label = paste('Slope without cutoff:', 
                              round(lm(I(log10(bin_value)) ~ I(log10(bin_midpoint)), data=dat)$coef[2], 2)),
                hjust = 0, vjust = -1.5)
  
    if (!is.na(cutoff)) {
      p <- p +
        stat_smooth(method = 'lm', se = FALSE, color = 'goldenrod', size = 2,
                    aes(x = bin_midpoint, y = bin_value), data = subset(dat, bin_midpoint <= cutoff)) +
        geom_text(x = -Inf, y = -Inf, 
                  label = paste('Slope with cutoff:', 
                                round(lm(I(log10(bin_value)) ~ I(log10(bin_midpoint)), data=subset(dat, bin_midpoint <= cutoff))$coef[2], 2)),
                  hjust = 0, vjust = -0.2)
    }
  }
  return(p)
}


# Run binning algorithm for density scaling visualization.
numbins <- 15 # Note we ran this for 10, 15, and 50 bins.
# Subset late successional species with slow life histories
slow <- subset(HFforestgeo, subset = sp %in% c('tsugca','fagugr', 'nyssyl','picrub','picabi'))
# Bin late successional species
slow_bin <-  with(slow,logbin(x=dbh, y=NULL, n = numbins))
# Subset early/mid successional species with fast life histories
fast <- subset(HFforestgeo, subset = sp %in% c('querru','acerru','pinust','betpop','betlen','betpap','queral','pinres','pruser','quervel','popgra'))
# Bin early/mid successional species 
fast_bin <-  with(fast,logbin(x=dbh, y=NULL, n = numbins))
# Bin all trees together 
allTrees_bin <-  with(HFforestgeo,logbin(x=dbh, y=NULL, n = numbins))

# Create plots
plotlogbin_cutoff(slow_bin, xl=log10(slow_bin$bin_midpoint)+1, yl=log10(slow_bin$bin_count), plottitle='Late Successional, Slow Trees', barfill="#276F02")
plotlogbin_cutoff(fast_bin, xl=log10(fast_bin$bin_midpoint)+1, yl=log10(fast_bin$bin_count), plottitle='Early/Mid Successional, Fast Trees', barfill="#B8D641")

# Obtain slope estimate of all trees (including fast and slow species)
summary(lm(log10(bin_count+1)~log10(bin_midpoint+1), data=allTrees_bin))
# Obtain slope estimate for fast trees
summary(lm(log10(bin_count+1)~log10(bin_midpoint+1), data=fast_bin))
s
# Obtain slope estimate for fast trees
summary(lm(log10(bin_count+1)~log10(bin_midpoint+1), data=slow_bin))
s

