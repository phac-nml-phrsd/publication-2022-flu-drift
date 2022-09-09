###
###   STATISTICAL ANALYSIS
###

library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(lubridate)
library(stringr)
library(patchwork)

if(!exists('df.final')) load('merged-epi.RData')


# Seasonal lag between severity index and antigenic distance
lag.si.ad = 1

# Exclude seasons that do not have
# enough data points to calculate 
# the severity index or the 
# interseason genetic distance.

dfa = df.final %>% 
  filter(n.inter > 20,  
         si.n >= 3) %>% 
  filter(!is.na(vax.eff))


# --- plot data ---

g =  dfa %>%
  filter(season.lag == lag.si.ad) %>% 
  ggplot(aes(x=dist.inter.m, y=si, 
             fill  = factor(distance.type))) + 
  geom_point(aes(size = vax.eff), shape=21, alpha=.8)+
  ggrepel::geom_text_repel(aes(label=yss), color='black', size=2, force = 70)+
  facet_grid(subtype ~ distance.type) + 
  theme(panel.grid.minor = element_blank())+
  labs(
    x = 'interseason distance',
    y = 'severity index'
  )
g

pdf('plots/fig-data.pdf', width=14)
plot(g)
dev.off()

# This one is prettier for the manuscript:

g.ms = dfa %>%
  # use season lag of 1, exclude full HA1 definition
  filter(season.lag==1, distance.type != "hamming-full.HA1") %>%
  mutate(distance.type.f = factor(distance.type,
                                  levels = c("hamming-narrow", 
                                             "hamming-broad3", 
                                             "hamming-full.HA12"),
                                  labels = c("Narrow", "Broad", "Full"))) %>%
  ggplot(aes(x=dist.inter.m, y=si, 
             fill  = factor(distance.type.f))) + 
  geom_point(aes(size = vax.eff), shape=21, alpha=.8)+
  scale_size_continuous(limits=c(-0.2,1), 
                        breaks=c(0, 0.25, 0.5, 0.75), 
                        labels=c("0%", "25%", "50%", "75%"),
                        #range=c(1,9)
  )+
  ggrepel::geom_text_repel(aes(label=yss), color="gray36", size=2.5, force = 70, 
                           box.padding = 0.6, segment.color = "gray25")+
  facet_grid(subtype ~ distance.type.f) +
  coord_cartesian(ylim = c(-0.75,1.25)) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        strip.text.x = element_text(margin = margin(10,0,15,0),
                                    face = "bold"),
        strip.text.y = element_text(angle = 0, 
                                    margin = margin(0,15,0,15), 
                                    face = "bold"),
        axis.text  = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 8, b = 2)),
        axis.title.y = element_text(margin = margin(l = 2, r = 8)),
        legend.margin = margin(5,5,5,5),
        legend.position = "top",
        legend.justification = "right",
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 12),
        legend.box.background = element_rect(colour = "black"))+
  labs(
    x = 'Interseasonal Distance',
    y = 'Severity Index'
  )+
  guides(fill = "none",
         size = guide_legend("Vaccine Effectiveness"))
g.ms


pdf('plots/fig-data-ms.pdf', width = 12)
plot(g.ms)
dev.off()


# Run all regressions and 
# generate the diagnostic plots

for(st in c('H1N1', 'H3N2', 'B')){
  d = dfa %>% 
    filter(subtype==st, season.lag == lag.si.ad)
  
  m = lm(data = d, formula = 'si ~ dist.inter.m * vax.eff')
  #summary(m)
  
  # diagnostic plots
  pdf(paste0('plots/diag-lm-',st,'.pdf'))
  par(mfrow=c(2,2))
  diagplot =   plot(m, ask = FALSE)
  dev.off()
}

# Generate simulated values centred on estimates
# to run the regression replicates:

nr = nrow(dfa)
n.repl = 100
tmp = list()

for(i in 1:n.repl){
  tmp[[i]] = dfa %>% 
    mutate(d = rnorm(n=nr, mean = dist.inter.m, sd=dist.inter.sd),
           v = rnorm(n=nr, mean = vax.eff, sd= (vax.eff.hi-vax.eff.lo)/2),
           iter = i)
}
dfr = bind_rows(tmp)


#' Run one single linear regression (will be used in a loop)
#' 
#' @param i Intger. Index.
#' @param st String. Subtype of the virus (e.g., "H1N1")
#' @param lag Integer. Season lag between clinical data and genetic data. 
#' @param dfr Dataframe of merged data (clinical and genetic)
#' @param dist.type String. Distance type (e.g., "hamming")
#' 
run_regr <- function(i, st, lag, dfr, dist.type) {
  
  if(i%%25==0) message(paste(i, ' '), appendLF = FALSE)
  
  d = dfr %>% 
    filter(subtype==st, 
           season.lag ==lag, 
           subtype == st, 
           distance.type == dist.type,
           iter == i)
  
  np = nrow(d)  # number of data points
  
  m = lm(data = d, formula = 'si ~ d * v')
  
  # SI = k.d * d + k.v * v + k.dv * d * v + error
  
  s = summary(m)
  
  res = data.frame(
    subtype = st, 
    season.lag = lag, 
    distance.type = dist.type,
    np   = np, 
    # retrieve slopes
    k.d  = s$coefficients[2,1],
    k.v  = s$coefficients[3,1],
    k.dv = s$coefficients[4,1],
    # retrieve p-values
    p.d  = s$coefficients[2,4],
    p.v  = s$coefficients[3,4],
    p.dv = s$coefficients[4,4],
    r2adj = s$adj.r.squared,
    iter = i
  )
  return(res)
}

sub.vec = unique(dfr$subtype)
lag.vec = unique(dfr$season.lag)
dst.vec = unique(dfr$distance.type)


# Loop through all subtypes

q=1 ; z = list()

for(st in sub.vec){
  for(lag in lag.vec){
    for(dist.type in dst.vec){
      
      message(paste('\n',st, lag, dist.type))
      
      z[[q]] = lapply(1:n.repl, run_regr, 
                      st=st, lag=lag, dfr=dfr, 
                      dist.type = dist.type)
      q = q+1
    }
  }
}

res.all = lapply(z, bind_rows) %>% bind_rows()

#' Plot regression results
#' 
plot_reg_results <- function(res.all, st, var.type, lag, do.guides) {
  
  if(0){
    st = 'H1N1'
    var.type = '^k'
    lag = 1
    do.guides = TRUE
  }
  
  df.plot = res.all %>%
    filter(subtype == st) %>% 
    pivot_longer(-c(iter, subtype, season.lag, distance.type)) %>% 
    filter(season.lag == lag) %>%
    filter(name != 'np') %>%
    filter(grepl(var.type, name)) %>% 
    #
    # Select distance
    filter(distance.type != 'hamming-full.HA1') %>%
    #
    # Rename for pretty plot
    mutate(dist.plot = case_when(
      grepl( 'broad$', distance.type)  ~ 'Broad1',
      grepl( 'broad3', distance.type)  ~ 'Broad',
      grepl( 'narrow', distance.type,) ~ 'Narrow',
      grepl( 'misc1', distance.type,)  ~ 'Misc1',
      grepl( 'misc2', distance.type,)  ~ 'Misc2',
      grepl( 'HA12', distance.type,)   ~ 'Full',
      grepl( 'HA1$', distance.type,)   ~ 'Full HA1',
      TRUE ~ distance.type
    )) %>% 
    mutate(name.plot = 
             case_when(
               name == 'k.v'  ~ 'Coeff Vax',
               name == 'k.d'  ~ 'Coeff Dist',
               name == 'k.dv' ~ 'Coeff Dist:Vax',
               name == 'p.v'  ~ 'p-value Vax',
               name == 'p.d'  ~ 'p-value Dist',
               name == 'p.dv' ~ 'p-value Dist:Vax'
             ))
  
  ci = c(0.5, 0.95)
  
  ss = df.plot %>% 
    group_by(subtype, dist.plot, name.plot) %>% 
    summarize(m = mean(value), 
              qlo  = quantile(value, probs = 0.5 - ci[1]/2),
              qhi  = quantile(value, probs = 0.5 + ci[1]/2),
              qvlo = quantile(value, probs = 0.5 - ci[2]/2),
              qvhi = quantile(value, probs = 0.5 + ci[2]/2),
              .groups = 'drop'
              )
  
  # Manual adjustment to display the distance types.
  # For now and for the supplementary material, 
  # only H3N2 has a `Misc2` distace type: 
  if(st != 'H3N2') ss = filter(ss, dist.plot != 'Misc2')
  
  # Create table (CSV file for appendix)
  tabcsv = ss
  names(tabcsv) <- c('Subtype', 'Distance type', 'Parameter',
                     'mean', 
                     '25% quantile',
                     '50% quantile',
                     '2.5% quantile',
                     '97.5% quantile'
                     )
  
  # --- plot
  
  tmp = ifelse(var.type=='^k', 'Coeff. estimates', 'P-values of coeff. estimates')
  thetitle = paste(tmp, st)
  
  # "Natural" reordering of distances
  df.plot$dist.plot <- factor(df.plot$dist.plot, 
                              levels = c('Narrow', 'Broad', 'Full'))
  
  a  = .4
  sz = 3
  
  g = ss %>% 
    ggplot() +
    geom_segment(aes(x=dist.plot, xend=dist.plot, 
                     y=qvlo, yend=qvhi, color=dist.plot), 
                 alpha = a*0.7, size = sz)+
    geom_segment(aes(x=dist.plot, xend=dist.plot, 
                     y=qlo, yend=qhi, color=dist.plot), 
                 alpha = a, size = sz)+
    geom_point(aes(x=dist.plot, y=m, color=dist.plot), size= 0.8*sz,
               shape=21, fill='white', stroke=1.2) +
    facet_wrap(~ name.plot, scales = 'free') +
    theme(panel.grid.minor   = element_blank(), 
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill='grey95'),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(y='', x='', title = thetitle, 
         color = 'Distance type', fill = 'Distance type' )
  
  if(var.type == '^k') 
    g = g + geom_hline(yintercept = 0, color='darkgrey', linetype='dashed')
  
  if(do.guides)
    g = g + guides(color='none', fill='none')
  
  res = list(plot = g, table = tabcsv)
  return(res)
}

# Generate all plots 

for(st in c('H1N1', 'H3N2', 'B')){
  for(lag in c(1,2)){
    
    tmp.k = plot_reg_results(res.all, st, var.type='^k', lag=lag, do.guides = TRUE)
    tmp.p = plot_reg_results(res.all, st, var.type='^p', lag=lag, do.guides = TRUE)
    
    g.k = tmp.k$plot
    g.p = tmp.p$plot
    
    g = wrap_plots(g.k, g.p, ncol=1, guides = 'collect')
    fname  = paste0('plots/fig-results-',st,'-lag',lag,'.pdf')
    pdf(fname, width=14)
    plot(g)
    dev.off()
  }
}

# Figures for manuscript

fig = list() ; j=1
tabcsv = list()

for(st in c('H1N1', 'H3N2', 'B')){
  
  print(st)
  
  do.guides = FALSE
  if(st == 'B') do.guides=TRUE
  
  tmp.k = plot_reg_results(res.all, st, var.type='^k', lag=lag.si.ad, do.guides = do.guides)
  tmp.p = plot_reg_results(res.all, st, var.type='^p', lag=lag.si.ad, do.guides = FALSE)
  
  f.k = tmp.k$plot
  f.p = tmp.p$plot
  
  th = theme(plot.margin = margin(10,25,10,3))
  f.k = f.k + th
  f.p = f.p+ th
  
  fig[[j]] = wrap_plots(f.k, f.p, ncol=1)
  
  # CSV tables
  tabcsv[[j]] = rbind(tmp.k$table, tmp.p$table)
  
  # increment counter
  j = j+1
}

# CSV table saved 
tabcsv2 = do.call('rbind', tabcsv)
write.csv(tabcsv2, 'plots/table_suppl_regression_results.csv', row.names = FALSE)

# plot saved 
fig.final = (wrap_plots(fig, guides = 'collect') & 
  theme(legend.position = 'bottom') ) +
  plot_annotation(caption = 'circle indicates mean estimates across 100 Monte Carlo iterations,\n
                  light and dark vertical bars show 95% and 50% quantiles.')

fig.final

pdf('plots/fig-results.pdf', width = 14, height=7)
plot(fig.final)
dev.off()









