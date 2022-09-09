###
###   STATISTICAL ANALYSIS
###

library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(lubridate)
library(stringr)
library(patchwork)

# https://www.r-bloggers.com/2018/05/regression-with-interaction-terms-how-centering-predictors-influences-main-effects/
# https://www.theanalysisfactor.com/interpreting-interactions-in-regression/


if(!exists('df.final')) load('merged-epi.RData')


# Exclude seasons that do not have
# enough data points to calculate 
# the severity index or the 
# interseason genetic distance.

dfa = df.final %>% 
  filter(n.inter > 20, #10 
         si.n >= 3) %>% # 3
  filter(!is.na(vax.eff))

h1n1 = df.final %>% 
  filter(subtype=='H1N1', season.lag == 1)

g =  dfa %>%
  filter(season.lag==2) %>% 
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

# redo plot with different parameters
# change axes of facet wrap, remove distance HA1
g =  dfa %>%
  filter(season.lag==1, distance.type != "hamming-full.HA1") %>% 
  ggplot(aes(x=dist.inter.m, y=si, 
             fill  = factor(subtype))) + 
  geom_point(aes(size = vax.eff), shape=21, alpha=.8)+
  ggrepel::geom_text_repel(aes(label=yss), color='black', size=2, force = 70)+
  facet_grid(distance.type ~ subtype) + 
  theme(panel.grid.minor = element_blank())+
  labs(
    x = 'interseason distance',
    y = 'severity index'
  )
g

pdf('plots/fig-data-zc.pdf', width=14)
plot(g)
dev.off()


# new plots, keep old orientation
g = dfa %>%
  # use season lag of 1, exclude full HA1 definition
  filter(season.lag==1, distance.type != "hamming-full.HA1") %>%
  mutate(distance.type.f = factor(distance.type,
                                  levels = c("hamming-narrow", 
                                             "hamming-broad", 
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
        strip.text.x = element_text(margin = margin(10,0,15,0)),
        strip.text.y = element_text(angle = 0, 
                                    margin = margin(0,15,0,15), 
                                    face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 8, b = 2)),
        axis.title.y = element_text(margin = margin(l = 2, r = 8)),
        legend.margin = margin(5,5,5,5),
        legend.position = "top",
        legend.justification = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box.background = element_rect(colour = "black"))+
  labs(
    x = 'Interseasonal Distance',
    y = 'Severity Index'
  )+
  guides(fill = "none",
         size = guide_legend("Vaccine Effectiveness"))
g


pdf('plots/fig-data-zc-new.pdf', width = 12)
plot(g)
dev.off()



# try with relative distances
num.residues = data.frame(subtype = c("H1N1", "H3N2", "B"),
                          hamming_narrow = c(50, 131, 78),
                          hamming_broad = c(68, 215, 96),
                          hamming_full.HA1 = c(327, 329, 345),
                          hamming_full.HA12 = c(549, 550, 570))

# transform to long format
num.residues.long = num.residues %>%
  pivot_longer(cols = !subtype, 
               names_to = "distance.type", 
               values_to = "res.count")

# replace "_" with "-" in distance type description
num.residues.long$distance.type = gsub("_", "-", num.residues.long$distance.type)

# new dataframe with relative interseasonal distance
dfa.v2 = dfa %>%
  left_join(num.residues.long, by = c("subtype", "distance.type")) %>%
  mutate(rel.dist.inter.m = dist.inter.m/res.count)

g =  dfa.v2 %>%
  filter(season.lag==1) %>% 
  ggplot(aes(x=rel.dist.inter.m, y=si, 
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

pdf('plots/fig-data-relative-zc.pdf', width=14)
plot(g)
dev.off()


# new plots (formatted) with relative distances
# year labels may need further formatting
g = dfa.v2 %>%
  # use season lag of 1, exclude full HA1 definition
  filter(season.lag==1, distance.type != "hamming-full.HA1") %>%
  mutate(distance.type.f = factor(distance.type,
                                  levels = c("hamming-narrow", 
                                             "hamming-broad", 
                                             "hamming-full.HA12"),
                                  labels = c("Narrow", "Broad", "Full"))) %>%
  ggplot(aes(x=rel.dist.inter.m, y=si, 
             fill  = factor(distance.type.f))) + 
  geom_point(aes(size = vax.eff), shape=21, alpha=.8)+
  scale_size_continuous(limits=c(-0.2,1), 
                        breaks=c(0, 0.25, 0.5, 0.75), 
                        labels=c("0%", "25%", "50%", "75%"))+
  ggrepel::geom_text_repel(aes(label=yss), color="gray36", size=2.5, force = 70, 
                           box.padding = 0.6, segment.color = "gray25")+
  facet_grid(subtype ~ distance.type.f) +
  coord_cartesian(ylim = c(-0.75,1.25)) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        strip.text.x = element_text(margin = margin(10,0,15,0)),
        strip.text.y = element_text(angle = 0, 
                                    margin = margin(0,15,0,15), 
                                    face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 8, b = 2)),
        axis.title.y = element_text(margin = margin(l = 2, r = 8)),
        legend.margin = margin(5,5,5,5),
        legend.position = "top",
        legend.justification = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box.background = element_rect(colour = "black"))+
  labs(
    x = 'Interseasonal Distance',
    y = 'Severity Index'
  )+
  guides(fill = "none",
         size = guide_legend("Vaccine Effectiveness"))
g

pdf('plots/fig-data-relative-zc-new.pdf', width=12)
plot(g)
dev.off()


st = 'H3N2'   # H1N1 H3N2
lag = 1

d = dfa %>% 
  filter(subtype==st, season.lag ==lag)

m = lm(data = d, formula = 'si ~ dist.inter.m * vax.eff')
summary(m)
# diagnistic plots
if(0) plot(m)


# m0 = lm(data = d, formula = 'si ~ dist.inter.m + vax.eff')
# summary(m0)


# different package for OLS
library(rms)
d.dd = datadist(d)
options(datadist="d.dd")

n = ols(si ~ dist.inter.m * vax.eff, data = d)
summary(n)
summary(n, dist.inter.m = c(1,2))

# Plot predicted values for different values of vax eff
plot(Predict(n, dist.inter.m, vax.eff = 0.1))
plot(Predict(n, dist.inter.m, vax.eff = 0.5))
plot(Predict(n, dist.inter.m, vax.eff = 0.9))


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


# st = 'H3N2' ; lag=1; dist.type = 'hamming-narrow'

sub.vec = unique(dfr$subtype)
lag.vec = unique(dfr$season.lag)
dst.vec = unique(dfr$distance.type)

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


plot_data <- function(res.all) {
  
}

plot_reg_results <- function(st, var.type, lag) {
  df.plot = res.all %>%
    filter(subtype == st) %>% 
    pivot_longer(-c(iter, subtype, season.lag, distance.type)) %>% 
    filter(season.lag == lag) %>%
    filter(name != 'np') %>%
    filter(grepl(var.type, name))
  
  
  g = df.plot %>% 
    ggplot() + 
    geom_hline(yintercept = 0, color='darkgrey', linetype='dashed')+
    geom_boxplot(aes(x=name, y=value, color = distance.type), outlier.shape = NA)+
    facet_wrap(~name, scales = 'free') +
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(y='', x='', title = paste(st) )
  return(g)
}


for(st in c('H1N1', 'H3N2', 'B')){
  for(lag in c(1,2)){
    g.k = plot_reg_results(st, var.type='^k', lag=lag)
    g.p = plot_reg_results(st, var.type='^p', lag=lag)
    
    g = wrap_plots(g.k, g.p, ncol=1, guides = 'collect')
    fname  = paste0('plots/fig-results-',st,'-lag',lag,'.pdf')
    pdf(fname, width=14)
    plot(g)
    dev.off()
  }
}



#### Generate table ####

# Function to round values to 2 decimal places
round_num = function(input_num) {
  output_num = format(round(input_num, 2), nsmall = 2)
  return(output_num)
}

# Format results into a table
res.all.sum = res.all %>%
  group_by(subtype, season.lag, distance.type) %>%
  # take mean and SD of coefficients and p-values
  # round all numbers to 2 decimal places
  summarise(k.d.mean = round_num(mean(k.d)),
            k.d.sd = round_num(sd(k.d)),
            p.d.mean = round_num(mean(p.d)),
            p.d.sd = round_num(sd(p.d)),
            k.v.mean = round_num(mean(k.v)),
            k.v.sd = round_num(sd(k.v)),
            p.v.mean = round_num(mean(p.v)),
            p.v.sd = round_num(sd(p.v)),
            k.dv.mean = round_num(mean(k.dv)),
            k.dv.sd = round_num(sd(k.dv)),
            p.dv.mean = round_num(mean(p.dv)),
            p.dv.sd = round_num(sd(p.dv)),
            .groups = "drop") %>%
  # combine coefficients and SDs,
  # recode distance type values so they can be rearranged
  mutate(k.d = paste0(k.d.mean, " (", k.d.sd, ")"),
         p.d = paste0(p.d.mean, " (", p.d.sd, ")"),
         k.v = paste0(k.v.mean, " (", k.v.sd, ")"),
         p.v = paste0(p.v.mean, " (", p.v.sd, ")"),
         k.dv = paste0(k.dv.mean, " (", k.dv.sd, ")"),
         p.dv = paste0(p.dv.mean, " (", p.dv.sd, ")"),
         distance.type = case_when(
           distance.type == "hamming-narrow" ~ "1.hamming.narrow",
           distance.type == "hamming-broad" ~ "2.hamming.broad",
           distance.type == "hamming-full.HA1" ~ "3.hamming.full.HA1",
           distance.type == "hamming-full.HA12" ~ "4.hamming.full.HA12"
         )) %>%
  # select variables to include in table
  select(subtype, season.lag, distance.type, k.d:p.dv) %>%
  # rearrange by distance type
  arrange(distance.type) %>%
  # change to long format
  pivot_longer(cols = k.d:p.dv, names_to = "coeff", 
               values_to = "value") %>%
  # create new column to categorize coefficients
  # create new column for coefficient category, subtype and lag
  # create new column for distance type and parameter type
  mutate(variable.type = case_when(
    coeff %in% c("k.d", "p.d") ~ "1.antigenic.dist",
    coeff %in% c("k.v", "p.v") ~ "2.vax.eff",
    coeff %in% c("k.dv", "p.dv") ~ "3.AD.VE.interaction"),
    var.subtype.lag = paste0(variable.type, "_", subtype, "_lag", season.lag),
    coeff.type = ifelse(grepl("^k", coeff), "k", "p"),
    dist.coeff.type = paste0(coeff.type, "_", distance.type)
  ) %>%
  # select variables of interest
  select(var.subtype.lag, dist.coeff.type, value) %>%
  pivot_wider(names_from = dist.coeff.type,
              values_from = value) %>%
  # rearrange row order
  arrange(var.subtype.lag)


# insert empty rows in dataframe
library(berryFunctions)
res.all.sum.formatted = insertRows(res.all.sum, 
                           c(1,2,5,8,11,12,15,18,21,22,25,28), 
                           new = NA)
# remove first column
res.all.sum.formatted = res.all.sum.formatted %>%
  select(-var.subtype.lag)

write.csv(res.all.sum.formatted, "coefficient_table.csv",
          row.names = FALSE)
