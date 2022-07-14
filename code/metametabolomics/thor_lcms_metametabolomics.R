setwd("~/workspace/20220428-clean_thor_lcms/")
set.seed(8675309)

library(tidyverse)
library(cowplot)
library(UpSetR)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(ggVennDiagram)
library(scales)
#remotes::install_github("YixingTT/MTBrush/MTBrush") # uncomment this line to install
library(MTBrush)

## Read in data
cd.raw <- read.table('data/Compound_table_all_norm_and_raw_peak_area.csv', sep="\t", header=T, quote="\"")
cd.raw$Compounds.ID <- paste0('c', cd.raw$Compounds.ID)

## Format lookups for later
mid.lookup <- cd.raw %>% 
  select(Compounds.ID, Name, Tags, Formula, RT..min., Calc..MW)
colnames(mid.lookup)[6] <- 'Calc_MW'
colnames(mid.lookup)[5] <- 'RT_min'
known.matches <- mid.lookup %>% 
  filter(Name!='')

## Pull out areas
area.qcnorm <- cd.raw %>% 
  select(Compounds.ID, matches("^Norm..Area..2"))

## Double check that Compounds.ID is a unique identifier
nrow(
  area.qcnorm %>% 
  group_by(Compounds.ID) %>% 
  summarize(n=n()) %>% 
  filter(n>1)
) ## should return 0

## Transform to tall and break out by condition
area.tall <- area.qcnorm %>% 
  gather(key='replicate', value='peak.area', matches("^Norm..Area..2"))
area.tall$replicate <- sub("Norm..Area..20210525_(\\w+\\d+).raw.+", "\\1", area.tall$replicate)
area.tall$condition <- sub("\\d+$", "", area.tall$replicate)

## Save some maps for later
conditions <- unique(area.tall$condition)
rep2cond <- area.tall %>% 
  select(replicate, condition) %>% 
  distinct()

## Filter based on presence in QC samples
area.tall.qcfilt <- area.tall %>%
  filter(!is.na(peak.area))
area.qcfilt <- area.tall.qcfilt %>% 
  select(-condition) %>% 
  spread(replicate, peak.area)

## Norm to run
area.runsum <- area.tall.qcfilt %>% 
  group_by(replicate) %>% 
  summarize(tic = sum(peak.area))
area.tall.qcfilt <- area.tall.qcfilt %>%
  left_join(area.runsum, by='replicate') %>% 
  mutate(peak.area.norm = peak.area/tic)
area.qcfilt.norm <- area.tall.qcfilt %>% 
  select(-condition, -tic, -peak.area) %>% 
  spread(replicate, peak.area.norm)

## Lokisin subnetwork
lok.ann <- cd.raw %>% 
  filter(Tags=='A;B') %>% 
  select(Compounds.ID, Name)
lok <- area.tall.qcfilt %>%
  filter(Compounds.ID %in% lok.ann$Compounds.ID) %>% 
  left_join(lok.ann, by='Compounds.ID') %>% 
  filter(condition %in% c('K', 'BK', 'FK', 'BFK'))
lok$condition <- factor(lok$condition, levels=c('B', 'F', 'K', 'BF', 'BK', 'FK', 'BFK'))
lok.desc <- cd.raw %>% 
  filter(Compounds.ID %in% lok.ann$Compounds.ID)

ggplot(lok, aes(y=peak.area, x=condition, color=condition))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  facet_wrap(~Name, scales='free') +
  theme_bw()+
  theme(
    legend.position = 'top',
    legend.title = element_blank(),
    strip.background = element_blank(),
    axis.title.x = element_blank()
  )+
  scale_y_continuous(labels = function(x) format(x, scientific = T))+
  expand_limits(y=0)+
  scale_color_brewer(palette='Set2')+
  ylab("QC-normalized Peak Area")
ggsave("lok_subnet/lok.pa.pdf", dpi=1600, height=8, width=7)

## Dig in to the lokisin subnetwork
lok.tic <- lok %>% 
  group_by(Compounds.ID, Name) %>% 
  summarize(tic=sum(peak.area), tic_l10=log10(tic))
write.table(lok.tic, 'lok_subnet/lok_tic.tsv', sep="\t", quote = F, row.names = F, col.names = F)

## Save some suppl tables
write.table(area.qcfilt, 'lcms_abundance_matrix.tsv', row.names = F, sep="\t")
write.table(mid.lookup, 'lcms_annotations.tsv', row.names = F, sep="\t")

## t-stats for WT
lcms <- area.qcfilt %>%
  pivot_longer(-Compounds.ID, names_to = "condition") %>%
  mutate(
    replicate = as.numeric(str_extract(condition, "[0-9]+")),
    condition = str_remove(condition, "[0-9]+"),
    B = str_detect(condition, "B"),
    F = str_detect(condition, "F"),
    K = str_detect(condition, "K"),
    D = str_detect(condition, "D"),
    M = str_detect(condition, "M")
  )
split_code <- function(z) { ifelse(z, 1, -1) }
split_df <- lcms %>%
  mutate(across(B:M, split_code)) %>%
  split_dataset(Compounds.ID)
lm_func <- function(x) {
  lm(log(value) ~ B * F * K, data = x)
}
stats_list <- list()
stats_list[["K"]] <- split_df %>%
  map(~ filter(., D != 1, condition != "QC")) %>%
  fit_statistics(lm_func, Compounds.ID) %>%
  mutate(term = factor(term, levels = c("(Intercept)", "B", "F", "K", "B:F", "B:K", "F:K", "B:F:K")))
group_list <- unique(lcms$Compounds.ID)
all_wt_lcms <- stats_list[["K"]] %>%
  select(Compounds.ID, term, statistic) %>%
  pivot_wider(names_from = "term", values_from = "statistic")
write.table(all_wt_lcms, "lcms_tstat_wt.tsv", sep = "\t", row.names = F, quote = F)

awl.t <- all_wt_lcms %>% 
  gather(key='variable', value = 'tstatistic', -Compounds.ID) %>% 
  filter(variable != '(Intercept)') %>% 
  mutate(variable = factor(variable, levels=c('B', 'F', 'K', 'B:F', 'B:K', 'F:K', 'B:F:K')))
int3way <- awl.t %>% 
  filter(variable=='B:F:K') %>% 
  mutate(prank = ntile(tstatistic, 1000)/10,
         bfkflag = case_when(
           prank > 97.5 ~ 'Upper',
           prank < 2.5 ~ 'Lower',
           T ~ 'All'
         )
         )
awl.t.hi <- awl.t %>% 
  filter(Compounds.ID %in% (int3way %>% filter(bfkflag=='Upper'))$Compounds.ID)
awl.t.lo <- awl.t %>% 
  filter(Compounds.ID %in% (int3way %>% filter(bfkflag=='Lower'))$Compounds.ID)

single_dist <- function(x){
  b <- 30
  ss <- awl.t %>% filter(variable==x)
  hi <- awl.t.hi %>% filter(variable==x)
  lo <- awl.t.lo %>% filter(variable==x)
  return(
    ggplot(ss, aes(x=tstatistic))+
      geom_histogram(bins=b, color='white', fill='gray50')+
      geom_histogram(data=hi, bins=b, fill=34, alpha=.7)+
      geom_histogram(data=lo, bins=b, fill=143, alpha=.7)+
      theme_bw()+theme(
        plot.title = element_text(hjust = 0.5)
      )+
      ylab("# MFs")+xlab("t-statistic")+ggtitle(x)
  )
}

plot_grid(
  plot_grid(
    single_dist('B'), single_dist('F'), single_dist('K'),
    nrow=1, labels=c('A', 'B', 'C')
  ),
  plot_grid(
    single_dist('B:F'), single_dist('B:K'), single_dist('F:K'),
    nrow=1, labels=c('D', 'E', 'F')
  ),  
  single_dist('B:F:K'),
  nrow=3, labels=c('', '', 'G')
)

ggsave("lcms_tstat_wt.pdf", height=8, width=10, dpi=1600)


source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

not.lohi <- awl.t %>% 
  filter(!(Compounds.ID %in% c(awl.t.hi$Compounds.ID, awl.t.lo$Compounds.ID)))
ggplot(awl.t, aes(x=variable, y=tstatistic))+
  geom_hline(yintercept = 0, color='black', linetype='dashed')+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),
                   adjust = 2, trim = TRUE, fill='gray90') +
  geom_boxplot(outlier.shape = NA, position = position_nudge(x=.25, y=0),
               width = .1, coef=0)+
  geom_point(data=not.lohi, position = position_jitter(width = .1), 
             size = .5, alpha=.3, color='gray40') +
  geom_point(data=awl.t.hi, position = position_jitter(width = .1), 
             size = .5, alpha=.5, color=34) +
  geom_point(data=awl.t.lo, position = position_jitter(width = .1, ), 
             size = .5, alpha=.5, color=143) +
  guides(fill=FALSE, color = FALSE)+
  ylab('t-statistic') + xlab('Term') +
  theme_bw()+theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face='bold'),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+
  facet_wrap(~variable, scales="free", nrow=1)
ggsave("lcms_tstat_wt_rain.pdf", height=3, width=14, dpi=1600)

int <- c('B:F', 'B:K', 'F:K', 'B:F:K')
ggplot(awl.t %>% filter(variable %in% int), aes(x=variable, y=tstatistic))+
  geom_hline(yintercept = 0, color='black', linetype='dashed')+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),
                   adjust = 2, trim = TRUE, fill='gray90') +
  geom_boxplot(outlier.shape = NA, position = position_nudge(x=.25, y=0),
               width = .1, coef=0)+
  geom_point(data=not.lohi %>% filter(variable %in% int), position = position_jitter(width = .1), 
             size = 1, alpha=.3, color='gray40') +
  geom_point(data=awl.t.hi %>% filter(variable %in% int), position = position_jitter(width = .1), 
             size = 1, alpha=.5, color=34) +
  geom_point(data=awl.t.lo %>% filter(variable %in% int), position = position_jitter(width = .1, ), 
             size = 1, alpha=.5, color=143) +
  guides(fill=FALSE, color = FALSE)+
  ylab('t-statistic') + xlab('Term') +
  theme_bw()+theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face='bold'),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+
  facet_wrap(~variable, scales="free", nrow=1)
ggsave("lcms_tstat_wt_rain_intonly.pdf", height=3, width=8, dpi=1600)


## Create the table for D
lm_func_d <- function(x) {
  lm(value ~ B * F * D, data = x)
}
stats_list[["D"]] <- split_df %>%
  map(~ filter(., K != 1, condition != "QC")) %>%
  fit_statistics(lm_func_d, Compounds.ID) %>%
  mutate(term = factor(term, levels = c("(Intercept)", "B", "F", "D", "B:F", "B:D", "F:D", "B:F:D")))
all_mut_lcms <- stats_list[["D"]] %>%
  select(Compounds.ID, term, statistic) %>%
  pivot_wider(names_from = "term", values_from = "statistic")
write.table(all_mut_lcms, "lcms_tstat_mut.tsv", sep = "\t", row.names = F, quote = F)

mid.lookup %>% 
  group_by(mf.type.wt) %>% summarize(n=n(), p=100*n/4904)



## Look for pet
pet2 <- area.qcfilt %>%
  gather(key='replicate', value='area', -Compounds.ID) %>% 
  left_join(rep2cond, by='replicate') %>% 
  filter(condition %in% c('BFK', 'BFD')) %>% 
  group_by(Compounds.ID, condition) %>% 
  summarize(ma=mean(area)) %>% 
  spread(condition, ma) %>%
  mutate(k.ratio=BFD/BFK, l2r=log(k.ratio, base=2)) %>% 
  left_join(mid.lookup, by='Compounds.ID') %>% 
  filter(Tags %in% c('A;B'))
ggplot(pet2, aes(x=l2r))+
  geom_histogram()

## Define a baseline based on koreenceine signal in non-K conditions
kor.ann <- cd.raw %>% 
  filter(Tags=='C') %>% 
  select(Compounds.ID, Name)
kor <- area.tall.qcfilt %>%
  filter(Compounds.ID %in% kor.ann$Compounds.ID) %>% 
  left_join(kor.ann, by='Compounds.ID')
notkkor <- kor %>% 
  filter(!(grepl("K", condition) | condition == 'QC'))
pres.cut <- max(notkkor$peak.area) + (1 * sd(notkkor$peak.area) )
pres.cut.sds <- (pres.cut - mean(notkkor$peak.area) ) / sd(notkkor$peak.area)
c.a <- ggplot(kor, aes(y=peak.area, x=condition, color=Name))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_dodge(.9), alpha=.5) +
  geom_hline(yintercept = pres.cut, linetype='dashed', color='red')+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    legend.position = 'top'
  )+
  ylab("Peak Area")+xlab("Condition")+
  scale_y_log10()
cp.all <- area.tall.qcfilt %>% 
  filter(!(condition %in% c('M', 'QC')))
cp.ss <- cp.all %>% 
  filter(peak.area < 3e7)
c.b <- ggplot(cp.all, aes(x=peak.area))+
  geom_freqpoly(bins=500)+
  geom_vline(xintercept = pres.cut, color='red', linetype='dashed')+
  theme_bw()+
  ylab("# MFs")+xlab("Peak Area")
c.c <- ggplot(cp.ss, aes(x=peak.area))+
  geom_freqpoly(bins=500)+
  geom_vline(xintercept = pres.cut, color='red', linetype='dashed')+
  theme_bw()+
  ylab("# MFs")+xlab("Peak Area")
c.d <- ggplot(cp.all, aes(x=peak.area))+
  geom_freqpoly(bins=500)+
  geom_vline(xintercept = pres.cut, color='red', linetype='dashed')+
  theme_bw()+
  ylab("# MFs")+xlab("Peak Area")+
  scale_x_log10()
plot_grid(c.a, 
          plot_grid(c.b, c.c, c.d, 
                    labels = c('B', 'C', 'D'),
                    nrow=1, axis = 'tb', align='hv'),
          labels = c('A', ''),
          nrow = 2)
ggsave("presence_cutoffs.pdf", height=8, width=8, dpi=1600)

## Set present in ea replicate by cutoff
area.tall.qcfilt <- area.tall.qcfilt %>% 
  mutate(peak.present = case_when(
    peak.area >= pres.cut ~ 1,
    T ~ 0)
  )

## Summarize sample to sample TIC variation
ggplot(area.tall.qcfilt, aes(x=condition, y=tic))+
  geom_boxplot()+
  theme_classic()+
  ylab("Total ion current")+
  xlab("Condition")

## Check area distributions
ggplot(area.tall.qcfilt %>% filter(condition != 'QC'), aes(x=peak.area.norm, color=condition))+
  geom_freqpoly()+
  theme_bw()+
  scale_x_log10()+
  ylab("Count")+xlab("TIC Normalized Peak Area")

## Calc occurance for each MF per condition
a.occur<- area.tall.qcfilt %>%
  group_by(Compounds.ID, condition) %>% 
  summarize(n.reps=sum(peak.present))
ggplot(a.occur, aes(x=n.reps))+
  geom_histogram()
rep.cutoff <- 4
a.tall.qcfilt <- area.tall.qcfilt %>%
  left_join(a.occur, by=c('Compounds.ID', 'condition'))

## Get a list of the media MFs for later
media.MF <- a.occur %>% 
  filter(condition=='M' & n.reps > 0)

## Calc standard deviation across all replicates, remove those with 0 variation or part of QC
## (req'd for downstream PCA)
area.sd <- a.tall.qcfilt %>%
  group_by(Compounds.ID) %>% 
  summarize(sd=sd(peak.area)) %>% 
  filter(sd>0)
a.cut.pca <- a.tall.qcfilt %>%
  filter(condition != 'QC') %>% 
  filter(Compounds.ID %in% area.sd$Compounds.ID)

## Break out into prcomp formatted matricies
a.cut <- a.cut.pca %>% 
  select(Compounds.ID, replicate, peak.area) %>% 
  spread(Compounds.ID, peak.area) %>% 
  column_to_rownames('replicate')

## Compute pres/abs maps
u.pres <- a.occur %>% 
  mutate(cond.present = case_when(
    n.reps >= rep.cutoff ~ 1,
    T ~ 0
  )) %>% 
  select(Compounds.ID, condition, cond.present) %>% 
  spread(condition, cond.present)
u.pres.nomedia <- u.pres %>% 
  filter(M==0) %>% 
  select(-matches('M')) %>% 
  column_to_rownames('Compounds.ID')

## types of metabolites
only.cocul <- u.pres.nomedia %>% 
  filter(B==0 & F==0 & K==0 & D==0) %>% 
  rownames_to_column('Compounds.ID')
only.trip <- only.cocul %>% 
  filter(BFK==1 | BFD==1) %>% 
  filter(BF==0 & BK==0 & BD==0 & FK==0 & FD==0)
other.cc <- only.cocul %>% 
  filter(!(Compounds.ID %in% only.trip$Compounds.ID))
mid.lookup <- mid.lookup %>% 
  mutate(mf.type = case_when(
    Compounds.ID %in% other.cc$Compounds.ID ~ 'Only found in coculture',
    Compounds.ID %in% only.trip$Compounds.ID ~ 'Only found in three-member community',
    T ~ 'none'
  ), in.media = case_when(
    Compounds.ID %in% media.MF$Compounds.ID ~ 1,
    T ~ 0
  ))


## PCA with unit norm
pca.u <- a.cut %>% 
  prcomp(scale = T, center = T)
pca.u.var <- data.frame(
  PC = paste0('PC', 1:length(pca.u$sdev)),
  PC.n = 1:length(pca.u$sdev),
  var.explained = (pca.u$sdev)^2/sum((pca.u$sdev)^2)*100
) %>% mutate(cumvar=cumsum(var.explained))
pca.u.var$PC <- factor(pca.u.var$PC,
                       levels = paste0('PC', 1:length(pca.u$sdev))
)
pca.u.vmap <- pca.u.var %>% 
  mutate(ve.label=paste0(PC, ' (', round(var.explained, 1), '%)')) %>% 
  select(PC, ve.label)

## Skee plot and cumulative variance for top PCs
n.pc <- 15
u.v <- ggplot(pca.u.var %>% filter(PC.n<=n.pc), aes(x=PC, y=var.explained))+
  geom_line(linetype='dashed', group=1)+
  geom_point()+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle=90)
  )+
  ylab("Variance explained (%)")+
  xlab("Principal component")
u.c <- ggplot(pca.u.var %>% filter(PC.n<=n.pc), aes(x=PC, y=cumvar))+
  geom_line(linetype='dashed', group=1)+
  geom_point()+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle=90)
  )+
  ylab("Cumulative variance explained (%)")+
  xlab("Principal component")

## Break out PC coords
pca.u.coord <- as.data.frame(pca.u$x) %>% 
  rownames_to_column('replicate') %>%
  left_join(rep2cond, by='replicate')

## Break out PC loadings coords
pca.u.load <- as.data.frame(pca.u$rotation) %>% 
  rownames_to_column('Compounds.ID') %>% 
  left_join(mid.lookup, by='Compounds.ID')

## Plot PC1-2
u.12 <- ggplot(pca.u.coord, aes(x=PC1, y=PC2, color=condition))+
  #geom_point()+
  geom_text(aes(label=replicate))+
  theme_bw()+
  theme(
    legend.position = 'none'
  )+
  xlab(pca.u.vmap[pca.u.vmap$PC=='PC1','ve.label']) +
  ylab(pca.u.vmap[pca.u.vmap$PC=='PC2','ve.label'])

## Plot PC1-2 loadings
l.12 <- ggplot(pca.u.load %>% filter(!(Compounds.ID %in% media.MF$Compounds.ID)), aes(x=PC1, y=PC2))+
  geom_point(alpha=.5, color='gray50')+
  geom_point(data=pca.u.load %>% filter(Compounds.ID %in% c(only.cocul$Compounds.ID)), aes(color=mf.type), alpha=.5)+
  theme_bw()+theme(
    legend.title = element_blank(),
    legend.position = 'bottom'
  )+
  xlab(pca.u.vmap[pca.u.vmap$PC=='PC1','ve.label']) +
  ylab(pca.u.vmap[pca.u.vmap$PC=='PC2','ve.label'])
  
plot_grid(plot_grid(u.v, u.c,
                    labels = c('A', 'B'),
                    nrow=1,
                    axis='lrtb', align='hv'),
          plot_grid(u.12, l.12,
                    labels = c('C', 'D'),
                    nrow=1,
                    axis='lrtb', align='hv'),
          nrow=2
)
ggsave("unitscale_area_pc.pdf", height=10, width=10, dpi=1600)


## Only WT

## types of metabolites
wonly.cocul <- u.pres.nomedia %>% 
  filter(B==0 & F==0 & K==0) %>% 
  rownames_to_column('Compounds.ID')
wonly.trip <- wonly.cocul %>% 
  filter(BFK==1) %>% 
  filter(BF==0 & BK==0 & FK==0)
wother.cc <- wonly.cocul %>% 
  filter(!(Compounds.ID %in% wonly.trip$Compounds.ID))
mid.lookup <- mid.lookup %>% 
  mutate(mf.type.wt = case_when(
    Compounds.ID %in% wother.cc$Compounds.ID ~ 'Only found in coculture',
    Compounds.ID %in% wonly.trip$Compounds.ID ~ 'Only found in three-member community',
    T ~ 'none'
  ))

## PCA with unit norm
wpca.u <- a.cut %>% 
  rownames_to_column('rep') %>% 
  filter(!(grepl("D", rep))) %>% 
  column_to_rownames(var='rep') %>% 
  prcomp(scale = T, center = T)
wpca.u.var <- data.frame(
  PC = paste0('PC', 1:length(wpca.u$sdev)),
  PC.n = 1:length(wpca.u$sdev),
  var.explained = (wpca.u$sdev)^2/sum((wpca.u$sdev)^2)*100
) %>% mutate(cumvar=cumsum(var.explained))
wpca.u.var$PC <- factor(wpca.u.var$PC,
                       levels = paste0('PC', 1:length(wpca.u$sdev))
)
wpca.u.vmap <- wpca.u.var %>% 
  mutate(ve.label=paste0(PC, ' (', round(var.explained, 1), '%)')) %>% 
  select(PC, ve.label)

## Skee plot and cumulative variance for top PCs
n.pc <- 15
wu.v <- ggplot(wpca.u.var %>% filter(PC.n<=n.pc), aes(x=PC, y=var.explained))+
  geom_line(linetype='dashed', group=1)+
  geom_point()+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle=90)
  )+
  ylab("Variance explained (%)")+
  xlab("Principal component")
wu.c <- ggplot(wpca.u.var %>% filter(PC.n<=n.pc), aes(x=PC, y=cumvar))+
  geom_line(linetype='dashed', group=1)+
  geom_point()+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle=90)
  )+
  ylab("Cumulative variance explained (%)")+
  xlab("Principal component")

## Break out PC coords
wpca.u.coord <- as.data.frame(wpca.u$x) %>% 
  rownames_to_column('replicate') %>%
  left_join(rep2cond, by='replicate')

## Break out PC loadings coords
wpca.u.load <- as.data.frame(wpca.u$rotation) %>% 
  rownames_to_column('Compounds.ID') %>% 
  left_join(mid.lookup, by='Compounds.ID')

## Plot PC1-2
wu.12 <- ggplot(wpca.u.coord, aes(x=PC1, y=PC2, color=condition))+
  #geom_point()+
  geom_text(aes(label=replicate))+
  theme_bw()+
  theme(
    legend.position = 'none'
  )+
  xlab(wpca.u.vmap[wpca.u.vmap$PC=='PC1','ve.label']) +
  ylab(wpca.u.vmap[wpca.u.vmap$PC=='PC2','ve.label'])

## Plot PC1-2 loadings
wl.12 <- ggplot(wpca.u.load %>% filter(!(Compounds.ID %in% media.MF$Compounds.ID)), aes(x=PC1, y=PC2))+
  geom_point(alpha=.5, color='gray50')+
  geom_point(data=wpca.u.load %>% filter(Compounds.ID %in% c(wonly.cocul$Compounds.ID)), aes(color=mf.type.wt), alpha=.5)+
  theme_bw()+theme(
    legend.title = element_blank(),
    legend.position = 'bottom'
  )+
  xlab(wpca.u.vmap[wpca.u.vmap$PC=='PC1','ve.label']) +
  ylab(wpca.u.vmap[wpca.u.vmap$PC=='PC2','ve.label'])

plot_grid(plot_grid(wu.v, wu.c,
                    labels = c('A', 'B'),
                    nrow=1,
                    axis='lrtb', align='hv'),
          plot_grid(wu.12, wl.12,
                    labels = c('C', 'D'),
                    nrow=1,
                    axis='lrtb', align='hv'),
          nrow=2
)
ggsave("wt_unitscale_area_pc.pdf", height=10, width=10, dpi=1600)






## Correlation matrix
a.cut.scale <- as.data.frame(scale(a.cut, scale=T, center=T))
a.cor <- cor(
  a.cut.scale %>% t(),
  method = 'spearman'
)

## Define annotations for the downstream plotting
annstrip <- rep2cond %>% 
  mutate(
    B = case_when(
      condition %in% c('B', 'BF', 'BK', 'BFK', 'BD', 'BFD') ~ 'Y',
      T ~ 'N'),
    F = case_when(
      condition %in% c('F', 'BF', 'FK', 'FD', 'BFD', 'BFK') ~ 'Y',
      T ~ 'N'),
    K = case_when(
      condition %in% c('K', 'BK', 'FK', 'BFK') ~ 'Y',
      T ~ 'N'),
    D = case_when(
      condition %in% c('D', 'BD', 'FD', 'BFD') ~ 'Y',
      T ~ 'N')
  ) %>% 
  select(replicate, B, F, K, D) %>% 
  column_to_rownames('replicate')
anncol <- list(
  B = c(Y = '#6ea8ca', N = '#ffffff'),
  F = c(Y = '#ca9c6e', N = '#ffffff'),
  K = c(Y = '#ca6e6e', N = '#ffffff'),
  D = c(Y = '#999999', N = '#ffffff')
) 

##Correlation heatmap
# pheatmap(
#   mat=a.cor,
#   border_color = NA,
#   color = colorRampPalette(
#     brewer.pal(n=7,name="RdBu")
#   )(255),
#   breaks = seq(-1,1,2/255),
#   fontsize = 5,
#   treeheight_row = 25,
#   treeheight_col = 25,
#   cutree_rows = 2,
#   cutree_cols = 7,
#   annotation_row = annstrip,
#   annotation_col = annstrip,
#   annotation_colors = anncol,
#   annotation_legend = F,
#   filename = 'lcms_correlation.pdf',
#   height=5,
#   width=5
# )

########################################
########################################
########################################

upset(u.pres.nomedia,
      sets = rev(c('K', 'D')),
      nintersects = NA,
      keep.order = T
)

upset(u.pres.nomedia,
      sets = rev(c('BK', 'BD')),
      nintersects = NA,
      keep.order = T
)

upset(u.pres.nomedia,
      sets = rev(c('FK', 'FD')),
      nintersects = NA,
      keep.order = T
)

upset(u.pres.nomedia,
      sets = rev(c('BFK', 'BFD')),
      nintersects = NA,
      keep.order = T
)




vk <- ggVennDiagram(
  list(
    (u.pres.nomedia %>%
      rownames_to_column('Compounds.ID') %>% 
      gather(key="condition", value="present", -Compounds.ID) %>% 
      filter(condition=='K') %>% 
      filter(present==1))$Compounds.ID,
    (u.pres.nomedia %>%
       rownames_to_column('Compounds.ID') %>% 
       gather(key="condition", value="present", -Compounds.ID) %>% 
       filter(condition=='D') %>% 
       filter(present==1))$Compounds.ID
  ),label_percent_digit=1
)
vbk <- ggVennDiagram(
  list(
    (u.pres.nomedia %>%
       rownames_to_column('Compounds.ID') %>% 
       gather(key="condition", value="present", -Compounds.ID) %>% 
       filter(condition=='BK') %>% 
       filter(present==1))$Compounds.ID,
    (u.pres.nomedia %>%
       rownames_to_column('Compounds.ID') %>% 
       gather(key="condition", value="present", -Compounds.ID) %>% 
       filter(condition=='BD') %>% 
       filter(present==1))$Compounds.ID
  ),label_percent_digit=1
)
vfk <- ggVennDiagram(
  list(
    (u.pres.nomedia %>%
       rownames_to_column('Compounds.ID') %>% 
       gather(key="condition", value="present", -Compounds.ID) %>% 
       filter(condition=='FK') %>% 
       filter(present==1))$Compounds.ID,
    (u.pres.nomedia %>%
       rownames_to_column('Compounds.ID') %>% 
       gather(key="condition", value="present", -Compounds.ID) %>% 
       filter(condition=='FD') %>% 
       filter(present==1))$Compounds.ID
  ),label_percent_digit=1
)
vbfk <- ggVennDiagram(
  list(
    (u.pres.nomedia %>%
       rownames_to_column('Compounds.ID') %>% 
       gather(key="condition", value="present", -Compounds.ID) %>% 
       filter(condition=='BFK') %>% 
       filter(present==1))$Compounds.ID,
    (u.pres.nomedia %>%
       rownames_to_column('Compounds.ID') %>% 
       gather(key="condition", value="present", -Compounds.ID) %>% 
       filter(condition=='BFD') %>% 
       filter(present==1))$Compounds.ID
  ),label_percent_digit=1
)

plot_grid(vk, vbk, vfk, vbfk,
          nrow=1)

vk
vbk
vfk
vbfk

####################
####################
####################
####################


## THOR WT
upset(u.pres.nomedia,
  sets = rev(c('B', 'F', 'K', 'BF', 'BK', 'FK', 'BFK')),
  nintersects = NA, 
  keep.order = T
)
## All of it
upset(u.pres.nomedia,
      sets = rev(c('B', 'F', 'K', 'BF', 'BK', 'FK', 'BFK', 'D', 'BD', 'FD', 'BFD')),
      nintersects = NA, 
      keep.order = T
)

## Total mfs
tots <- u.pres.nomedia %>%
  summarize(across(everything(), sum)) %>% 
  t() %>% 
  as.data.frame() %>% 
  rename(total.MFs = V1)

## WT-only PCA with unit norm
wt <- c('B', 'F', 'K', 'BF', 'BK', 'FK', 'BFK', 'M')
wt.cut.tall <- a.cut.tall %>% 
  filter(condition %in% wt) 

## Calc standard deviation across all replicates, remove those with 0 variation
## (req'd for downstream PCA)
wt.sd <- wt.cut.tall %>% 
  group_by(Compounds.ID) %>% 
  summarize(sd=sd(peak.area.norm.cut)) %>% 
  filter(sd>0)
wt.cut.pca <- wt.cut.tall %>% 
  filter(Compounds.ID %in% wt.sd$Compounds.ID)

## Break out into prcomp formatted matricies
wt.cut <- wt.cut.pca %>% 
  select(Compounds.ID, replicate, peak.area.norm.cut) %>% 
  spread(Compounds.ID, peak.area.norm.cut) %>% 
  column_to_rownames('replicate')

pca.wt <- wt.cut %>% 
  prcomp(scale = T, center = T)
pca.wt.var <- data.frame(
  PC = paste0('PC', 1:length(pca.wt$sdev)),
  PC.n = 1:length(pca.wt$sdev),
  var.explained = (pca.wt$sdev)^2/sum((pca.wt$sdev)^2)*100
) %>% mutate(cumvar=cumsum(var.explained))
pca.wt.var$PC <- factor(pca.wt.var$PC,
                       levels = paste0('PC', 1:length(pca.wt$sdev))
)
pca.wt.vmap <- pca.wt.var %>% 
  mutate(ve.label=paste0(PC, ' (', round(var.explained, 1), '%)')) %>% 
  select(PC, ve.label)

## Skee plot and cumulative variance for top PCs
n.pc <- 10
wt.v <- ggplot(pca.wt.var %>% filter(PC.n<=n.pc), aes(x=PC, y=var.explained))+
  geom_line(linetype='dashed', group=1)+
  geom_point()+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle=90)
  )+
  ylab("Variance explained (%)")+
  xlab("Principal component")
wt.c <- ggplot(pca.wt.var %>% filter(PC.n<=n.pc), aes(x=PC, y=cumvar))+
  geom_line(linetype='dashed', group=1)+
  geom_point()+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle=90)
  )+
  ylab("Cumulative variance explained (%)")+
  xlab("Principal component")

## Break out PC coords
pca.wt.coord <- as.data.frame(pca.wt$x) %>% 
  rownames_to_column('replicate') %>%
  left_join(rep2cond, by='replicate')
pca.wt.coord$condition <- factor(pca.wt.coord$condition, levels=wt)

## Break out PC loadings coords
pca.wt.load <- as.data.frame(pca.wt$rotation) %>% 
  rownames_to_column('Compounds.ID')

## Plot PC1-2
wt.12 <- ggplot(pca.wt.coord, aes(x=PC1, y=PC2, fill=condition))+
  geom_point(shape=23, color='black', alpha=.8, size=3)+
  #geom_text(aes(label=replicate))+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    legend.position = c(.75,.8),
    legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')
  )+
  xlab(pca.wt.vmap[pca.wt.vmap$PC=='PC1','ve.label'])+
  ylab(pca.wt.vmap[pca.wt.vmap$PC=='PC2','ve.label'])+
  scale_fill_brewer(palette='Set2')+
  guides(fill=guide_legend(ncol=2))

## Plot PC1-2 loadings
k.lab <- pca.wt.load %>% 
  filter(Compounds.ID %in% knowns$Compounds.ID) %>% 
  left_join(knowns, by='Compounds.ID')
wtl.12 <- ggplot(pca.wt.load %>% filter(!(Compounds.ID %in% media.MF$Compounds.ID)), aes(x=PC1, y=PC2))+
  geom_point(alpha=.1)+
  geom_point(data=k.lab, color='red')+
  geom_label_repel(data=k.lab, aes(label = mf.name),
                   segment.color = 'red')+
  theme_bw()+
  xlab(pca.wt.vmap[pca.wt.vmap$PC=='PC1','ve.label']) +
  ylab(pca.wt.vmap[pca.wt.vmap$PC=='PC2','ve.label'])

plot_grid(wt.v, wt.c, wt.12, wtl.12,
          labels=LETTERS[1:4],
          nrow=2
)
ggsave("wt_normunitscale_area_pc.pdf", height=8, width=8, dpi=1600)

## Zoom in of the Pk samples
ggplot(pca.wt.coord %>% filter(PC1>0), aes(x=PC1, y=PC2, color=condition))+
  #geom_point()+
  geom_text(aes(label=replicate))+
  theme_bw()+
  theme(
    legend.position = 'none'
  )+
  xlab(pca.wt.vmap[pca.u.vmap$PC=='PC1','ve.label']) +
  ylab(pca.wt.vmap[pca.u.vmap$PC=='PC2','ve.label'])

## Kor loadings dive
## kor absolute
kor.l <- pca.wt.load %>% 
  filter(Compounds.ID %in% knowns$Compounds.ID) %>% 
  gather("PC", "load", -Compounds.ID) %>% 
  mutate(abs=abs(load)) %>% arrange(-abs)
kor.l.sum <- kor.l %>% 
  filter(Compounds.ID %in% c('c61', 'c96')) %>% 
  select(-abs) %>% 
  spread(key=Compounds.ID, value=load) %>% 
  mutate(dif=abs(c61-c96))

p56c <- ggplot(pca.wt.coord, aes(x=PC2, y=PC5, fill=condition))+
  geom_point(shape=23, color='black', alpha=.8, size=3)+
  #geom_text(aes(label=replicate))+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    legend.position = c(.8,.2),
    legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')
  )+
  xlab(pca.wt.vmap[pca.wt.vmap$PC=='PC2','ve.label'])+
  ylab(pca.wt.vmap[pca.wt.vmap$PC=='PC5','ve.label'])+
  scale_fill_brewer(palette='Set2')+
  guides(fill=guide_legend(ncol=2))
p56l <- ggplot(pca.wt.load %>% filter(!(Compounds.ID %in% media.MF$Compounds.ID)), aes(x=PC2, y=PC5))+
  geom_point(alpha=.1)+
  geom_point(data=k.lab, color='red')+
  geom_label_repel(data=k.lab, aes(label = mf.name),
                   segment.color = 'red')+
  theme_bw()+
  xlab(pca.wt.vmap[pca.wt.vmap$PC=='PC2','ve.label']) +
  ylab(pca.wt.vmap[pca.wt.vmap$PC=='PC5','ve.label'])

k <- c('K', 'BK', 'FK', 'BFK')
kor.abu <- area.norm.tall %>% 
  filter(Compounds.ID %in% knowns$Compounds.ID) %>% 
  filter(condition %in% k)
kor.abu$Compounds.ID <- as.factor(kor.abu$Compounds.ID)
kor.abu <- kor.abu %>% 
  left_join(knowns, by='Compounds.ID')

mA <- aov(peak.area.norm ~ condition,
          data = kor.abu %>% filter(mf.name=='korA'))
tA <- as.data.frame(TukeyHSD(mA, conf.level=.95)$condition) %>% 
  select(`p adj`) %>% 
  rownames_to_column('comparison')
mB <- aov(peak.area.norm ~ condition,
          data = kor.abu %>% filter(mf.name=='korB'))
tB <- as.data.frame(TukeyHSD(mB, conf.level=.95)$condition) %>% 
  select(`p adj`) %>% 
  rownames_to_column('comparison')
mC <- aov(peak.area.norm ~ condition,
          data = kor.abu %>% filter(mf.name=='korC'))
tC <- as.data.frame(TukeyHSD(mC, conf.level=.95)$condition) %>% 
  select(`p adj`) %>% 
  rownames_to_column('comparison')
colnames(tA)[2] <- 'korA.p.adj'
colnames(tB)[2] <- 'korB.p.adj'
colnames(tC)[2] <- 'korC.p.adj'
t.all <- tA %>% 
  left_join(tB, by='comparison') %>% 
  left_join(tC, by='comparison')
write.table(t.all, 'kor_tukey.tsv', sep="\t", row.names = F)

kor.sg <- data.frame(
  condition = c(k, k, k),
  mf.name = c('korA', 'korA', 'korA', 'korA', 
              'korB', 'korB', 'korB', 'korB', 
              'korC', 'korC', 'korC', 'korC'),
  sig.grp = c('a', 'b', 'ab', 'ab',
              'a', 'b', 'a', 'b',
              'a', 'a', 'a', 'a')
)

kor.abu$condition <- factor(kor.abu$condition, levels=k)

k.abu.p <- ggplot(kor.abu, aes(x=condition, y=peak.area.norm, color=mf.name))+
  geom_boxplot(alpha=.8, outlier.shape = NA)+
  geom_jitter(alpha=.8, width = 0.25, size=2)+
  geom_label(data=kor.sg, aes(y=0.014, label=sig.grp), fontface='bold', label.size = NA, color='black')+
  facet_wrap(~mf.name)+
  theme_bw()+theme(
    legend.position = 'none',
    axis.title.x = element_blank()
  )+
  ylab("TIC Normalized Peak Area")+
  scale_color_brewer(palette='Set1')+
  scale_y_continuous(limits=c(0,0.015))

plot_grid(
  plot_grid(p56c, p56l,
            nrow=1,
            labels=LETTERS[1:2]),
  k.abu.p,
  nrow=2,
  labels=c('', LETTERS[3]),
  rel_heights = c(1,0.75)
)
ggsave("wt_kor.pdf", dpi=1600, height=8, width=9)

## Calc uniqueness
only.c <- u.pres.nomedia %>%
  select(B, F, K, BF, BK, FK, BFK) %>% 
  mutate(tot=rowSums(.)) %>% 
  filter(tot==1) %>% 
  select(-tot) %>% 
  gather(key='condition',value='present') %>%
  group_by(condition) %>% 
  summarize(only=sum(present)) %>% 
  left_join(tots %>% rownames_to_column("condition"), by='condition') %>% 
  mutate(pct.only=100*only/total.MFs) %>% 
  select(condition, only, pct.only)
only2.c <- u.pres.nomedia %>%
  select(B, F, K, BF, BK, FK, BFK) %>% 
  mutate(tot=rowSums(.)) %>% 
  filter(tot==2) %>% 
  select(-tot) %>% 
  gather(key='condition',value='present') %>%
  group_by(condition) %>% 
  summarize(only2=sum(present)) %>% 
  left_join(tots %>% rownames_to_column("condition"), by='condition') %>% 
  mutate(pct.only2=100*only2/total.MFs) %>% 
  select(condition, only2, pct.only2)
  
l <- data.frame(
  type = c('only', 'only2', 'not'),
  lab = c('Unique to this condition', 'Present in exactly 2 conditions', 'Present in >2 conditions')
)
ov <- only.c %>% 
  left_join(only2.c, by='condition') %>% 
  left_join(tots %>% rownames_to_column("condition"), by='condition') %>% 
  mutate(not = total.MFs-only-only2, pct.not = 100*not/total.MFs) %>% 
  select(condition, only, only2, not) %>% 
  gather(key='type', value='mfs', -condition) %>% 
  left_join(l, by='type')

ov$condition <- factor(ov$condition, levels = rev(c('B', 'F', 'K', 'BF', 'BK', 'FK', 'BFK')))
op <- ggplot(ov, aes(x=condition, y=mfs, fill=lab))+
  geom_bar(stat='identity', size=.5, color='black', alpha=.6)+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    legend.position = 'bottom',
    axis.text.y = element_text(face='bold'),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face='bold')
  )+
  coord_flip() +
  ylab("Molecular features")+
  scale_fill_brewer(palette='Set2', guide=guide_legend(reverse=T))

b.u.f <- u.pres.nomedia %>% 
  rownames_to_column("Compounds.ID") %>% 
  filter(B>0 | F>0)
b.u.k <- u.pres.nomedia %>% 
  rownames_to_column("Compounds.ID") %>% 
  filter(B>0 | K>0) 
f.u.k <- u.pres.nomedia %>% 
  rownames_to_column("Compounds.ID") %>% 
  filter(F>0 | K>0) 
u.sing <- u.pres.nomedia %>% 
  rownames_to_column("Compounds.ID") %>% 
  filter(B>0 | F>0 | K>0)
u.pairs <- u.pres.nomedia %>% 
  rownames_to_column("Compounds.ID") %>% 
  filter(BF>0 | BK>0 | FK>0)

pair <- u.pres.nomedia %>% 
  rownames_to_column("Compounds.ID") %>% 
  filter(BF>0 | BK>0 | FK>0) %>% 
  select(Compounds.ID, BF, BK, FK) %>% 
  gather(key='condition', value='present', -Compounds.ID) %>% 
  filter(present>0)

triple <- u.pres.nomedia %>%
  rownames_to_column("Compounds.ID") %>% 
  filter(BFK>0) %>% 
  select(Compounds.ID) %>% 
  mutate(condition='BFKsing', present=1)
triple <- triple %>% 
  bind_rows(
    u.pres.nomedia %>%
      rownames_to_column("Compounds.ID") %>% 
      filter(BFK>0) %>% 
      select(Compounds.ID) %>% 
      mutate(condition='BFKpair', present=1)
  )

pair <- pair %>% bind_rows(triple)

pair$union <- NA
pair[pair$condition == 'BF' & pair$Compounds.ID %in% b.u.f$Compounds.ID,]$union <- 'Y'
pair[pair$condition == 'BF' & !(pair$Compounds.ID %in% b.u.f$Compounds.ID),]$union <- 'N'
pair[pair$condition == 'BK' & pair$Compounds.ID %in% b.u.k$Compounds.ID,]$union <- 'Y'
pair[pair$condition == 'BK' & !(pair$Compounds.ID %in% b.u.k$Compounds.ID),]$union <- 'N'
pair[pair$condition == 'FK' & pair$Compounds.ID %in% f.u.k$Compounds.ID,]$union <- 'Y'
pair[pair$condition == 'FK' & !(pair$Compounds.ID %in% f.u.k$Compounds.ID),]$union <- 'N'
pair[pair$condition == 'BFKsing' & pair$Compounds.ID %in% u.sing$Compounds.ID,]$union <- 'Y'
pair[pair$condition == 'BFKsing' & !(pair$Compounds.ID %in% u.sing$Compounds.ID),]$union <- 'N'
pair[pair$condition == 'BFKpair' & pair$Compounds.ID %in% u.pairs$Compounds.ID,]$union <- 'Y'
pair[pair$condition == 'BFKpair' & !(pair$Compounds.ID %in% u.pairs$Compounds.ID),]$union <- 'N'
u <- data.frame(condition=c('BF', 'BK', 'FK', 'BFKsing', 'BFKpair'),
                union= c('U', 'U', 'U', 'U', 'U'),
                n=c(
                  nrow(b.u.f %>% filter(!(Compounds.ID %in% (pair %>% filter(condition=='BF'))$Compounds.ID))),
                  nrow(b.u.k %>% filter(!(Compounds.ID %in% (pair %>% filter(condition=='BK'))$Compounds.ID))),
                  nrow(f.u.k %>% filter(!(Compounds.ID %in% (pair %>% filter(condition=='FK'))$Compounds.ID))),
                  nrow(u.sing %>% filter(!(Compounds.ID %in% (pair %>% filter(condition=='BFKsing'))$Compounds.ID))),
                  nrow(u.pairs %>% filter(!(Compounds.ID %in% (pair %>% filter(condition=='BFKpair'))$Compounds.ID)))
                )
  )

pair.c <- pair %>% 
  group_by(condition, union) %>% 
  summarize(n=n()) %>% 
  bind_rows(u)

pn <- c('Individual components only', 'Intersect', 'Joint condition only')
plab <- data.frame(
  union=c('U', 'Y', 'N'),
  lab= pn
)
pair.c <- pair.c %>% 
  left_join(plab, by='union')
pair.c$condition <- factor(pair.c$condition, levels=rev(c('BF', 'BK', 'FK', 'BFKsing', 'BFKpair')))
pair.c$lab <- factor(pair.c$lab, levels=rev(pn))

cl <- c("(B \u222a F) vs. BF", "(B \u222a K) vs. BK", "(F \u222a K) vs. FK", "(B \u222a F \u222a K)\nvs. BFK", "(BF \u222a BK \u222a FK)\nvs. BFK")
clab <- data.frame(
  condition=c('BF', 'BK', 'FK', 'BFKsing', 'BFKpair'),
  clab=cl
)
pair.c <- pair.c %>% 
  left_join(clab, by='condition')
pair.c$clab <- factor(pair.c$clab, levels=rev(cl))

pp <- ggplot(pair.c, aes(x=clab, y=n, fill=lab))+
  geom_bar(stat='identity', size=.5, color='black', alpha=.6)+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    legend.position = 'bottom',
    axis.text.y = element_text(face='bold'),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face='bold')
  )+
  coord_flip() +
  ylab("Molecular features")+
  scale_fill_brewer(palette='Set1', guide=guide_legend(reverse=T))

b.v <- list(
  B = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(B>0))$Compounds.ID,
  BF = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BF>0))$Compounds.ID,
  BK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BK>0))$Compounds.ID,
  BFK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BFK>0))$Compounds.ID
)
f.v <- list(
  F = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(F>0))$Compounds.ID,
  BF = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BF>0))$Compounds.ID,
  FK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(FK>0))$Compounds.ID,
  BFK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BFK>0))$Compounds.ID
)
k.v <- list(
  K = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(K>0))$Compounds.ID,
  BK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BK>0))$Compounds.ID,
  FK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(FK>0))$Compounds.ID,
  BFK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BFK>0))$Compounds.ID
)
t.v <- list(
  BF = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BF>0))$Compounds.ID,
  BK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BK>0))$Compounds.ID,
  FK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(FK>0))$Compounds.ID,
  BFK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BFK>0))$Compounds.ID
)

{
lsz <- unit(3.5, 'mm')
lpo <- c(.95,.21)
csz <- 2.5
bvdata <- process_data(Venn(b.v))
bp2 <- ggplot() +
  geom_sf(aes(fill=count), data = venn_region(bvdata)) +
  geom_sf(size = 1, color = "#ffffff", data = venn_setedge(bvdata), show.legend = F) +
  geom_sf_text(aes(label = name), fontface='bold', data = venn_setlabel(bvdata)) +
  geom_sf_label(aes(label=count), size=csz, alpha=.9, data = venn_region(bvdata)) +
  theme_void()+theme(
    legend.title=element_text(face='bold'),
    legend.key.size = lsz,
    legend.position = lpo
  )+
  guides(fill=guide_colorbar(title="MFs"))
fvdata <- process_data(Venn(f.v))
fp2 <- ggplot() +
  geom_sf(aes(fill=count), data = venn_region(fvdata)) +
  geom_sf(size = 1, color = "#ffffff", data = venn_setedge(fvdata), show.legend = F) +
  geom_sf_text(aes(label = name), fontface='bold', data = venn_setlabel(fvdata)) +
  geom_sf_label(aes(label=count), size=csz, alpha=.9, data = venn_region(fvdata)) +
  theme_void()+theme(
    legend.title=element_text(face='bold'),
    legend.key.size = lsz,
    legend.position = lpo
  )+
  guides(fill=guide_colorbar(title="MFs"))
kvdata <- process_data(Venn(k.v))
kp2 <- ggplot() +
  geom_sf(aes(fill=count), data = venn_region(kvdata)) +
  geom_sf(size = 1, color = "#ffffff", data = venn_setedge(kvdata), show.legend = F) +
  geom_sf_text(aes(label = name), fontface='bold', data = venn_setlabel(kvdata)) +
  geom_sf_label(aes(label=count), size=csz, alpha=.9, data = venn_region(kvdata)) +
  theme_void()+theme(
    legend.title=element_text(face='bold'),
    legend.key.size = lsz,
    legend.position = lpo
  )+
  guides(fill=guide_colorbar(title="MFs"))
tvdata <- process_data(Venn(t.v))
tp2 <- ggplot() +
  geom_sf(aes(fill=count), data = venn_region(tvdata)) +
  geom_sf(size = 1, color = "#ffffff", data = venn_setedge(tvdata), show.legend = F) +
  geom_sf_text(aes(label = name), fontface='bold', data = venn_setlabel(tvdata)) +
  geom_sf_label(aes(label=count), size=csz, alpha=.9, data = venn_region(tvdata)) +
  theme_void()+theme(
    legend.title=element_text(face='bold'),
    legend.key.size = lsz,
    legend.position = lpo
  )+
  guides(fill=guide_colorbar(title="MFs"))

plot_grid(op,
          plot_grid(bp2, fp2,
                    kp2, tp2,
                    nrow=2),
          pp,
          nrow=3,
          rel_heights = c(0.5, 1, .5))
ggsave("venn4s_heat.pdf", dpi=1600, height=12, width=9, device=cairo_pdf)
}


## 2-member D comparisons
d.v <- list(
  K = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(K>0))$Compounds.ID,
  D = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(D>0))$Compounds.ID
)
db.v <- list(
  BK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BK>0))$Compounds.ID,
  BD = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BD>0))$Compounds.ID
)
df.v <- list(
  FK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(FK>0))$Compounds.ID,
  FD = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(FD>0))$Compounds.ID
)
dt.v <- list(
  BFK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BFK>0))$Compounds.ID,
  BFD = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BFD>0))$Compounds.ID
)

lsz2 <- unit(5, 'mm')
lpo2 <- c(1, .21)
csz2 <- 7.5
ssz2 <- 8
ddata <- process_data(Venn(d.v))
dv2 <- ggplot() +
  geom_sf(aes(fill=count), data = venn_region(ddata)) +
  geom_sf(size = 1, color = "#ffffff", data = venn_setedge(ddata), show.legend = F) +
  geom_sf_text(aes(label = name), size=ssz2, fontface='bold', data = venn_setlabel(ddata)) +
  geom_sf_label(aes(label=count), size=csz2, alpha=.9, data = venn_region(ddata)) +
  theme_void()+theme(
    legend.title=element_text(face='bold'),
    legend.key.size = lsz2,
    legend.position = lpo2
  )+
  guides(fill=guide_colorbar(title="MFs"))
dbdata <- process_data(Venn(db.v))
dbv2 <- ggplot() +
  geom_sf(aes(fill=count), data = venn_region(dbdata)) +
  geom_sf(size = 1, color = "#ffffff", data = venn_setedge(dbdata), show.legend = F) +
  geom_sf_text(aes(label = name), size=ssz2, fontface='bold', data = venn_setlabel(dbdata)) +
  geom_sf_label(aes(label=count), size=csz2, alpha=.9, data = venn_region(dbdata)) +
  theme_void()+theme(
    legend.title=element_text(face='bold'),
    legend.key.size = lsz2,
    legend.position = lpo2
  )+
  guides(fill=guide_colorbar(title="MFs"))
dfdata <- process_data(Venn(df.v))
dfv2 <- ggplot() +
  geom_sf(aes(fill=count), data = venn_region(dfdata)) +
  geom_sf(size = 1, color = "#ffffff", data = venn_setedge(dfdata), show.legend = F) +
  geom_sf_text(aes(label = name), size=ssz2, fontface='bold', data = venn_setlabel(dfdata)) +
  geom_sf_label(aes(label=count), size=csz2, alpha=.9, data = venn_region(dfdata)) +
  theme_void()+theme(
    legend.title=element_text(face='bold'),
    legend.key.size = lsz2,
    legend.position = lpo2
  )+
  guides(fill=guide_colorbar(title="MFs"))
dtdata <- process_data(Venn(dt.v))
dtv2 <- ggplot() +
  geom_sf(aes(fill=count), data = venn_region(dtdata)) +
  geom_sf(size = 1, color = "#ffffff", data = venn_setedge(dtdata), show.legend = F) +
  geom_sf_text(aes(label = name), size=ssz2, fontface='bold', data = venn_setlabel(dtdata)) +
  geom_sf_label(aes(label=count), size=csz2, alpha=.9, data = venn_region(dtdata)) +
  theme_void()+theme(
    legend.title=element_text(face='bold'),
    legend.key.size = lsz2,
    legend.position = lpo2
  )+
  guides(fill=guide_colorbar(title="MFs"))

wfv <- list(
  K = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(K>0))$Compounds.ID,
  FK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(FK>0))$Compounds.ID,
  FD = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(FD>0))$Compounds.ID,
  D = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(D>0))$Compounds.ID
)
wfdata <- process_data(Venn(wfv))
wfp <- ggplot() +
  geom_sf(aes(fill=count), data = venn_region(wfdata)) +
  geom_sf(size = 1, color = "#ffffff", data = venn_setedge(wfdata), show.legend = F) +
  geom_sf_text(aes(label = name), size=ssz2, fontface='bold', data = venn_setlabel(wfdata)) +
  geom_sf_label(aes(label=count), size=csz2, alpha=.9, data = venn_region(wfdata)) +
  theme_void()+theme(
    legend.title=element_text(face='bold'),
    legend.key.size = lsz2,
    legend.position = lpo2
  )+
  guides(fill=guide_colorbar(title="MFs"))

wbv <- list(
  K = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(K>0))$Compounds.ID,
  BK = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BK>0))$Compounds.ID,
  BD = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(BD>0))$Compounds.ID,
  D = (u.pres.nomedia %>% rownames_to_column('Compounds.ID') %>% filter(D>0))$Compounds.ID
)
wbdata <- process_data(Venn(wbv))
wbp <- ggplot() +
  geom_sf(aes(fill=count), data = venn_region(wbdata)) +
  geom_sf(size = 1, color = "#ffffff", data = venn_setedge(wbdata), show.legend = F) +
  geom_sf_text(aes(label = name), size=ssz2, fontface='bold', data = venn_setlabel(wbdata)) +
  geom_sf_label(aes(label=count), size=csz2, alpha=.9, data = venn_region(wbdata)) +
  theme_void()+theme(
    legend.title=element_text(face='bold'),
    legend.key.size = lsz2,
    legend.position = lpo2
  )+
  guides(fill=guide_colorbar(title="MFs"))

plot_grid(dv2, dbv2, dfv2, dtv2, wbp, wfp,
          nrow=3
)
ggsave("kecmutv.pdf", dpi=1600, height=12, width=14)


## Explore Pseudo vs kec del conditions PCA
## K-only PCA with unit norm
kkeep <- c('K', 'BK', 'FK', 'BFK', 'D', 'BD', 'FD', 'BFD')
kk.cut.tall <- a.cut.tall %>% 
  filter(condition %in% kkeep) 

## Calc standard deviation across all replicates, remove those with 0 variation
## (req'd for downstream PCA)
kk.sd <- kk.cut.tall %>% 
  group_by(Compounds.ID) %>% 
  summarize(sd=sd(peak.area.norm.cut)) %>% 
  filter(sd>0)
kk.cut.pca <- kk.cut.tall %>% 
  filter(Compounds.ID %in% kk.sd$Compounds.ID)

## Break out into prcomp formatted matricies
kk.cut <- kk.cut.pca %>% 
  select(Compounds.ID, replicate, peak.area.norm.cut) %>% 
  spread(Compounds.ID, peak.area.norm.cut) %>% 
  column_to_rownames('replicate')

pca.kk <- kk.cut %>% 
  prcomp(scale = T, center = T)
pca.kk.var <- data.frame(
  PC = paste0('PC', 1:length(pca.kk$sdev)),
  PC.n = 1:length(pca.kk$sdev),
  var.explained = (pca.kk$sdev)^2/sum((pca.kk$sdev)^2)*100
) %>% mutate(cumvar=cumsum(var.explained))
pca.kk.var$PC <- factor(pca.kk.var$PC,
                        levels = paste0('PC', 1:length(pca.kk$sdev))
)
pca.kk.vmap <- pca.kk.var %>% 
  mutate(ve.label=paste0(PC, ' (', round(var.explained, 1), '%)')) %>% 
  select(PC, ve.label)

## Skee plot and cumulative variance for top PCs
n.pc <- 10
kk.v <- ggplot(pca.kk.var %>% filter(PC.n<=n.pc), aes(x=PC, y=var.explained))+
  geom_line(linetype='dashed', group=1)+
  geom_point()+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle=90)
  )+
  ylab("Variance explained (%)")+
  xlab("Principal component")
kk.c <- ggplot(pca.kk.var %>% filter(PC.n<=n.pc), aes(x=PC, y=cumvar))+
  geom_line(linetype='dashed', group=1)+
  geom_point()+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle=90)
  )+
  ylab("Cumulative variance explained (%)")+
  xlab("Principal component")

## Break out PC coords
pca.kk.coord <- as.data.frame(pca.kk$x) %>% 
  rownames_to_column('replicate') %>%
  left_join(rep2cond, by='replicate')
pca.kk.coord$condition <- factor(pca.kk.coord$condition, levels=kkeep)

## Break out PC loadings coords
pca.kk.load <- as.data.frame(pca.kk$rotation) %>% 
  rownames_to_column('Compounds.ID')

## Plot PC1-2
kk.12 <- ggplot(pca.kk.coord, aes(x=PC1, y=PC2, fill=condition))+
  geom_point(shape=23, color='black', alpha=.8, size=3)+
  #geom_text(aes(label=replicate))+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    legend.position = c(.8,.15),
    legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
    legend.key.size = unit(2.5, 'mm')
  )+
  xlab(pca.kk.vmap[pca.kk.vmap$PC=='PC1','ve.label'])+
  ylab(pca.kk.vmap[pca.kk.vmap$PC=='PC2','ve.label'])+
  scale_fill_brewer(palette='Set2')+
  guides(fill=guide_legend(ncol=2))

## Plot PC1-2 loadings
kk.lab <- pca.kk.load %>% 
  filter(Compounds.ID %in% knowns$Compounds.ID) %>% 
  left_join(knowns, by='Compounds.ID')
kkl.12 <- ggplot(pca.kk.load %>% filter(!(Compounds.ID %in% media.MF$Compounds.ID)), aes(x=PC1, y=PC2))+
  geom_point(alpha=.1)+
  geom_point(data=kk.lab, color='red')+
  geom_label_repel(data=kk.lab, aes(label = mf.name),
                   segment.color = 'red')+
  theme_bw()+
  xlab(pca.kk.vmap[pca.kk.vmap$PC=='PC1','ve.label']) +
  ylab(pca.kk.vmap[pca.kk.vmap$PC=='PC2','ve.label'])

plot_grid(kk.v, kk.c, kk.12, kkl.12,
          labels=LETTERS[1:4],
          nrow=2
)
ggsave("kk_normunitscale_area_pc.pdf", height=8, width=8, dpi=1600)

## Explore PCs
kk13 <- ggplot(pca.kk.coord, aes(x=PC1, y=PC3, fill=condition))+
  geom_point(shape=23, color='black', alpha=.8, size=3)+
  #geom_text(aes(label=replicate))+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    legend.position = c(.2,.4),
    legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
    legend.key.size = unit(2.5, 'mm')
  )+
  xlab(pca.kk.vmap[pca.kk.vmap$PC=='PC1','ve.label'])+
  ylab(pca.kk.vmap[pca.kk.vmap$PC=='PC3','ve.label'])+
  scale_fill_brewer(palette='Set2')+
  guides(fill=guide_legend(ncol=2))

kk13l <- ggplot(pca.kk.load %>% filter(!(Compounds.ID %in% media.MF$Compounds.ID)), aes(x=PC1, y=PC3))+
  geom_point(alpha=.1)+
  geom_point(data=kk.lab, color='red')+
  geom_label_repel(data=kk.lab, aes(label = mf.name),
                   segment.color = 'red')+
  theme_bw()+
  xlab(pca.kk.vmap[pca.kk.vmap$PC=='PC1','ve.label']) +
  ylab(pca.kk.vmap[pca.kk.vmap$PC=='PC3','ve.label'])

plot_grid(kk13, kk13l,
          labels=LETTERS[1:2],
          nrow=1
)
ggsave("kk_normunitscale_area_pc13.pdf", height=4, width=8, dpi=1600)


## Calculate an all vs all Tukey test for each Compounds.ID / conditions pair
tuk.per.row <- function(x){
  all.aov <- aov(peak.area.norm ~ condition,
                 data = area.norm.tall %>% filter(Compounds.ID==x))
  all.tuk <- as.data.frame(TukeyHSD(all.aov, conf.level=.95)$condition) %>%
    select(`p adj`) %>% 
    rownames_to_column('comparison') %>% 
    mutate(Compounds.ID=x)
  return(all.tuk)
}

# all.v.all.tuk <- map_dfr(area.raw$Compounds.ID, tuk.per.row) %>% 
#   spread(key=comparison, value=`p adj`)
# all.v.all.tuk <- mid.lookup %>% 
#   right_join(all.v.all.tuk, by='Compounds.ID')
# write.table(all.v.all.tuk, "all_v_all_tuk.tsv", sep="\t", row.names = F)

all.v.all.tuk <- read.table("all_v_all_tuk.tsv", sep="\t", header=T)

ava.tall <- all.v.all.tuk %>% 
  select(-Name, -Calc_MW, -RT_min, -Formula) %>% 
  gather(key="comparison", value='p.adj', -Compounds.ID)
comps <- unique(ava.tall$comparison)

## Example...BFK signif
sig.cut <- 0.01
bfk.sig <- mid.lookup %>% 
  right_join(ava.tall, by="Compounds.ID") %>%
  filter(grepl("^BFK-", comparison) | grepl("-BFK$", comparison)) %>% 
  filter(p.adj < sig.cut) %>% 
  spread(key=comparison, value=p.adj) %>% 
  na.omit()
write.table(bfk.sig, "bfk_signif_mf_list.csv", sep=',', quote=F, row.names = F)

bfd.sig <- mid.lookup %>% 
  right_join(ava.tall, by="Compounds.ID") %>%
  filter(grepl("^BFD-", comparison) | grepl("-BFD$", comparison)) %>% 
  filter(p.adj < sig.cut) %>% 
  spread(key=comparison, value=p.adj) %>% 
  na.omit()
write.table(bfd.sig, "bfd_signif_mf_list.csv", sep=',', quote=F, row.names = F)

## transat
mwcut <- sum(117.1, 165.2, 133.1, 132.1, 75.1, 147.1)
tat <- mid.lookup %>% 
  right_join(ava.tall, by="Compounds.ID") %>%
  filter(comparison %in% c('BK-B','BFK-B','BK-BF','BFK-BF','K-B','F-B','K-BF','F-BF','FK-B','FK-BF')) %>% 
  filter(p.adj < sig.cut) %>% 
  spread(key=comparison, value=p.adj) %>% 
  na.omit() %>% 
  arrange(-Calc_MW) %>% 
  filter(Calc_MW>=mwcut)
write.table(tat, "transat_signif_mf_list.csv", sep=',', quote=F, row.names = F)


pot_kur <- c('c3213', 'c2205', 'c348', 'c466', 'c1769')
tmp <- all.v.all.tuk %>% 
  filter(Compounds.ID %in% pot_kur)
tmp <- u.pres %>% 
  filter(Compounds.ID %in% pot_kur)
  