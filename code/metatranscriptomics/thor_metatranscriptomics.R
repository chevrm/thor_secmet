set.seed(8675309)

library(tidyverse)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(cowplot)
library(lfc)
library(openxlsx)
library(gggenes)
library(ggridges)


## Define base colors
b.col <- '#6ea8ca'
f.col <- '#ca9c6e'
k.col <- '#ca6e6e'

## Read in htseq count files
wt.htseq <- read.table("htseq_all.tsv", sep="\t", header=T) %>% 
  mutate(organism = case_when(
    organism == 'Bc' ~ 'Bcer',
    organism == 'Fj' ~ 'Fjoh',
    organism == 'Pk' ~ 'Pkor',
    T ~ 'ERROR'
  ), condition = case_when(
    condition == 'B' ~ 'Bcer',
    condition == 'F' ~ 'Fjoh',
    condition == 'K' ~ 'Pkor',
    condition == 'BF' ~ 'Bcer.Fjoh',
    condition == 'BK' ~ 'Bcer.Pkor',
    condition == 'FK' ~ 'Fjoh.Pkor',
    condition == 'BFK' ~ 'Bcer.Fjoh.Pkor',
    T ~ 'ERROR'
  ))
delkec.htseq <- read.table("htseq_delkec.tsv", sep="\t", header=T) %>% 
  mutate(organism = case_when(
    organism == 'Bc' ~ 'Bcer',
    organism == 'Fj' ~ 'Fjoh',
    organism == 'Pk' ~ 'Pkor',
    T ~ 'ERROR'
  ), condition = case_when(
    condition == 'BdelkecK' ~ 'Bcer.Pkor*',
    condition == 'BFdelkecK' ~ 'Bcer.Fjoh.Pkor*',
    condition == 'BFK' ~ 'Bcer.Fjoh.Pkor',
    condition == 'BK' ~ 'Bcer.Pkor',
    condition == 'delkecK' ~ 'Pkor*',
    condition == 'FdelkecK' ~ 'Fjoh.Pkor*',
    condition == 'FK' ~ 'Fjoh.Pkor',
    condition == 'FloK' ~ 'Fjoh.Pkor_li',
    condition == 'K' ~ 'Pkor',
    T ~ 'ERROR'
  ))

## Read in gene annotations
g2a <- read.table("gene2bgc.tsv", sep="\t", header=T, na.strings = 'not_in_bgc')
g2c <- read.table("gene2cog.tsv", sep="\t", header=T, na.strings = 'No_COG')
g2gk <- read.table("all.as_gk.tsv", sep="\t", header=T)
bg <- read.xlsx("bgc_ann.xlsx", sheet=2)

## Read in cluster annotations
as <- read.table("antismash_summary.tsv", sep="\t", header=T)
man_bgc <- read.xlsx("bgc_ann.xlsx", sheet=1)

## Read in other maps
cogmap <- read.table("cog_map.csv", sep=',', header=T, col.names = c('cog_large_category', 'cog_category', 'cog_description'))

## Make master gene annotation df
gene_ann <- g2c %>% 
  mutate(organism = sub("^(\\w\\w).+", "\\1", gene)) %>% 
  mutate(is.q = case_when(
    grepl("Q", cog_category) ~ T,
    T ~ F
  )) %>% 
  mutate(organism = case_when(
    organism == 'Bc' ~ 'Bcer',
    organism == 'Fj' ~ 'Fjoh',
    organism == 'Pk' ~ 'Pkor',
    T ~ as.character(NA)
  )) %>% 
  left_join(g2a %>% select(gene, bgc), by='gene') %>% 
  left_join(g2gk, by='gene') %>% 
  mutate(gene_kind = case_when(
    is.na(gene_kind) & !(is.na(bgc)) ~ 'other',
    T ~ as.character(gene_kind)
  )) %>% 
  mutate(gene_kind = case_when(
    gene_kind == 'biosynthetic' ~ 'Core biosynthesis',
    gene_kind == 'biosynthetic-additional' ~ 'Ancillary biosynthesis',
    gene_kind == 'transport' ~ 'Transport',
    gene_kind == 'regulatory' ~ 'Regulation',
    gene_kind == 'other' ~ 'Other',
    T ~ as.character(NA)
  ))

## Make maps for later
wt.samp2cond <- wt.htseq %>% 
  select(sample, condition) %>% 
  distinct()
delkec.samp2cond <- delkec.htseq %>% 
  select(sample, condition) %>% 
  distinct()
gene2org <- gene_ann %>% 
  select(gene, organism) %>% 
  distinct()

## Define WT condition groups
bwt.cond <- c('Bcer', 'Bcer.Fjoh', 'Bcer.Pkor', 'Bcer.Fjoh.Pkor')
fwt.cond <- c('Fjoh', 'Bcer.Fjoh', 'Fjoh.Pkor', 'Bcer.Fjoh.Pkor')
kwt.cond <- c('Pkor', 'Bcer.Pkor', 'Fjoh.Pkor', 'Bcer.Fjoh.Pkor')

## Create matricies for edgeR processing
bwt.e <- wt.htseq %>%
  filter(organism == 'Bcer' & condition %in% bwt.cond) %>% 
  select(gene, sample, htseq_count) %>% 
  spread(key=sample, value=htseq_count) %>% 
  column_to_rownames('gene')
fwt.e <- wt.htseq %>%
  filter(organism == 'Fjoh' & condition %in% fwt.cond) %>% 
  select(gene, sample, htseq_count) %>% 
  spread(key=sample, value=htseq_count) %>% 
  column_to_rownames('gene')
kwt.e <- wt.htseq %>%
  filter(organism == 'Pkor' & condition %in% kwt.cond) %>% 
  select(gene, sample, htseq_count) %>% 
  spread(key=sample, value=htseq_count) %>% 
  column_to_rownames('gene')

## Calculate CPMs
bwt.grp <- sub("\\d$", "", colnames(bwt.e))
bwt.d <- DGEList(counts=bwt.e, group=factor(bwt.grp))
bwt.d <- calcNormFactors(bwt.d) ## TMM norm
bwt.cpm <- as.data.frame(cpm(bwt.d)) %>% 
  rownames_to_column('gene') %>% 
  gather(key='sample', value='cpm', -gene)
fwt.grp <- sub("\\d$", "", colnames(fwt.e))
fwt.d <- DGEList(counts=fwt.e, group=factor(fwt.grp))
fwt.d <- calcNormFactors(fwt.d) ## TMM norm
fwt.cpm <- as.data.frame(cpm(fwt.d)) %>% 
  rownames_to_column('gene') %>% 
  gather(key='sample', value='cpm', -gene)
kwt.grp <- sub("\\d$", "", colnames(kwt.e))
kwt.d <- DGEList(counts=kwt.e, group=factor(kwt.grp))
kwt.d <- calcNormFactors(kwt.d) ## TMM norm
kwt.cpm <- as.data.frame(cpm(kwt.d)) %>% 
  rownames_to_column('gene') %>% 
  gather(key='sample', value='cpm', -gene)
wt.cpm <- rbind(bwt.cpm, fwt.cpm, kwt.cpm) %>% 
  left_join(wt.samp2cond, by='sample') %>% 
  left_join(gene_ann, by='gene')

## Calculate each gene's maximum per condition mean
## to define lower cutoff
wt.ptile.cutoff <- 0.025
wt.cond.max <- wt.cpm %>% 
  group_by(gene, condition, organism) %>% 
  summarize(gmean=mean(cpm)) %>% 
  group_by(gene, organism) %>% 
  summarize(gmax=max(gmean)) %>% 
  group_by(organism) %>% 
  mutate(ptile = percent_rank(gmax))
wt.line <- wt.cond.max %>% 
  filter(ptile > wt.ptile.cutoff) %>% 
  arrange(ptile) %>% 
  distinct(organism, .keep_all = T) %>% 
  select(organism, gmax) %>% 
  rename(cutoff=gmax)

## Plot the distributions and cutoff
bwt.cp <- ggplot(wt.cond.max %>% filter(organism=='Bcer'), aes(x=gmax))+
  geom_vline(xintercept = as.numeric(wt.line[wt.line$organism=='Bcer','cutoff']), linetype='dashed') +
  geom_density(color='#333333', alpha=0.3, size=1, fill=b.col)+
  scale_x_log10(limits=c(1e-3, 1e6))+
  theme_bw()+theme(
    legend.title = element_blank(),
    legend.position = 'none'
  )+
  ylab('Frequency of genes')+xlab("Maximum per condition mean CPM")+
  ggtitle('Bcer')+
  scale_y_continuous(limits=c(0,.5))
fwt.cp <- ggplot(wt.cond.max %>% filter(organism=='Fjoh'), aes(x=gmax))+
  geom_vline(xintercept = as.numeric(wt.line[wt.line$organism=='Fjoh','cutoff']), linetype='dashed') +
  geom_density(color='#333333', alpha=0.3, size=1, fill=f.col)+
  scale_x_log10(limits=c(1e-3, 1e6))+
  theme_bw()+theme(
    legend.title = element_blank(),
    legend.position = 'none'
  )+
  ylab('Frequency of genes')+xlab("Maximum per condition mean CPM")+
  ggtitle('Fjoh')+
  scale_y_continuous(limits=c(0,.5))
kwt.cp <- ggplot(wt.cond.max %>% filter(organism=='Pkor'), aes(x=gmax))+
  geom_vline(xintercept = as.numeric(wt.line[wt.line$organism=='Pkor','cutoff']), linetype='dashed') +
  geom_density(color='#333333', alpha=0.3, size=1, fill=k.col)+
  scale_x_log10(limits=c(1e-3, 1e6))+
  theme_bw()+theme(
    legend.title = element_blank(),
    legend.position = 'none'
  )+
  ylab('Frequency of genes')+xlab("Maximum per condition mean CPM")+
  ggtitle('Pkor')+
  scale_y_continuous(limits=c(0,.5))
plot_grid(bwt.cp, fwt.cp, kwt.cp,
          labels=LETTERS[1:3],
          nrow = 1)
ggsave("wt_max_cpm.pdf", dpi=1600, height=4, width=12)

## Filter for genes above the cutoff
wt.keep <- wt.cond.max %>% 
  left_join(wt.line, by='organism') %>% 
  filter(gmax >= cutoff)
wt.cpm.cut <- wt.cpm %>% 
  filter(gene %in% wt.keep$gene)
write.table(wt.cpm.cut, "wt_cpm.tsv", sep="\t", row.names = F, quote = F)

## Calculate fold changes
bwt.fc <- wt.cpm.cut %>% 
  filter(organism=='Bcer') %>% 
  group_by(gene, condition) %>% 
  summarize(mean.cpm=mean(cpm)) %>% 
  spread(key = condition, value=mean.cpm) %>% 
  select(gene, all_of(bwt.cond)) %>% 
  ungroup() %>% 
  mutate(Bcer.Fjoh_LFC = PsiLFC(Bcer.Fjoh, Bcer),
         Bcer.Pkor_LFC = PsiLFC(Bcer.Pkor, Bcer),
         Bcer.Fjoh.Pkor_LFC = PsiLFC(Bcer.Fjoh.Pkor, Bcer))
fwt.fc <- wt.cpm.cut %>% 
  filter(organism=='Fjoh') %>% 
  group_by(gene, condition) %>% 
  summarize(mean.cpm=mean(cpm)) %>% 
  spread(key = condition, value=mean.cpm) %>% 
  select(gene, all_of(fwt.cond)) %>% 
  ungroup() %>% 
  mutate(Bcer.Fjoh_LFC = PsiLFC(Bcer.Fjoh, Fjoh),
         Fjoh.Pkor_LFC = PsiLFC(Fjoh.Pkor, Fjoh),
         Bcer.Fjoh.Pkor_LFC = PsiLFC(Bcer.Fjoh.Pkor, Fjoh))
kwt.fc <- wt.cpm.cut %>% 
  filter(organism=='Pkor') %>% 
  group_by(gene, condition) %>% 
  summarize(mean.cpm=mean(cpm)) %>% 
  spread(key = condition, value=mean.cpm) %>% 
  select(gene, all_of(kwt.cond)) %>% 
  ungroup() %>% 
  mutate(Bcer.Pkor_LFC = PsiLFC(Bcer.Pkor, Pkor),
         Fjoh.Pkor_LFC = PsiLFC(Fjoh.Pkor, Pkor),
         Bcer.Fjoh.Pkor_LFC = PsiLFC(Bcer.Fjoh.Pkor, Pkor))

## Make comparison plots
bwt.comp <- rbind(
  bwt.fc %>% 
    mutate(set='Genome-wide'),
  bwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$is.q==T,'gene']) %>% 
    mutate(set='COG-Q'),
  bwt.fc %>% 
    filter(gene %in% gene_ann[!(is.na(gene_ann$bgc)),'gene']) %>% 
    mutate(set='All antiSMASH BGC genes'),
  bwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Core biosynthesis','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Core biosynthesis'),
  bwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Ancillary biosynthesis','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Ancillary biosynthesis'),
  bwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Regulation','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Regulation'),
  bwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Transport','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Transport'),
  bwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Other','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Other')
) %>% 
  select(gene, set, Bcer.Fjoh_LFC, Bcer.Pkor_LFC, Bcer.Fjoh.Pkor_LFC) %>% 
  gather(key='comparison', value='log2FC', -gene, -set) %>% 
  mutate(comparison = sub("_LFC$", "", comparison))
bwt.comp$set <- factor(bwt.comp$set, levels = c('Genome-wide', 'COG-Q', 'All antiSMASH BGC genes', 'antiSMASH BGC genes - Core biosynthesis', 'antiSMASH BGC genes - Ancillary biosynthesis', 'antiSMASH BGC genes - Regulation', 'antiSMASH BGC genes - Transport', 'antiSMASH BGC genes - Other'))
bwt.comp$comparison <- factor(bwt.comp$comparison, levels = c('Bcer.Fjoh', 'Bcer.Pkor', 'Bcer.Fjoh.Pkor'))
ggplot(bwt.comp, aes(x=comparison, y=log2FC))+
  geom_hline(yintercept = 0, color='#bbbbbb')+
  geom_jitter(alpha=.3, width=.25, aes(color='0'))+
  geom_boxplot(outlier.shape = NA, alpha=0, width=0.2)+
  theme_bw()+theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.title.y = element_text(face='bold'),
    axis.text.x = element_text(face='bold')
  )+
  scale_color_manual(values=c(b.col))+
  facet_wrap(~set)+
  ylab(expression(Psi*' log2FC vs. Bcer alone'))
ggsave("wt_lfc.bcer.pdf", dpi=1600, height=10, width=10)

fwt.comp <- rbind(
  fwt.fc %>% 
    mutate(set='Genome-wide'),
  fwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$is.q==T,'gene']) %>% 
    mutate(set='COG-Q'),
  fwt.fc %>% 
    filter(gene %in% gene_ann[!(is.na(gene_ann$bgc)),'gene']) %>% 
    mutate(set='All antiSMASH BGC genes'),
  fwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Core biosynthesis','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Core biosynthesis'),
  fwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Ancillary biosynthesis','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Ancillary biosynthesis'),
  fwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Regulation','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Regulation'),
  fwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Transport','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Transport'),
  fwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Other','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Other')
) %>% 
  select(gene, set, Bcer.Fjoh_LFC, Fjoh.Pkor_LFC, Bcer.Fjoh.Pkor_LFC) %>% 
  gather(key='comparison', value='log2FC', -gene, -set) %>% 
  mutate(comparison = sub("_LFC$", "", comparison))
fwt.comp$set <- factor(fwt.comp$set, levels = c('Genome-wide', 'COG-Q', 'All antiSMASH BGC genes', 'antiSMASH BGC genes - Core biosynthesis', 'antiSMASH BGC genes - Ancillary biosynthesis', 'antiSMASH BGC genes - Regulation', 'antiSMASH BGC genes - Transport', 'antiSMASH BGC genes - Other'))
fwt.comp$comparison <- factor(fwt.comp$comparison, levels = c('Bcer.Fjoh', 'Fjoh.Pkor', 'Bcer.Fjoh.Pkor'))
ggplot(fwt.comp, aes(x=comparison, y=log2FC))+
  geom_hline(yintercept = 0, color='#bbbbbb')+
  geom_jitter(alpha=.3, width=.25, aes(color='0'))+
  geom_boxplot(outlier.shape = NA, alpha=0, width=0.2)+
  theme_bw()+theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.title.y = element_text(face='bold'),
    axis.text.x = element_text(face='bold')
  )+
  scale_color_manual(values=c(f.col))+
  facet_wrap(~set)+
  ylab(expression(Psi*' log2FC vs. Fjoh alone'))
ggsave("wt_lfc.fjoh.pdf", dpi=1600, height=10, width=10)

kwt.comp <- rbind(
  kwt.fc %>% 
    mutate(set='Genome-wide'),
  kwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$is.q==T,'gene']) %>% 
    mutate(set='COG-Q'),
  kwt.fc %>% 
    filter(gene %in% gene_ann[!(is.na(gene_ann$bgc)),'gene']) %>% 
    mutate(set='All antiSMASH BGC genes'),
  kwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Core biosynthesis','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Core biosynthesis'),
  kwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Ancillary biosynthesis','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Ancillary biosynthesis'),
  kwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Regulation','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Regulation'),
  kwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Transport','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Transport'),
  kwt.fc %>% 
    filter(gene %in% gene_ann[gene_ann$gene_kind=='Other','gene']) %>% 
    mutate(set='antiSMASH BGC genes - Other')
) %>% 
  select(gene, set, Bcer.Pkor_LFC, Fjoh.Pkor_LFC, Bcer.Fjoh.Pkor_LFC) %>% 
  gather(key='comparison', value='log2FC', -gene, -set) %>% 
  mutate(comparison = sub("_LFC$", "", comparison))
kwt.comp$set <- factor(kwt.comp$set, levels = c('Genome-wide', 'COG-Q', 'All antiSMASH BGC genes', 'antiSMASH BGC genes - Core biosynthesis', 'antiSMASH BGC genes - Ancillary biosynthesis', 'antiSMASH BGC genes - Regulation', 'antiSMASH BGC genes - Transport', 'antiSMASH BGC genes - Other'))
kwt.comp$comparison <- factor(kwt.comp$comparison, levels = c('Bcer.Pkor', 'Fjoh.Pkor', 'Bcer.Fjoh.Pkor'))
ggplot(kwt.comp, aes(x=comparison, y=log2FC))+
  geom_hline(yintercept = 0, color='#bbbbbb')+
  geom_jitter(alpha=.3, width=.25, aes(color='0'))+
  geom_boxplot(outlier.shape = NA, alpha=0, width=0.2)+
  theme_classic()+theme(
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.title.y = element_text(face='bold'),
    axis.text.x = element_text(face='bold')
  )+
  scale_color_manual(values=c(k.col))+
  facet_wrap(~set)+
  ylab(expression(Psi*' log2FC vs. Pkor alone'))
ggsave("wt_lfc.pkor.pdf", dpi=1600, height=10, width=10)

## Summary table for zoomed out LFCs 
bst <- bwt.comp %>%
  group_by(set, comparison) %>% 
  summarize(n=n(),
            LFCmean = mean(log2FC),
            LFCmedian = median(log2FC),
            LFCstdev = sd(log2FC),
            LFC25pctile = quantile(log2FC, probs = 0.25),
            LFC75pctile = quantile(log2FC, probs = 0.75),
            LFCiqr = LFC75pctile - LFC25pctile,
            LFCmin = min(log2FC),
            LFCmax = max(log2FC)
            ) %>% 
  mutate(organism='Bcer') %>% 
  select(12, 1:11)
fst <- fwt.comp %>%
  group_by(set, comparison) %>% 
  summarize(n=n(),
            LFCmean = mean(log2FC),
            LFCmedian = median(log2FC),
            LFCstdev = sd(log2FC),
            LFC25pctile = quantile(log2FC, probs = 0.25),
            LFC75pctile = quantile(log2FC, probs = 0.75),
            LFCiqr = LFC75pctile - LFC25pctile,
            LFCmin = min(log2FC),
            LFCmax = max(log2FC)
  ) %>% 
  mutate(organism='Fjoh') %>% 
  select(12, 1:11)
kst <- kwt.comp %>%
  group_by(set, comparison) %>% 
  summarize(n=n(),
            LFCmean = mean(log2FC),
            LFCmedian = median(log2FC),
            LFCstdev = sd(log2FC),
            LFC25pctile = quantile(log2FC, probs = 0.25),
            LFC75pctile = quantile(log2FC, probs = 0.75),
            LFCiqr = LFC75pctile - LFC25pctile,
            LFCmin = min(log2FC),
            LFCmax = max(log2FC)
  ) %>%   mutate(organism='Pkor') %>% 
  select(12, 1:11)
wtst <- rbind(bst, fst, kst)
wtst$comparison <- factor(wtst$comparison, levels=c('Bcer.Fjoh', 'Bcer.Pkor', 'Fjoh.Pkor', 'Bcer.Fjoh.Pkor'))
stp <- ggplot(wtst, aes(x=organism, y=LFCiqr))+
  geom_bar(stat='identity', position='dodge', color='gray30', aes(fill=comparison))+
  theme_bw()+
  theme(
    axis.text.x = element_text(face='italic'),
    legend.text = element_text(face='italic'),
    legend.position = 'top'
  )+
  ylab(expression('IQR of '*Psi*' log2FC vs. organism alone'))+
  xlab('Organism')+
  scale_fill_brewer(palette='Set2', name='Comparison')+
  facet_wrap(~set)

## Summarize the different LFC data
btm <- bwt.comp %>% 
  mutate(organism='Bcer')
ftm <- fwt.comp %>% 
  mutate(organism='Fjoh')
ktm <- kwt.comp %>% 
  mutate(organism='Pkor')
wttm <- rbind(btm, ftm, ktm)
wttm$comparison <- factor(wttm$comparison, levels=c('Bcer.Fjoh', 'Bcer.Pkor', 'Fjoh.Pkor', 'Bcer.Fjoh.Pkor'))
ttp <- ggplot(wttm, aes(x=organism, y=log2FC))+
  geom_point(position=position_jitterdodge(), alpha=.3, aes(color=comparison))+
  geom_hline(yintercept = 0, size=0.5)+
  geom_boxplot(outlier.shape = NA, color='gray30', aes(fill=comparison))+
  theme_bw()+
  theme(
    axis.text.x = element_text(face='italic'),
    legend.text = element_text(face='italic'),
    legend.position = 'top'
  )+
  ylab(expression(Psi*' log2FC vs. organism alone'))+
  xlab('Organism')+
  #scale_y_continuous(limits=c(-3.5,3.5))+
  scale_color_brewer(palette='Set2', name='Comparison')+
  scale_fill_brewer(palette='Set2', name='Comparison')+
  facet_wrap(~set)
plot_grid(ttp, stp,
          labels=LETTERS[1:2],
          nrow=2)
ggsave('allwtLFCsum.pdf', width=9, height=14, dpi=1600)

## Gather all gene info
gff.head <- c('contig', 'ann_method', 'locus_type', 'start', 'end', 'score', 'strand', 'phase', 'ann')
all.gff <- rbind(
  read.table("Bc.gff", header=F, sep="\t", col.names = gff.head) %>% 
    mutate(gene=sub("^ID=(.+);$", "Bc_\\1", ann)),
  read.table("Fj.gff", header=F, sep="\t", col.names = gff.head) %>% 
    mutate(gene=sub("^ID=(.+);$", "Fj_\\1", ann)),
  read.table("Pk.gff", header=F, sep="\t", col.names = gff.head) %>% 
    mutate(gene=sub("^ID=(.+);$", "Pk_\\1", ann))
) %>% left_join(gene_ann, by='gene') %>% 
  mutate(strand=case_when(
    strand=='+' ~ 1,
    strand=='-' ~ -1
  ))


## Let's make some gene maps!
wt.c.order <- c('Bcer', 'Fjoh', 'Pkor',
                'Bcer.Fjoh', 'Bcer.Pkor', 'Fjoh.Pkor', 'Bcer.Fjoh.Pkor')

genemap.plot <- function(bgc_locus){
  ## Determine organism and pull proper FC df
  if(grepl("^Bc", bgc_locus)){
    bgc.wt.fc <- bwt.fc  
    bgc.org <- 'Bcer'
  }else if(grepl("^Fj", bgc_locus)){
    bgc.wt.fc <- fwt.fc
    bgc.org <- 'Fjoh'
  }else if(grepl("^Pk", bgc_locus)){
    bgc.wt.fc <- kwt.fc
    bgc.org <- 'Pkor'
  }else{
    bgc.wt.fc <- 'ERROR'
    bgc.org <- 'ERROR'
  }
  ## Subset down to the BGC
  if(bgc_locus %in% bg$replace_antismash_bgc){
    bgc_genes <- bg %>% 
      filter(replace_antismash_bgc==bgc_locus) %>%
      select(gene, ann) %>% 
      distinct(gene, .keep_all=T) %>% 
      rename(lab=ann)
  }else if(bgc_locus %in% gene_ann$bgc){
    bgc_genes <- gene_ann %>% 
      filter(bgc==bgc_locus) %>%
      select(gene) %>% 
      distinct() %>% 
      mutate(lab=gene)
  }else{
    bgc_genes = 'ERROR'
  }
  bgc.wt.fc <- bgc_genes %>%
    left_join(bgc.wt.fc, by='gene')
  ## Get the expression alone df
  if(bgc.org=='Bcer'){
    bgc.alone <- bgc.wt.fc %>% 
      select(gene, lab, Bcer)
  }else if(bgc.org=='Fjoh'){
    bgc.alone <- bgc.wt.fc %>% 
      select(gene, lab, Fjoh)
  }else if(bgc.org=='Pkor'){
    bgc.alone <- bgc.wt.fc %>% 
      select(gene, lab, Pkor)
  }else{
    bgc.alone <- 'ERROR'
  }
  bgc.alone <- bgc.alone %>% 
    left_join(all.gff, by='gene')
  ## Adjust x axis to start counting at BGC start
  bgc.min <- min(bgc.alone$start)
  bgc.alone <- bgc.alone %>% 
    mutate(start.adj = start-bgc.min+1, end.adj=end-bgc.min+1)
  ## Plot alone expression
  if(bgc.org=='Bcer'){
    bap <- ggplot(bgc.alone, aes(xmin=start.adj, xmax=end.adj, y=bgc.org)) +
      geom_gene_arrow(aes(fill=Bcer, forward=strand),
                      arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm"))+
      geom_gene_label(align='left', aes(label=lab), color='#FFFFFF')+
      theme_genes()+
      theme(
        axis.title.y = element_blank(),
        legend.key.size = unit(3, 'mm')
      )+
      scale_fill_gradient(name='CPM')
  }else if(bgc.org=='Fjoh'){
    bap <- ggplot(bgc.alone, aes(xmin=start.adj, xmax=end.adj, y=bgc.org)) +
      geom_gene_arrow(aes(fill=Fjoh, forward=strand),
                      arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm"))+
      geom_gene_label(align='left', aes(label=lab), color='#FFFFFF')+
      theme_genes()+
      theme(
        axis.title.y = element_blank(),
        legend.key.size = unit(3, 'mm')
      )+
      scale_fill_gradient(name='CPM')
  }else if(bgc.org=='Pkor'){
    bap <- ggplot(bgc.alone, aes(xmin=start.adj, xmax=end.adj, y=bgc.org)) +
      geom_gene_arrow(aes(fill=Pkor, forward=strand),
                      arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm"))+
      geom_gene_label(align='left', aes(label=lab), color='#FFFFFF')+
      theme_genes()+
      theme(
        axis.title.y = element_blank(),
        legend.key.size = unit(3, 'mm')
      )+
      scale_fill_gradient(name='CPM')
  }else{
    bap <- 'ERROR'
  }
  ## Create lfc df  
  bgc.lfc.t <- bgc.wt.fc %>% 
    gather(key='comparison', value='lfc', -gene, -lab) %>% 
    filter(grepl("_LFC$", comparison)) %>% 
    mutate(comparison=sub("_LFC$", "", comparison)) %>% 
    left_join(all.gff, by='gene') %>% 
    mutate(start.adj = start-bgc.min+1, end.adj=end-bgc.min+1)
  bgc.lfc.t$comparison <- factor(bgc.lfc.t$comparison, levels=wt.c.order)
  ## Plot LFC
  blp <- ggplot(bgc.lfc.t, aes(xmin=start.adj, xmax=end.adj, y=comparison)) +
    geom_gene_arrow(aes(fill=lfc, forward=strand),
                    arrowhead_height = unit(3, "mm"),
                    arrowhead_width = unit(1, "mm"))+
    theme_genes()+
    theme(
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.key.size = unit(3, 'mm')
    )+
    facet_wrap(~ comparison, scales='free', ncol=1) +
    scale_fill_gradient2(midpoint=0,
                         low = RColorBrewer::brewer.pal(11, "BrBG")[1:5],
                         mid = RColorBrewer::brewer.pal(11, "BrBG")[6],
                         high = RColorBrewer::brewer.pal(11, "BrBG")[7:11],
                         #low = RColorBrewer::brewer.pal(11, "BrBG")[2],
                         #high = RColorBrewer::brewer.pal(11, "BrBG")[10],
                         #high='seagreen', low='hotpink',
                         name=expression(Psi*' log2FC')
                        )
  ## Generate multipanel plot and return
  bgc.grd <- plot_grid(bap, blp,
                       ncol=1, align = 'v',
                       labels = c('A', 'B'),
                       rel_heights = c(1, 1.25))
  return(bgc.grd)
}

get_wt_fc_genes <- function(gene_list){
  ## Determine organism and pull proper FC df
  if(grepl("^Bc", gene_list[1])){
    bgc.wt.fc <- bwt.fc
  }else if(grepl("^Fj", gene_list[1])){
    bgc.wt.fc <- fwt.fc
  }else if(grepl("^Pk", gene_list[1])){
    bgc.wt.fc <- kwt.fc
  }else{
    bgc.wt.fc <- 'ERROR'
  }
  return(bgc.wt.fc %>% filter(gene %in% gene_list))
}

genemap.grid <- function(bgc_locus){
  ## Determine organism and pull proper FC df
  if(grepl("^Bc", bgc_locus)){
    bgc.wt.fc <- bwt.fc  
    bgc.org <- 'Bcer'
  }else if(grepl("^Fj", bgc_locus)){
    bgc.wt.fc <- fwt.fc
    bgc.org <- 'Fjoh'
  }else if(grepl("^Pk", bgc_locus)){
    bgc.wt.fc <- kwt.fc
    bgc.org <- 'Pkor'
  }else{
    bgc.wt.fc <- 'ERROR'
    bgc.org <- 'ERROR'
  }
  ## Subset down to the BGC
  if(bgc_locus %in% bg$replace_antismash_bgc){
    bgc_genes <- bg %>% 
      filter(replace_antismash_bgc==bgc_locus) %>%
      select(gene, ann) %>% 
      distinct(gene, .keep_all=T) %>% 
      rename(lab=ann)
  }else if(bgc_locus %in% gene_ann$bgc){
    bgc_genes <- gene_ann %>% 
      filter(bgc==bgc_locus) %>%
      select(gene) %>% 
      distinct() %>% 
      mutate(lab=gene)
  }else{
    bgc_genes = 'ERROR'
  }
  bgc.wt.fc <- bgc_genes %>%
    left_join(bgc.wt.fc, by='gene')
  ## Get the expression alone df
  if(bgc.org=='Bcer'){
    bgc.alone <- bgc.wt.fc %>% 
      select(gene, lab, Bcer)
  }else if(bgc.org=='Fjoh'){
    bgc.alone <- bgc.wt.fc %>% 
      select(gene, lab, Fjoh)
  }else if(bgc.org=='Pkor'){
    bgc.alone <- bgc.wt.fc %>% 
      select(gene, lab, Pkor)
  }else{
    bgc.alone <- 'ERROR'
  }
  bgc.alone <- bgc.alone %>% 
    left_join(all.gff, by='gene')
  bgc.alone[is.na(bgc.alone$lab),]$lab <- bgc.alone[is.na(bgc.alone$lab),]$gene
  ## Adjust x axis to start counting at BGC start
  bgc.min <- min(bgc.alone$start)
  bgc.alone <- bgc.alone %>% 
    mutate(start.adj = start-bgc.min+1, end.adj=end-bgc.min+1)
  ## Plot alone expression
  if(bgc.org=='Bcer'){
    bap <- ggplot(bgc.alone, aes(xmin=start.adj, xmax=end.adj, y=bgc.org)) +
      geom_gene_arrow(aes(fill=Bcer, forward=strand),
                      arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm"))+
      geom_gene_label(align='left', aes(label=lab), color='#FFFFFF')+
      theme_genes()+
      theme(
        axis.title.y = element_blank(),
        legend.key.size = unit(3, 'mm')
      )+
      scale_fill_gradient(name='CPM')
  }else if(bgc.org=='Fjoh'){
    bap <- ggplot(bgc.alone, aes(xmin=start.adj, xmax=end.adj, y=bgc.org)) +
      geom_gene_arrow(aes(fill=Fjoh, forward=strand),
                      arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm"))+
      geom_gene_label(align='left', aes(label=lab), color='#FFFFFF')+
      theme_genes()+
      theme(
        axis.title.y = element_blank(),
        legend.key.size = unit(3, 'mm')
      )+
      scale_fill_gradient(name='CPM')
  }else if(bgc.org=='Pkor'){
    bap <- ggplot(bgc.alone, aes(xmin=start.adj, xmax=end.adj, y=bgc.org)) +
      geom_gene_arrow(aes(fill=Pkor, forward=strand),
                      arrowhead_height = unit(3, "mm"),
                      arrowhead_width = unit(1, "mm"))+
      geom_gene_label(align='left', aes(label=lab), color='#FFFFFF')+
      theme_genes()+
      theme(
        axis.title.y = element_blank(),
        legend.key.size = unit(3, 'mm')
      )+
      scale_fill_gradient(name='CPM')
  }else{
    bap <- 'ERROR'
  }
  ## Create lfc df  
  bgc.lfc.t <- bgc.wt.fc %>% 
    gather(key='comparison', value='lfc', -gene, -lab) %>% 
    filter(grepl("_LFC$", comparison)) %>% 
    mutate(comparison=sub("_LFC$", "", comparison)) %>% 
    left_join(all.gff, by='gene') %>% 
    mutate(start.adj = start-bgc.min+1, end.adj=end-bgc.min+1)
  bgc.lfc.t$comparison <- factor(bgc.lfc.t$comparison, levels=rev(wt.c.order))
  bgc.lfc.t[is.na(bgc.lfc.t$lab),]$lab <- bgc.lfc.t[is.na(bgc.lfc.t$lab),]$gene
  ## Plot LFC
  blg <- ggplot(bgc.lfc.t, aes(y=comparison, x=reorder(lab, start)))+
    geom_tile(aes(fill=lfc), color='white')+
    scale_fill_gradient2(midpoint=0,
                         low = RColorBrewer::brewer.pal(11, "BrBG")[1:5],
                         mid = RColorBrewer::brewer.pal(11, "BrBG")[6],
                         high = RColorBrewer::brewer.pal(11, "BrBG")[7:11],
                         #low = RColorBrewer::brewer.pal(11, "BrBG")[2],
                         #high = RColorBrewer::brewer.pal(11, "BrBG")[10],
                         #high='seagreen', low='hotpink',
                         name=expression(Psi*' '*log[2]*'FC')
    )+
    theme_bw()+
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle=60, hjust=-.1),
      #axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.key.size = unit(3, 'mm')
    )+
    coord_equal()+
    scale_x_discrete(position='top')
  ## Generate multipanel plot and return
  bgc.grd <- plot_grid(bap, blg,
                       ncol=1, align = 'v',
                       labels = c('A', 'B'))#, rel_heights = c(1, 1))
  return(bgc.grd)
}


## Testing
genemap.plot('Bc_Ga0417192_02_rn03')
genemap.grid('Bc_Ga0417192_02_rn03')

# Generate gene maps for every BGC locus
# for(bgc_locus in as$bgc){
#   #pdf.fn <- paste0('bgc_genemaps/', bgc_locus, '.pdf')
#   png.fn <- paste0('bgc_genemaps/', bgc_locus, '.png')
#   plt <- genemap.plot(bgc_locus)
#   #ggsave(plot=plt, filename = pdf.fn, dpi=1600, height=3, width=7)
#   ggsave(plot=plt, filename = png.fn, dpi=1600, height=3, width=7)
# }

# Generate gene grids for every BGC locus
# for(bgc_locus in as$bgc){
#    pdf.fn <- paste0('bgc_genegrids/', bgc_locus, '.pdf')
#    png.fn <- paste0('bgc_genegrids/', bgc_locus, '.png')
#    plt <- genemap.grid(bgc_locus)
#    ggsave(plot=plt, filename = pdf.fn, dpi=1600, height=3, width=7)
#    #ggsave(plot=plt, filename = png.fn, dpi=1600, height=3, width=7)
# }


## BGC summary stats
bgc.cpm <- wt.cpm %>% 
    filter(!(bgc %in% bg$replace_antismash_bgc)) %>% 
  rbind(
    wt.cpm %>% filter(gene %in% bg$gene)
  ) %>% 
  filter(!(is.na(bgc)))
  
bgc.ss <- bgc.cpm %>% 
  group_by(bgc, sample) %>% 
  summarize(bgc.cpm=sum(cpm)) %>% 
  left_join(
    man_bgc %>%
      select(antiSMASH_bgc, organism, size) %>% rename(bgc=antiSMASH_bgc), 
    by='bgc'
  ) %>% 
  mutate(bgc.rpkm=bgc.cpm/size*1000) %>%
  left_join(wt.samp2cond, by='sample') %>% 
  group_by(organism, bgc, condition) %>% 
  summarize(bgc.cpm.avg=mean(bgc.cpm), bgc.cpm.sd=sd(bgc.cpm),
            bgc.rpkm.avg=mean(bgc.rpkm), bgc.rpkm.sd=sd(bgc.rpkm))

fccut <- 1
bwt.ss <- bgc.ss %>% 
  select(organism, bgc, condition, bgc.cpm.avg) %>% 
  filter(organism=='Bacillus cereus UW85') %>% 
  spread(key=condition, value=bgc.cpm.avg) %>% 
  ungroup() %>% 
  mutate(Bcer.Fjoh_LFC = PsiLFC(Bcer.Fjoh, Bcer),
         Bcer.Pkor_LFC = PsiLFC(Bcer.Pkor, Bcer),
         Bcer.Fjoh.Pkor_LFC = PsiLFC(Bcer.Fjoh.Pkor, Bcer))
bwt.ss.c <- bwt.ss %>%
  gather(key='comparison', value='lfc', -organism, -bgc) %>% 
  filter(grepl("LFC$", comparison)) %>% 
  filter(abs(lfc) >= fccut)
fwt.ss <- bgc.ss %>% 
  select(organism, bgc, condition, bgc.cpm.avg) %>% 
  filter(organism=='Flavobacterium johnsoniae UW101') %>% 
  spread(key=condition, value=bgc.cpm.avg) %>% 
  ungroup() %>% 
  mutate(Bcer.Fjoh_LFC = PsiLFC(Bcer.Fjoh, Fjoh),
         Fjoh.Pkor_LFC = PsiLFC(Fjoh.Pkor, Fjoh),
         Bcer.Fjoh.Pkor_LFC = PsiLFC(Bcer.Fjoh.Pkor, Fjoh))
fwt.ss.c <- fwt.ss %>%
  gather(key='comparison', value='lfc', -organism, -bgc) %>% 
  filter(grepl("LFC$", comparison)) %>% 
  filter(abs(lfc) >= fccut)
kwt.ss <- bgc.ss %>% 
  select(organism, bgc, condition, bgc.cpm.avg) %>% 
  filter(organism=='Pseudomonas koreensis CI12') %>% 
  spread(key=condition, value=bgc.cpm.avg) %>% 
  ungroup() %>% 
  mutate(Bcer.Pkor_LFC = PsiLFC(Bcer.Pkor, Pkor),
         Fjoh.Pkor_LFC = PsiLFC(Fjoh.Pkor, Pkor),
         Bcer.Fjoh.Pkor_LFC = PsiLFC(Bcer.Fjoh.Pkor, Pkor))
kwt.ss.c <- kwt.ss %>%
  gather(key='comparison', value='lfc', -organism, -bgc) %>% 
  filter(grepl("LFC$", comparison)) %>% 
  filter(abs(lfc) >= fccut)

## Bcer, all
bwt.mat <- bwt.ss %>% 
  select(bgc, Bcer.Fjoh_LFC, Bcer.Pkor_LFC, Bcer.Fjoh.Pkor_LFC) %>% 
  column_to_rownames('bgc') %>% 
  as.matrix()
# pheatmap(bwt.mat,
#          breaks=seq(-max(abs(bwt.mat)), max(abs(bwt.mat)), length.out = 12), ## len should be pal len + 1
#          color=RColorBrewer::brewer.pal(11, "BrBG"),
#          filename = 'wt.bcer_bgc.pdf',
#          border_color = 'white',
#          cellwidth = 15,
#          cellheight = 15
# )
bwt.mat.ss <- bwt.ss %>%
  filter(bgc %in% bwt.ss.c$bgc) %>% 
  select(bgc, Bcer.Fjoh_LFC, Bcer.Pkor_LFC, Bcer.Fjoh.Pkor_LFC) %>% 
  column_to_rownames('bgc') %>% 
  as.matrix()
# pheatmap(bwt.mat.ss,
#          filename = 'wt.bcer_bgc_ss.pdf',
#          border_color = 'white',
#          breaks=seq(-max(abs(bwt.mat.ss)), max(abs(bwt.mat.ss)), length.out = 12), ## len = pal len + 1
#          color = RColorBrewer::brewer.pal(11, "BrBG"),
#          cellwidth = 15, cellheight = 15,
#          treeheight_col = 2, treeheight_row = 12
# )

## Fjoh
fwt.mat <- fwt.ss %>% 
  select(bgc, Bcer.Fjoh_LFC, Fjoh.Pkor_LFC, Bcer.Fjoh.Pkor_LFC) %>% 
  column_to_rownames('bgc') %>% 
  as.matrix()
# pheatmap(fwt.mat,
#          filename = 'wt.fjoh_bgc.pdf',
#          border_color = 'white',
#          breaks=seq(-max(abs(fwt.mat)), max(abs(fwt.mat)), length.out = 12), ## len should be pal len + 1
#          color = RColorBrewer::brewer.pal(11, "BrBG"),
#          cellwidth = 15,
#          cellheight = 15
# )
fwt.mat.ss <- fwt.ss %>%
  filter(bgc %in% fwt.ss.c$bgc) %>% 
  select(bgc, Bcer.Fjoh_LFC, Fjoh.Pkor_LFC, Bcer.Fjoh.Pkor_LFC) %>% 
  column_to_rownames('bgc') %>% 
  as.matrix()
# pheatmap(fwt.mat.ss,
#          filename = 'wt.fjoh_bgc_ss.pdf',
#          border_color = 'white',
#          breaks=seq(-max(abs(fwt.mat.ss)), max(abs(fwt.mat.ss)), length.out = 12), ## len should be pal len + 1
#          color = RColorBrewer::brewer.pal(11, "BrBG"),
#          cellwidth = 15, cellheight = 15,
#          treeheight_col = 2, treeheight_row = 12
# )
## Pkor
kwt.mat <- kwt.ss %>% 
  select(bgc, Bcer.Pkor_LFC, Fjoh.Pkor_LFC, Bcer.Fjoh.Pkor_LFC) %>% 
  column_to_rownames('bgc') %>% 
  as.matrix()
# pheatmap(kwt.mat,
#          filename = 'wt.pkor_bgc.pdf',
#          border_color = 'white',
#          breaks=seq(-max(abs(kwt.mat)), max(abs(kwt.mat)), length.out = 12), ## len should be pal len + 1
#          color = RColorBrewer::brewer.pal(11, "BrBG"),
#          cellwidth = 15,
#          cellheight = 15
# )
k.wt.keeps <- kwt.mat %>%  
  as_tibble(rownames = 'bgc') %>% 
  group_by(bgc) %>% 
  mutate(amax=max(abs(Bcer.Pkor_LFC), abs(Fjoh.Pkor_LFC), abs(Bcer.Fjoh.Pkor_LFC))) %>% 
  arrange(-amax) %>% 
  head(n=12)
kwt.mat.ss <- kwt.ss %>%
  filter(bgc %in% k.wt.keeps$bgc) %>% 
  select(bgc, Bcer.Pkor_LFC, Fjoh.Pkor_LFC, Bcer.Fjoh.Pkor_LFC) %>% 
  column_to_rownames('bgc') %>% 
  as.matrix()
# pheatmap(kwt.mat.ss,
#          filename = 'wt.pkor_bgc_ss.pdf',
#          border_color = 'white',
#          breaks=seq(-2, 2, length.out = 12), ## len should be pal len + 1
#          color = RColorBrewer::brewer.pal(11, "BrBG"),
#          cellwidth = 15, cellheight = 15,
#          treeheight_col = 2, treeheight_row = 12
# )


#max(abs(bwt.mat))
#max(abs(fwt.mat))
#max(abs(kwt.mat))

## Make summary table
bmt <- bwt.mat %>% 
  as.data.frame() %>% 
  rownames_to_column('bgc') %>% 
  mutate(org='Bcer') %>% 
  gather(key='comparison', value='LFC', -bgc, -org) %>% 
  mutate(comparison=sub("_LFC$", "", comparison),
         coculture=case_when(
           comparison == 'Bcer.Fjoh.Pkor' ~ 'Three-member',
           comparison == 'Bcer.Fjoh' ~ 'Pairwise with Fjoh',
           comparison == 'Bcer.Pkor' ~ 'Pairwise with Pkor',
           T ~ 'ERROR'))
fmt <- fwt.mat %>% 
  as.data.frame() %>% 
  rownames_to_column('bgc') %>% 
  mutate(org='Fjoh') %>% 
  gather(key='comparison', value='LFC', -bgc, -org) %>% 
  mutate(comparison=sub("_LFC$", "", comparison),
         coculture=case_when(
           comparison == 'Bcer.Fjoh.Pkor' ~ 'Three-member',
           comparison == 'Bcer.Fjoh' ~ 'Pairwise with Bcer',
           comparison == 'Fjoh.Pkor' ~ 'Pairwise with Pkor',
           T ~ 'ERROR'))
kmt <- kwt.mat %>% 
  as.data.frame() %>% 
  rownames_to_column('bgc') %>% 
  mutate(org='Pkor') %>% 
  gather(key='comparison', value='LFC', -bgc, -org) %>% 
  mutate(comparison=sub("_LFC$", "", comparison),
         coculture=case_when(
           comparison == 'Bcer.Fjoh.Pkor' ~ 'Three-member',
           comparison == 'Bcer.Pkor' ~ 'Pairwise with Bcer',
           comparison == 'Fjoh.Pkor' ~ 'Pairwise with Fjoh',
           T ~ 'ERROR'))

wt.lfc.sum <- rbind(bmt, fmt, kmt) %>% 
  select(bgc, org, coculture, LFC) %>% 
  rename(`BGC of` = org)
write.table(wt.lfc.sum, "wt_lfc_summary.tsv", sep="\t", row.names = F, quote = F)

ggplot(wt.lfc.sum, aes(y=LFC, x=coculture, color=`BGC of`))+
  geom_hline(yintercept = 0, linetype=1, size=1, color='gray50')+
  geom_boxplot(outlier.shape = NA, alpha=.7)+
  geom_jitter(alpha=.7, position=position_jitterdodge(jitter.width=.15))+
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(face='bold'),
    legend.title = element_text(face='bold'),
    axis.title.y = element_text(face='bold'),
    legend.position = 'top'
  )+
  scale_color_manual(values=c(b.col, f.col, k.col))+
  ylab(expression(Psi*' '*log[2]*'FC'))
ggsave("wt_bgc_sum.pdf", dpi=1600, height=3.5, width=7)

pair_b <- wt.lfc.sum %>%
  filter(coculture=='Pairwise with Bcer') %>% 
  ggplot(aes(y=LFC, x=reorder(bgc, -LFC)))+
    geom_hline(yintercept = 0, linetype=1, size=.5, color='gray50')+
    geom_segment(aes(x=reorder(bgc, -LFC),
                     xend=reorder(bgc, -LFC), 
                     y=0, 
                     yend=LFC),
                 color='gray50', size=1)+
    geom_point(aes(color=`BGC of`), size=4, alpha=.9)+
    theme_bw()+
    theme(
      axis.title.y = element_blank(),
      legend.title = element_text(face='bold'),
      axis.title.x = element_text(face='bold'),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = 'none'
    )+
    scale_color_manual(values=c(f.col, k.col))+
    ylab(expression(Psi*' '*log[2]*'FC'))+
    coord_flip()
pair_f <- wt.lfc.sum %>%
  filter(coculture=='Pairwise with Fjoh') %>% 
  ggplot(aes(y=LFC, x=reorder(bgc, -LFC)))+
  geom_hline(yintercept = 0, linetype=1, size=.5, color='gray50')+
  geom_segment(aes(x=reorder(bgc, -LFC),
                   xend=reorder(bgc, -LFC), 
                   y=0, 
                   yend=LFC),
               color='gray50', size=1)+
  geom_point(aes(color=`BGC of`), size=4, alpha=.9)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(),
    legend.title = element_text(face='bold'),
    axis.title.x = element_text(face='bold'),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  )+
  scale_color_manual(values=c(b.col, k.col))+
  ylab(expression(Psi*' '*log[2]*'FC'))+
  coord_flip()
pair_k <- wt.lfc.sum %>%
  filter(coculture=='Pairwise with Pkor') %>% 
  ggplot(aes(y=LFC, x=reorder(bgc, -LFC)))+
  geom_hline(yintercept = 0, linetype=1, size=.5, color='gray50')+
  geom_segment(aes(x=reorder(bgc, -LFC),
                   xend=reorder(bgc, -LFC), 
                   y=0, 
                   yend=LFC),
               color='gray50', size=1)+
  geom_point(aes(color=`BGC of`), size=4, alpha=.9)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(),
    legend.title = element_text(face='bold'),
    axis.title.x = element_text(face='bold'),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  )+
  scale_color_manual(values=c(b.col, f.col))+
  ylab(expression(Psi*' '*log[2]*'FC'))+
  coord_flip()
three_mem <- wt.lfc.sum %>%
  filter(coculture=='Three-member') %>% 
  ggplot(aes(y=LFC, x=reorder(bgc, -LFC)))+
  geom_hline(yintercept = 0, linetype=1, size=.5, color='gray50')+
  geom_segment(aes(x=reorder(bgc, -LFC),
                   xend=reorder(bgc, -LFC), 
                   y=0, 
                   yend=LFC),
               color='gray50', size=1)+
  geom_point(aes(color=`BGC of`), size=4, alpha=.9)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(),
    legend.title = element_text(face='bold'),
    axis.title.x = element_text(face='bold'),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  )+
  scale_color_manual(values=c(b.col, f.col, k.col))+
  ylab(expression(Psi*' '*log[2]*'FC'))+
  coord_flip()
lolly_leg <- wt.lfc.sum %>%
  filter(coculture=='Three-member') %>% 
  ggplot(aes(y=LFC, x=reorder(bgc, -LFC)))+
  geom_hline(yintercept = 0, linetype=1, size=.5, color='gray50')+
  geom_segment(aes(x=reorder(bgc, -LFC),
                   xend=reorder(bgc, -LFC), 
                   y=0, 
                   yend=LFC),
               color='gray50', size=1)+
  geom_point(aes(color=`BGC of`), size=4, alpha=.9)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(),
    legend.title = element_text(face='bold'),
    axis.title.x = element_text(face='bold'),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'top'
  )+
  scale_color_manual(values=c(b.col, f.col, k.col))+
  ylab(expression(Psi*' '*log[2]*'FC'))+
  coord_flip()

plot_grid(lolly_leg,
          plot_grid(pair_b, pair_f, pair_k, three_mem, nrow=1),
          nrow=2,
          rel_heights = c(.2,1)
          )
ggsave("bgc_lollypop.pdf", dpi=1600, height = 8.5, width = 14)

## kecA-F mut
bd.cond <- c('Bcer.Pkor', 'Bcer.Pkor*', 'Bcer.Fjoh.Pkor', 'Bcer.Fjoh.Pkor*')
fd.cond <- c('Fjoh.Pkor', 'Fjoh.Pkor*', 'Bcer.Fjoh.Pkor', 'Bcer.Fjoh.Pkor*')
kd.cond <- c('Pkor', 'Pkor*',
             'Bcer.Pkor', 'Bcer.Pkor*',
             'Fjoh.Pkor', 'Fjoh.Pkor*',
             'Bcer.Fjoh.Pkor', 'Bcer.Fjoh.Pkor*')

## Create matricies for edgeR processing
bd.e <- delkec.htseq %>%
  filter(organism == 'Bcer' & condition %in% bd.cond) %>% 
  select(gene, sample, htseq_count) %>% 
  spread(key=sample, value=htseq_count) %>% 
  column_to_rownames('gene')
fd.e <- delkec.htseq %>%
  filter(organism == 'Fjoh' & condition %in% fd.cond) %>% 
  select(gene, sample, htseq_count) %>% 
  spread(key=sample, value=htseq_count) %>% 
  column_to_rownames('gene')
kd.e <- delkec.htseq %>%
  filter(organism == 'Pkor' & condition %in% kd.cond) %>% 
  select(gene, sample, htseq_count) %>% 
  spread(key=sample, value=htseq_count) %>% 
  column_to_rownames('gene')

## Calculate CPMs
bd.grp <- sub("\\d$", "", colnames(bd.e))
bd.d <- DGEList(counts=bd.e, group=factor(bd.grp))
bd.d <- calcNormFactors(bd.d) ## TMM norm
bd.cpm <- as.data.frame(cpm(bd.d)) %>% 
  rownames_to_column('gene') %>% 
  gather(key='sample', value='cpm', -gene)
fd.grp <- sub("\\d$", "", colnames(fd.e))
fd.d <- DGEList(counts=fd.e, group=factor(fd.grp))
fd.d <- calcNormFactors(fd.d) ## TMM norm
fd.cpm <- as.data.frame(cpm(fd.d)) %>% 
  rownames_to_column('gene') %>% 
  gather(key='sample', value='cpm', -gene)
kd.grp <- sub("\\d$", "", colnames(kd.e))
kd.d <- DGEList(counts=kd.e, group=factor(kd.grp))
kd.d <- calcNormFactors(kd.d) ## TMM norm
kd.cpm <- as.data.frame(cpm(kd.d)) %>% 
  rownames_to_column('gene') %>% 
  gather(key='sample', value='cpm', -gene)
d.cpm <- rbind(bd.cpm, fd.cpm, kd.cpm) %>% 
  left_join(delkec.samp2cond, by='sample') %>% 
  left_join(gene_ann, by='gene')

#############################################
## Calculate each gene's maximum per condition mean
## to define lower cutoff
d.ptile.cutoff <- 0.025
d.cond.max <- d.cpm %>% 
  group_by(gene, condition, organism) %>% 
  summarize(gmean=mean(cpm)) %>% 
  group_by(gene, organism) %>% 
  summarize(gmax=max(gmean)) %>% 
  group_by(organism) %>% 
  mutate(ptile = percent_rank(gmax))
d.line <- d.cond.max %>% 
  filter(ptile > d.ptile.cutoff) %>% 
  arrange(ptile) %>% 
  distinct(organism, .keep_all = T) %>% 
  select(organism, gmax) %>% 
  rename(cutoff=gmax)

## Plot the distributions and cutoff
bd.cp <- ggplot(d.cond.max %>% filter(organism=='Bcer'), aes(x=gmax))+
  geom_vline(xintercept = as.numeric(d.line[wt.line$organism=='Bcer','cutoff']), linetype='dashed') +
  geom_density(color='#333333', alpha=0.3, size=1, fill=b.col)+
  scale_x_log10(limits=c(1e-3, 1e6))+
  theme_bw()+theme(
    legend.title = element_blank(),
    legend.position = 'none'
  )+
  ylab('Frequency of genes')+xlab("Maximum per condition mean CPM")+
  ggtitle('Bcer')+
  scale_y_continuous(limits=c(0,.6))
fd.cp <- ggplot(d.cond.max %>% filter(organism=='Fjoh'), aes(x=gmax))+
  geom_vline(xintercept = as.numeric(d.line[wt.line$organism=='Fjoh','cutoff']), linetype='dashed') +
  geom_density(color='#333333', alpha=0.3, size=1, fill=f.col)+
  scale_x_log10(limits=c(1e-3, 1e6))+
  theme_bw()+theme(
    legend.title = element_blank(),
    legend.position = 'none'
  )+
  ylab('Frequency of genes')+xlab("Maximum per condition mean CPM")+
  ggtitle('Fjoh')+
  scale_y_continuous(limits=c(0,.6))
kd.cp <- ggplot(d.cond.max %>% filter(organism=='Pkor'), aes(x=gmax))+
  geom_vline(xintercept = as.numeric(d.line[wt.line$organism=='Pkor','cutoff']), linetype='dashed') +
  geom_density(color='#333333', alpha=0.3, size=1, fill=k.col)+
  scale_x_log10(limits=c(1e-3, 1e6))+
  theme_bw()+theme(
    legend.title = element_blank(),
    legend.position = 'none'
  )+
  ylab('Frequency of genes')+xlab("Maximum per condition mean CPM")+
  ggtitle('Pkor')+
  scale_y_continuous(limits=c(0,.6))
plot_grid(bd.cp, fd.cp, kd.cp,
          labels=LETTERS[1:3],
          nrow = 1)
ggsave("d_max_cpm.pdf", dpi=1600, height=4, width=12)

## Filter for genes above the cutoff
d.keep <- d.cond.max %>% 
  left_join(wt.line, by='organism') %>% 
  filter(gmax >= cutoff)
d.cpm.cut <- d.cpm %>% 
  filter(gene %in% d.keep$gene)
#############################################



dbgc.cpm <- d.cpm.cut %>% 
  filter(!(bgc %in% bg$replace_antismash_bgc)) %>% 
  rbind(
    d.cpm.cut %>% filter(gene %in% bg$gene)
  ) %>% 
  filter(!(is.na(bgc)))

dbgc.ss <- dbgc.cpm %>% 
  group_by(bgc, sample) %>% 
  summarize(bgc.cpm=sum(cpm)) %>% 
  left_join(
    man_bgc %>%
      select(antiSMASH_bgc, organism, size) %>% rename(bgc=antiSMASH_bgc), 
    by='bgc'
  ) %>% 
  mutate(bgc.rpkm=bgc.cpm/size*1000) %>%
  left_join(delkec.samp2cond, by='sample') %>% 
  group_by(organism, bgc, condition) %>% 
  summarize(bgc.cpm.avg=mean(bgc.cpm), bgc.cpm.sd=sd(bgc.cpm),
            bgc.rpkm.avg=mean(bgc.rpkm), bgc.rpkm.sd=sd(bgc.rpkm))

####
bd.ss <- dbgc.ss %>% 
  select(organism, bgc, condition, bgc.cpm.avg) %>% 
  filter(organism=='Bacillus cereus UW85') %>% 
  spread(key=condition, value=bgc.cpm.avg) %>% 
  ungroup() %>% 
  mutate(Bcer.Pkor_LFC = PsiLFC(`Bcer.Pkor*`, Bcer.Pkor),
         Bcer.Fjoh.Pkor_LFC = PsiLFC(`Bcer.Fjoh.Pkor*`, Bcer.Fjoh.Pkor))
fd.ss <- dbgc.ss %>% 
  select(organism, bgc, condition, bgc.cpm.avg) %>% 
  filter(organism=='Flavobacterium johnsoniae UW101') %>% 
  spread(key=condition, value=bgc.cpm.avg) %>% 
  ungroup() %>% 
  mutate(Fjoh.Pkor_LFC = PsiLFC(`Fjoh.Pkor*`, Fjoh.Pkor),
         Bcer.Fjoh.Pkor_LFC = PsiLFC(`Bcer.Fjoh.Pkor*`, Bcer.Fjoh.Pkor))
kd.ss <- dbgc.ss %>% 
  select(organism, bgc, condition, bgc.cpm.avg) %>% 
  filter(organism=='Pseudomonas koreensis CI12') %>% 
  spread(key=condition, value=bgc.cpm.avg) %>% 
  ungroup() %>% 
  mutate(Pkor_LFC = PsiLFC(`Pkor*`, Pkor),
         Bcer.Pkor_LFC = PsiLFC(`Bcer.Pkor*`, Bcer.Pkor),
         Fjoh.Pkor_LFC = PsiLFC(`Fjoh.Pkor*`, Fjoh.Pkor),
         Bcer.Fjoh.Pkor_LFC = PsiLFC(`Bcer.Fjoh.Pkor*`, Bcer.Fjoh.Pkor))


write.table(d.cpm.cut, "d_cpm.tsv", sep="\t", row.names = F, quote = F)


###########
fccut <- 2
bd.pair.cut <- bd.ss %>% filter(abs(Bcer.Pkor_LFC)>=fccut)
bd.trip.cut <- bd.ss %>% filter(abs(Bcer.Fjoh.Pkor_LFC)>=fccut)
bd.plt.ss <- dbgc.ss %>%
  filter(condition %in% c('Bcer.Pkor', 'Bcer.Pkor*', 'Bcer.Fjoh.Pkor', 'Bcer.Fjoh.Pkor*') & bgc %in% c(bd.pair.cut$bgc, bd.trip.cut$bgc)) %>% 
  left_join(man_bgc %>% select(bgc, antiSMASH_bgc) %>% rename(vname=bgc, bgc=antiSMASH_bgc), by='bgc')
bd.plt.ss$condition <- factor(bd.plt.ss$condition,
                              levels=c('Bcer.Pkor', 'Bcer.Pkor*', 'Bcer.Fjoh.Pkor', 'Bcer.Fjoh.Pkor*'))
ggplot(bd.plt.ss, aes(x=vname, y=bgc.cpm.avg, fill=condition))+
  geom_errorbar(aes(ymax=bgc.cpm.avg+bgc.cpm.sd, 
                    ymin=bgc.cpm.avg-bgc.cpm.sd),
                position=position_dodge(.9),
                width=.2)+
  geom_bar(stat='identity',
           position=position_dodge())+
  scale_fill_brewer(palette = 'Set2')+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face='bold'),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank()
  )+
  ylab("BGC CPM")
ggsave("kecmut.bcer.cut.pdf", dpi=1600, height=3, width=6)


dbgc.ss.box <- dbgc.cpm %>% 
  group_by(bgc, sample) %>% 
  summarize(bgc.cpm=sum(cpm)) %>% 
  left_join(
    man_bgc %>%
      select(antiSMASH_bgc, organism, size) %>% rename(bgc=antiSMASH_bgc), 
    by='bgc'
  ) %>% 
  mutate(bgc.rpkm=bgc.cpm/size*1000) %>%
  left_join(delkec.samp2cond, by='sample') %>% 
  left_join(man_bgc %>% select(bgc, antiSMASH_bgc) %>% rename(vname=bgc, bgc=antiSMASH_bgc), by='bgc') %>% 
  filter(vname %in% c('Petrobactin', 'Bacillibactin')) %>% 
  mutate(ccc = case_when(
    grepl('Bcer.Pkor', condition) ~ 'Pairwise',
    T ~ 'Three-member'),
    pk = case_when(
      grepl('\\*', condition) ~ 'Mutant',
      T ~ 'Wildtype'),
    pk = factor(pk, levels=c('Wildtype', 'Mutant'))
    )


ggplot(dbgc.ss.box, aes(y=bgc.cpm, x=ccc, color=pk))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=.7, position=position_jitterdodge(jitter.width=.15))+
  facet_grid(vname~ccc, scales='free') +
  theme_bw()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(face='bold'),
    axis.title.y = element_text(face='bold'),
    strip.background = element_blank(),
    strip.text.y = element_text(face='bold'),
    strip.text.x = element_blank(),
    legend.title = element_blank(),
    legend.position = 'top'
  )+
  ylab("BGC expression (CPM)")+
  scale_color_brewer(palette='Set2')
ggsave("pet_bac_mut.pdf", dpi=1600, width=5, height = 6)  
  

bp.wm.rat <- dbgc.ss.box %>% 
  group_by(vname, ccc, pk) %>% 
  summarize(mean=mean(bgc.cpm)) %>% 
  spread(pk, mean) %>% 
  mutate(inc = Mutant - Wildtype,
         pctinc = inc/Wildtype*100)

fd.pair.cut <- fd.ss %>% filter(abs(Fjoh.Pkor_LFC)>=fccut)
fd.trip.cut <- fd.ss %>% filter(abs(Bcer.Fjoh.Pkor_LFC)>=fccut)
fd.plt.ss <- dbgc.ss %>%
  filter(condition %in% c('Fjoh.Pkor', 'Fjoh.Pkor*', 'Bcer.Fjoh.Pkor', 'Bcer.Fjoh.Pkor*') & bgc %in% c(fd.pair.cut$bgc, fd.trip.cut$bgc)) %>% 
  left_join(man_bgc %>% select(bgc, antiSMASH_bgc) %>% rename(vname=bgc, bgc=antiSMASH_bgc), by='bgc')
fd.plt.ss$condition <- factor(fd.plt.ss$condition,
                              levels=c('Fjoh.Pkor', 'Fjoh.Pkor*', 'Bcer.Fjoh.Pkor', 'Bcer.Fjoh.Pkor*'))
ggplot(fd.plt.ss, aes(x=vname, y=bgc.cpm.avg, fill=condition))+
  geom_errorbar(aes(ymax=bgc.cpm.avg+bgc.cpm.sd, 
                    ymin=bgc.cpm.avg-bgc.cpm.sd),
                position=position_dodge(.9),
                width=.2)+
  geom_bar(stat='identity',
           position=position_dodge())+
  scale_fill_brewer(palette = 'Set2')+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face='bold'),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank()
  )+
  ylab("BGC CPM")
ggsave("kecmut.fjoh.cut.pdf", dpi=1600, height=3, width=8)
ggplot(fd.plt.ss %>% filter( grepl("^F", vname)), aes(x=vname, y=bgc.cpm.avg, fill=condition))+
  geom_errorbar(aes(ymax=bgc.cpm.avg+bgc.cpm.sd, 
                    ymin=bgc.cpm.avg-bgc.cpm.sd),
                position=position_dodge(.9),
                width=.2)+
  geom_bar(stat='identity',
           position=position_dodge())+
  scale_fill_brewer(palette = 'Set2')+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face='bold'),
    legend.position = 'right',
    panel.grid.major.x = element_blank()
  )+
  ylab("BGC CPM")
ggsave("kecmut.fjoh.sid.pdf", dpi=1600, height=3, width=4)

kd.sing.cut <- kd.ss %>% filter(abs(Pkor_LFC)>=fccut)
kd.pairb.cut <- kd.ss %>% filter(abs(Bcer.Pkor_LFC)>=fccut)
kd.pairf.cut <- kd.ss %>% filter(abs(Fjoh.Pkor_LFC)>=fccut)
kd.trip.cut <- kd.ss %>% filter(abs(Bcer.Fjoh.Pkor_LFC)>=fccut)
kd.plt.ss <- dbgc.ss %>%
  filter(condition %in% c('Pkor', 'Pkor*',
                          'Bcer.Pkor', 'Bcer.Pkor*',
                          'Fjoh.Pkor', 'Fjoh.Pkor*',
                          'Bcer.Fjoh.Pkor', 'Bcer.Fjoh.Pkor*') & bgc %in% c(kd.sing.cut$bgc,
                                                                            kd.pairb.cut$bgc,
                                                                            kd.pairf.cut$bgc,
                                                                            kd.trip.cut$bgc)) %>% 
  left_join(man_bgc %>% select(bgc, antiSMASH_bgc) %>% rename(vname=bgc, bgc=antiSMASH_bgc), by='bgc')
kd.plt.ss$condition <- factor(kd.plt.ss$condition,
                              levels=c('Pkor', 'Pkor*',
                                       'Bcer.Pkor', 'Bcer.Pkor*',
                                       'Fjoh.Pkor', 'Fjoh.Pkor*',
                                       'Bcer.Fjoh.Pkor', 'Bcer.Fjoh.Pkor*'))
ggplot(kd.plt.ss, aes(x=vname, y=bgc.cpm.avg, fill=condition))+
  geom_errorbar(aes(ymax=bgc.cpm.avg+bgc.cpm.sd, 
                    ymin=bgc.cpm.avg-bgc.cpm.sd),
                position=position_dodge(.9),
                width=.2)+
  geom_bar(stat='identity',
           position=position_dodge())+
  scale_fill_brewer(palette = 'Set2')+
  theme_bw()+
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face='bold'),
    legend.position = 'bottom',
    panel.grid.major.x = element_blank()
  )+
  ylab("BGC CPM")
ggsave("kecmut.pkor.cut.pdf", dpi=1600, height=3, width=6)
