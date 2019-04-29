# Rebuttal figure 1c analysis of data and perform Test for Significance on EPSC
# amplitude across conditions
# Requires installation of packages: 'Skillings.Mack'
library(Skillings.Mack)
library(ggplot2)
df = read.csv('data/data_rebuttal_fig1.csv')

# subset for all DPDPE experiments with a baseline value
df_DPDPE = df[rep(!is.nan(df$effect_amplitude[grep('baseNaltrindole', df$conditionDPDPE)]),3),c('effect_amplitude', 'conditionDPDPE','cellDPDPE','expDPDPE')]
print('Naltrindole-DPDPE effect on baseline EPSC Skillings-Mack analysis')
Ski.Mack(df_DPDPE$effect_amplitude,df_DPDPE$conditionDPDPE,df_DPDPE$cellDPDPE)
print('Paired Wilcoxon signed rank test baseline vs Naltrindole')
# define here the two groups to compare
first_group = df_DPDPE$effect_amplitude[df_DPDPE$conditionDPDPE=='baseNaltrindole']
second_group = df_DPDPE$effect_amplitude[df_DPDPE$conditionDPDPE=='Naltrindole']
# Below is required to calculate n and statistic in a sequence independent way
# (first vs second group and second vs first group)
wilcox_antero = wilcox.test(first_group, second_group, paired = TRUE)
wilcox_retro = wilcox.test(second_group, first_group, paired = TRUE)
statistic = min(c(wilcox_antero$statistic, wilcox_retro$statistic))
diff = c(first_group-second_group)
n = length(na.omit(diff[diff != 0]))
print(paste('W(', n, ') = ', statistic, '. p = ', wilcox_antero$p.value))

print('Paired Wilcoxon signed rank test Naltrindole vs DPDPENaltrindole')
# define here the two groups to compare
first_group = df_DPDPE$effect_amplitude[df_DPDPE$conditionDPDPE=='Naltrindole']
second_group = df_DPDPE$effect_amplitude[df_DPDPE$conditionDPDPE=='DPDPENaltrindole']
# Below is required to calculate n and statistic in a sequence independent way
# (first vs second group and second vs first group)
wilcox_antero = wilcox.test(first_group, second_group, paired = TRUE)
wilcox_retro = wilcox.test(second_group, first_group, paired = TRUE)
statistic = min(c(wilcox_antero$statistic, wilcox_retro$statistic))
diff = c(first_group-second_group)
n = length(na.omit(diff[diff != 0]))
print(paste('W(', n, ') = ', statistic, '. p = ', wilcox_antero$p.value))

# Rebuttal figure 2 Analysis of data and perform plotting for DAMGO effect on
# thalamostriatal projections, to check whether there is a potential
# heterogeneity in the data that would suggest that inputs to D1 and D2 MSNs in
# the striatum may be modulated differently.
#import packages
library(ggplot2)
library(reshape2)
library(cowplot)
require(cowplot)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(heatmap.plus) # clustering visualization
library(colorspace) # colors
library(RColorBrewer)
library(colorRamps)
library(gplots)
library(plyr)
library(Skillings.Mack)

df=read.csv('data/effectdataset.csv')
df$X = NULL
msn_df=subset(df, cellType=='MSN' & drugGroup==2)
setHook(packageEvent("grDevices", "onLoad"),
function(...) grDevices::X11.options(type='cairo'))
options(device='x11')
# subset data for thalamostriatal projections only (including MD and AM)
thal_df = subset(msn_df, (stimSource=='MD' | stimSource=='AM'))

# Rebuttal figure 2a, plot distribution of DAMGO effect on thalamostriatal
# projections.

ggplot(data=subset(thal_df, parameter=='amplitude'), aes(agonistEffect)) +
    geom_histogram(binwidth=25/4)+
    ylim(0,10)+
    xlim(0,100)

print('damgoEffect normality test')
shapiro.test(subset(thal_df, parameter=='amplitude')$agonistEffect)

# First re-organize the df to allow for subsetting
sthal_df = thal_df[,c(1:12, 20, 15)]
tthal_df <- subset(sthal_df, parameter=='amplitude')
colnames(tthal_df)[colnames(tthal_df)=='baseValue'] <- 'amplitude'
colnames(tthal_df)[colnames(tthal_df)=='agonistEffect'] <- 'damgoEffect'
tthal_df$parameter <- NULL
tthal_df$onset <- subset(sthal_df, parameter=='onset')$baseValue
tthal_df$risetime <- subset(sthal_df, parameter=='risetime')$baseValue
tthal_df$slope <- subset(sthal_df, parameter=='slope')$baseValue
tthal_df$onset2peakTime <- subset(sthal_df, parameter=='onset2peakTime')$baseValue
tthal_df$halfwidth <- subset(sthal_df, parameter=='halfwidth')$baseValue
tthal_df$decay <- subset(sthal_df, parameter=='decay')$baseValue
tthal_df$area <- subset(sthal_df, parameter=='area')$baseValue
tthal_df$chargeOAmplitude <- subset(sthal_df, parameter=='chargeOAmplitude')$baseValue
# fill out NA for onset of will20180705_732c002 with median of onset
tthal_df[tthal_df$cellID == 'will20180705_732c002',]$onset = median(tthal_df$onset, na.rm=TRUE)
rownames(tthal_df) = 1:nrow(tthal_df)
tthal_df_ref = tthal_df
tthal_df = tthal_df[,-12]

# Rebuttal figure 2b, plot distribution of baseline EPSC risetime vs DAMGO
# effect for thalamostriatal inputs. Perform linear regression to check for
# correlation between effect of DAMGO and risetime of the baseline EPSC.
ggplot(data=tthal_df_ref, aes(risetime, damgoEffect))+
    geom_point(size=5)+
    geom_smooth(method=lm, se=TRUE)+
    ylim(0,100)+
    xlim(0,3)
model = lm(risetime ~ damgoEffect, data=tthal_df_ref)
summary(model)

# Rebuttal figure 2c-d, define two subpopulations based on the baseline EPSC
# parameters. Based on those two subpopulations, test whether DAMGO effect on
# baseline is significantly different between the two populations.

#Generate a matrix with all parameters and make sure that no constant or zero
# columns are present
input_thal = tthal_df[12:ncol(tthal_df)]
which(apply(input_thal, 2, var)==0) # some variables have a constant or zero column, those need to be removed
input_thal = input_thal[, apply(input_thal, 2, var)!=0]
input_thal = scale(input_thal, center=TRUE, scale=TRUE) # center to the mean value and divide by standard deviation
# rownames(input_thal) = 1:nrow(input_thal)

# generate distance matrixes for the cells and parameters:
d_HCA_thal_param = get_dist(t(input_thal), method="pearson")
d_HCA_thal_cell = get_dist((input_thal), method="euclidean")
# Check that for both matrixes which clustering methods produces the strongest
# clusters
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
ac <- function(x) {
  agnes(d_HCA_thal_cell, method = x)$ac
}
map_dbl(m, ac)
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
ac <- function(x) {
  agnes(d_HCA_thal_param, method = x)$ac
}
map_dbl(m, ac)

# Based on the comparison of methods it turns out that Ward provides the strongest
# clusters for both cell and parameter clustering. So continue with Ward.
HCA_thal_cell = agnes(d_HCA_thal_cell, method = "ward")
HCA_thal_param = agnes(d_HCA_thal_param, method = "ward")
heatmap(t(input_thal), Colv=as.dendrogram(HCA_thal_cell), Rowv=as.dendrogram(HCA_thal_param), scale='none', labRow = TRUE, labCol = TRUE, xlab="cells", ylab="parameters")
HCA_thal_cell$order[1:15]
HCA_thal_param$order.lab

# Plot and test whether the two putative subpopulations have a difference in DAMGO effect on Baseline
tthal_df_ref$cluster = rep(1,nrow(tthal_df_ref))
tthal_df_ref[clust1,]$cluster = 1
tthal_df_ref[clust2,]$cluster = 2
ggplot(tthal_df_ref, aes(damgoEffect)) +
    geom_histogram(data=subset(tthal_df_ref,cluster==1), alpha = 0.5, fill='red', binwidth=25/4)+
    geom_histogram(data=subset(tthal_df_ref,cluster==2), alpha = 0.5, fill='blue',  binwidth=25/4)+
    ylim(0,10)
print('Warning, if alpha in the above plot is not allowed, this will stall the plotting and subsequent stats')
print('please remove the alpha=0.5 from the above ggplot section')
clust1 = HCA_thal_cell$order[1:16]
clust2 = HCA_thal_cell$order[17:nrow(input_thal)]
thal_clust1 = tthal_df_ref$damgoEffect[clust1]
thal_clust2 = tthal_df_ref$damgoEffect[clust2]
wilcox.test(x=thal_clust1,
            y=thal_clust2,
            paired=FALSE)
