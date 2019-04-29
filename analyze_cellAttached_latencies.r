# Plot latency data and perform Test for Significance on spike probability
# across conditions
# Requires installation of packages: 'Skillings.Mack', 'ggplot2', 'cowplot'
library(Skillings.Mack)
library(ggplot2)
library(cowplot)
df = read.csv('data/data_cellAttached_latencies.csv')

# subset for all DAMGO experiments with a baseline value
df_DAMGO = df[rep(!is.nan(df$spike_latency_DAMGO[grep('baseDAMGO', df$conditionDAMGO)]),3),c('spike_latency_DAMGO', 'conditionDAMGO','cellDAMGO','expDAMGO')]
print('DAMGO AP latency Skillings-Mack analysis')
Ski.Mack(df_DAMGO$spike_latency_DAMGO,df_DAMGO$conditionDAMGO,df_DAMGO$cellDAMGO)
print('Paired Wilcoxon signed rank test baseline vs Damgo')
# define here the two groups to compare
first_group = df_DAMGO$spike_latency_DAMGO[df_DAMGO$conditionDAMGO=='baseDAMGO']
second_group = df_DAMGO$spike_latency_DAMGO[df_DAMGO$conditionDAMGO=='DAMGO']
# Below is required to calculate n and statistic in a sequence independent way
# (first vs second group and second vs first group)
wilcox_antero = wilcox.test(first_group, second_group, paired = TRUE)
wilcox_retro = wilcox.test(second_group, first_group, paired = TRUE)
statistic = min(c(wilcox_antero$statistic, wilcox_retro$statistic))
diff = c(first_group-second_group)
n = length(na.omit(diff[diff != 0]))
print(paste('W(', n, ') = ', statistic, '. p = ', wilcox_antero$p.value))

print('Paired Wilcoxon signed rank test Damgo vs CTAP')
# define here the two groups to compare
first_group = df_DAMGO$spike_latency_DAMGO[df_DAMGO$conditionDAMGO=='DAMGO']
second_group = df_DAMGO$spike_latency_DAMGO[df_DAMGO$conditionDAMGO=='DAMGOCTAP']
# Below is required to calculate n and statistic in a sequence independent way
# (first vs second group and second vs first group)
wilcox_antero = wilcox.test(first_group, second_group, paired = TRUE)
wilcox_retro = wilcox.test(second_group, first_group, paired = TRUE)
statistic = min(c(wilcox_antero$statistic, wilcox_retro$statistic))
diff = c(first_group-second_group)
n = length(na.omit(diff[diff != 0]))
print(paste('W(', n, ') = ', statistic, '. p = ', wilcox_antero$p.value))

# subset for all DPDPE experiments with a baseline value
df_DPDPE = df[rep(!is.nan(df$spike_latency_DPDPE[grep('baseDPDPE', df$conditionDPDPE)]),3),c('spike_latency_DPDPE', 'conditionDPDPE','cellDPDPE','expDPDPE')]
print('DPDPE AP latency Skillings-Mack analysis')
Ski.Mack(df_DPDPE$spike_latency_DPDPE,df_DPDPE$conditionDPDPE,df_DPDPE$cellDPDPE)
print('Paired Wilcoxon signed rank test baseline vs DPDPE')
# define here the two groups to compare
first_group = df_DPDPE$spike_latency_DPDPE[df_DPDPE$conditionDPDPE=='baseDPDPE']
second_group = df_DPDPE$spike_latency_DPDPE[df_DPDPE$conditionDPDPE=='DPDPE']
# Below is required to calculate n and statistic in a sequence independent way
# (first vs second group and second vs first group)
wilcox_antero = wilcox.test(first_group, second_group, paired = TRUE)
wilcox_retro = wilcox.test(second_group, first_group, paired = TRUE)
statistic = min(c(wilcox_antero$statistic, wilcox_retro$statistic))
diff = c(first_group-second_group)
n = length(na.omit(diff[diff != 0]))
print(paste('W(', n, ') = ', statistic, '. p = ', wilcox_antero$p.value))

print('Paired Wilcoxon signed rank test DPDPE vs Naltrindole')
# define here the two groups to compare
first_group = df_DPDPE$spike_latency_DPDPE[df_DPDPE$conditionDPDPE=='DPDPE']
second_group = df_DPDPE$spike_latency_DPDPE[df_DPDPE$conditionDPDPE=='DPDPENaltrindole']
# Below is required to calculate n and statistic in a sequence independent way
# (first vs second group and second vs first group)
wilcox_antero = wilcox.test(first_group, second_group, paired = TRUE)
wilcox_retro = wilcox.test(second_group, first_group, paired = TRUE)
statistic = min(c(wilcox_antero$statistic, wilcox_retro$statistic))
diff = c(first_group-second_group)
n = length(na.omit(diff[diff != 0]))
print(paste('W(', n, ') = ', statistic, '. p = ', wilcox_antero$p.value))

#for plotting combine condition and spike_latency columns.
spike_latencies = c(df_DAMGO$spike_latency_DAMGO,df_DPDPE$spike_latency_DPDPE)
conditions = factor(c(as.character(df_DAMGO$conditionDAMGO),as.character(df_DPDPE$conditionDPDPE)), levels=c('baseDAMGO', 'DAMGO','DAMGOCTAP','baseDPDPE','DPDPE','DPDPENaltrindole'))
df_plot = data.frame('condition'=conditions, 'spike_latency'=spike_latencies)
ggplot(data=df_plot, aes(x=condition, y=spike_latencies))+
    geom_boxplot()+
    coord_cartesian(ylim=c(0, 20))
