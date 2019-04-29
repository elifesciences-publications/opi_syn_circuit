% Read data
load('data/data_rebuttal_experiments.mat');

% Rebuttal figure 1c
% calculate percentage affect Naltrindole, DPDPE relative to baseline
plot([zeros(1,rec_cells)+0.5;ones(1,rec_cells)],[[df.nal_norm];[df.dpdpe_norm]]*100,'-o', 'color',[0.7 0.7 0.7])
hold on
errorbar([0.5,1],[nanmean([df.nal_norm]),nanmean([df.dpdpe_norm])]*100, [nanstd([df.nal_norm]*100)/sqrt(7),nanstd([df.dpdpe_norm]*100)/sqrt(7)],'color',[0 0 0])
xlim([0,1.5])
ylim([0, 2])
xticks([0.5 1])
xticklabels(['Nal','dpdpe'])
ylabel('% baseline EPSC')

%For stats: Export the data to a .csv
rebuttalfig1_export = [];

rebuttalfig1_export.effect_amplitude = [[df.base_amplitude],[df.nal_amplitude], [df.dpdpe_amplitude]];
numExp = length([df.base_amplitude]);
rebuttalfig1_export.cellDPDPE = repmat({df.cellID},1,3)';
rebuttalfig1_export.conditionDPDPE = [repmat({'baseNaltrindole'},1,numExp),repmat({'Naltrindole'},1,numExp),repmat({'DPDPENaltrindole'},1,numExp)]';
rebuttalfig1_export.expDPDPE = repmat({'DPDPE'},numExp*3,1);

export_table = table(rebuttalfig1_export.effect_amplitude', rebuttalfig1_export.conditionDPDPE, rebuttalfig1_export.cellDPDPE, rebuttalfig1_export.expDPDPE, 'VariableNames', {'effect_amplitude','conditionDPDPE', 'cellDPDPE', 'expDPDPE'});
writetable(export_table, 'data/data_rebuttal_fig1.csv')

% Use the below script to perform a Skillings-Mack test with post-hoc
% paired Wilcoxon signed-rank test. (uncomment the below in R) Also
% described in analyze_rebuttal_experiments.r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # Rebuttal figure 1c analysis of data and perform Test for Significance on EPSC
% # amplitude across conditions
% # Requires installation of packages: 'Skillings.Mack'
% library(Skillings.Mack)
% library(ggplot2)
% df = read.csv('data/data_rebuttal_fig1.csv')
%
% # subset for all DPDPE experiments with a baseline value
% df_DPDPE = df[rep(!is.nan(df$effect_amplitude[grep('baseNaltrindole', df$conditionDPDPE)]),3),c('effect_amplitude', 'conditionDPDPE','cellDPDPE','expDPDPE')]
% print('Naltrindole-DPDPE effect on baseline EPSC Skillings-Mack analysis')
% Ski.Mack(df_DPDPE$effect_amplitude,df_DPDPE$conditionDPDPE,df_DPDPE$cellDPDPE)
% print('Paired Wilcoxon signed rank test baseline vs Naltrindole')
% # define here the two groups to compare
% first_group = df_DPDPE$effect_amplitude[df_DPDPE$conditionDPDPE=='baseNaltrindole']
% second_group = df_DPDPE$effect_amplitude[df_DPDPE$conditionDPDPE=='Naltrindole']
% # Below is required to calculate n and statistic in a sequence independent way
% # (first vs second group and second vs first group)
% wilcox_antero = wilcox.test(first_group, second_group, paired = TRUE)
% wilcox_retro = wilcox.test(second_group, first_group, paired = TRUE)
% statistic = min(c(wilcox_antero$statistic, wilcox_retro$statistic))
% diff = c(first_group-second_group)
% n = length(na.omit(diff[diff != 0]))
% print(paste('W(', n, ') = ', statistic, '. p = ', wilcox_antero$p.value))
%
% print('Paired Wilcoxon signed rank test Naltrindole vs DPDPENaltrindole')
% # define here the two groups to compare
% first_group = df_DPDPE$effect_amplitude[df_DPDPE$conditionDPDPE=='Naltrindole']
% second_group = df_DPDPE$effect_amplitude[df_DPDPE$conditionDPDPE=='DPDPENaltrindole']
% # Below is required to calculate n and statistic in a sequence independent way
% # (first vs second group and second vs first group)
% wilcox_antero = wilcox.test(first_group, second_group, paired = TRUE)
% wilcox_retro = wilcox.test(second_group, first_group, paired = TRUE)
% statistic = min(c(wilcox_antero$statistic, wilcox_retro$statistic))
% diff = c(first_group-second_group)
% n = length(na.omit(diff[diff != 0]))
% print(paste('W(', n, ') = ', statistic, '. p = ', wilcox_antero$p.value))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # Till here copy the above to R for stats. and uncomment the above in R


% Subset data based on type of cell:  corticostriatal projecting cells are 1
rec_cells = length(df);
df_sub = [1:rec_cells;df.corticostriatal;df.EPSC_charge];
df_sub_neg = df_sub([1,3], df_sub(2, :) == 0);
df_sub_pos = df_sub([1,3], df_sub(2, :) == 1);
[temp, ind_pos] = sort(df_sub_pos(2,:),2,'ascend');
[temp, ind_neg] = sort(df_sub_neg(2,:),2,'ascend');

% Rebuttal figure 4d
figure
boxplot([df(df_sub_pos(1,:)).EPSC_charge], 'orientation','vertical');
ylim([0,5]);

% Rebuttal figure 4e
figure
boxplot([df(df_sub_pos(1,:)).EPSC_onset], 'orientation','vertical');
ylim([0,5]);

% BONUS:Check whether there is any correlation between onset and amplitude
% or charge. Plot summary graphs for EPSC amplitude and charge for
% corticostriatal neurons

% amplitude vs onset
figure(2)
subplot(2,2,1);
boxplot([df(df_sub_pos(1,:)).EPSC_amplitude]*-1, 'orientation','vertical');
figure(2)
ylim([0,500]);ylabel('amplitude (pA)');
figure(2)
subplot(2,2,2);
plot([df(df_sub_pos(1,:)).EPSC_onset], [df(df_sub_pos(1,:)).EPSC_amplitude]*-1,'ok');
xlim([0,5]);
ylim([0,500]);
figure(2)
subplot(2,2,4);
boxplot([df(df_sub_pos(1,:)).EPSC_onset], 'orientation','horizontal');
xlim([0,5]);
xlabel('onset (ms)');

% charge vs onset
figure(3)
subplot(2,2,1);
boxplot([df(df_sub_pos(1,:)).EPSC_charge], 'orientation','vertical');
ylim([0,5]);
ylabel('charge (pC)');
figure(3)
subplot(2,2,2);
plot([df(df_sub_pos(1,:)).EPSC_onset], [df(df_sub_pos(1,:)).EPSC_charge],'ok');
ylim([0,5]);
xlim([0,5]);
figure(3)
subplot(2,2,4);
boxplot([df(df_sub_pos(1,:)).EPSC_onset], 'orientation','horizontal');
xlim([0,5]);
xlabel('onset (ms)');
