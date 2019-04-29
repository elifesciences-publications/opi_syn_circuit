% This script create Figure 6 panel C graphs, showing the spike probability fold change in presence of opioids DPDPDE and Damgo
load('data/data_cellAttached_experiments.mat') %load data

%Prepare data for plotting:
%Subset all experiment with DPDPE as agonist
DPDPEPap = [dataStruc(:).baseDPDPEspikeProb]'; %baseline spike probability
n_DPDPE = length(DPDPEPap); % determine number of cells tested
DPDPEPap(1:n_DPDPE,2) = [dataStruc(:).DPDPEspikeProb]'; % DPDPE spike probability
DPDPEPap(:,3) = log2(DPDPEPap(:,2)./DPDPEPap(:,1)); % log2 of baseline effect of DPDPE on spike probability
mDPDPEPap = DPDPEPap(:,2)./DPDPEPap(:,1); % baseline effect ratio of DPDPE on spike probability
meanFCDPDPE = nanmean(mDPDPEPap); % averaged effect ratio of DPDPE on spike prob.
semFCDPDPE = std(mDPDPEPap)/sqrt(n_DPDPE); % SEM

%Subset all experiment with DPDPE as agonist
DamgoPap = [dataStruc(:).baseDamgospikeProb]';
n_Damgo = length(DamgoPap); % determine number of cells tested
DamgoPap(1:n_Damgo,2) = [dataStruc(:).DamgospikeProb]'; % Damgo spike probability
DamgoPap(:,3) = log2(DamgoPap(:,2)./DamgoPap(:,1)); % log2 of baseline effect of Damgo on spike probability
mDamgoPap = DamgoPap(:,2)./DamgoPap(:,1); % baseline effect ratio of Damgo on spike probability
meanFCDamgo = nanmean(mDamgoPap); % averaged effect ratio of Damgo on spike prob.
semFCDamgo = nanstd(mDamgoPap)/sqrt(n_Damgo); % SEM


% Figure6c: plot ratio of spike probability for DAMGO and DPDPE. Values are
% log2 transformed.

% For plotting-purposes; append 'NaN' to index to make both DPDPEPap and DamgoPap
% matrixes equal in length (number of observations, note since we do not
% infer any stats from this matrix this is ok)
if n_DPDPE ~= n_Damgo
    if n_DPDPE > n_Damgo
        DamgoPap(n_Damgo+1 : n_DPDPE,:) = NaN;
    else
        DPDPEPap(n_DPDPE+1 : n_Damgo,:) = NaN;
    end
end
figure
plot([ones(size(DPDPEPap,1),1),ones(size(DPDPEPap,1),1)+1],[DamgoPap(:,3),DPDPEPap(:,3)],'o','MarkerSize',10,'Color',[0.7 0.7 0.7], 'markerfacecolor', [0.7 0.7 0.7])
hold on
posSEM = log2([meanFCDamgo,meanFCDPDPE]) - log2([meanFCDamgo,meanFCDPDPE] + [semFCDamgo,semFCDPDPE]);
negSEM = log2([meanFCDamgo,meanFCDPDPE]) - log2([meanFCDamgo,meanFCDPDPE] - [semFCDamgo,semFCDPDPE]);
errorbar([1,2],log2([meanFCDamgo,meanFCDPDPE]),negSEM,posSEM,'k-', 'linewidth', 2)
xlim([0.5 2.5])
ylim([-4 5])
text(0.6,4,strcat('meanDAMGO - SEM',num2str((meanFCDamgo)),' - ',num2str((semFCDamgo))))
text(0.6,3,strcat('meanDPDPE - SEM',num2str((meanFCDPDPE)),' - ',num2str((semFCDPDPE))))
%T-test stats on the baseline-agonist
disp('paired t-test on spike probabilities of baseline and DAMGO conditions')
[H,P,CI,STATS] = ttest(DamgoPap(:,1), DamgoPap(:,2))
disp('paired t-test on spike probabilities of baseline and DPDPE conditions')
[H,P,CI,STATS] = ttest(DPDPEPap(:,1), DPDPEPap(:,2))

%Figure 6 - figure supplement 1: Plot the AP latencies of baseline, agonist,
%and antagoinist conditions.
for exp = 1:size(dataStruc,1);
  latencyStimSpikes = dataStruc(exp).latencyClosestSpike;
  latencyStimSpikes(latencyStimSpikes == 0) = NaN;
  latencyStimSpikes(dataStruc(exp).stimTrigSpikes == 0)= NaN;
  for x = 1:6
    if ~isnan(dataStruc(exp).startEpochs(x))
      data = reshape(latencyStimSpikes(dataStruc(exp).startEpochs(x):dataStruc(exp).endEpochs(x),:),1,numel(latencyStimSpikes(dataStruc(exp).startEpochs(x):dataStruc(exp).endEpochs(x),:)));
      meanLatency(exp,x)=nanmean(data);
      minLatency(exp,x)=min(data);
    elseif isnan(dataStruc(exp).startEpochs(x))
      meanLatency(exp,x) = nan;
    end
  end
end
% Some quick plotting to check data. The
meanLatency = meanLatency.*1000;


%For plotting and stats: Export the data to a .csv 
latency_export = [];

latency_export.spike_latency_Damgo = [meanLatency(:,1);meanLatency(:,2);meanLatency(:,3)];
numExp = length(meanLatency(:,1));
latency_export.cellDamgo = repmat(1:numExp,1,3)';
latency_export.conditionDamgo = [repmat({'baseDAMGO'},1,numExp),repmat({'DAMGO'},1,numExp),repmat({'DAMGOCTAP'},1,numExp)]';
latency_export.expDamgo = repmat({'DAMGO'},numExp*3,1);

latency_export.spike_latency_DPDPE = [meanLatency(:,4);meanLatency(:,5);meanLatency(:,6)];
numExp = length(meanLatency(:,4));
latency_export.cellDPDPE = repmat(1:numExp,1,3)';
latency_export.conditionDPDPE = [repmat({'baseDPDPE'},1,numExp),repmat({'DPDPE'},1,numExp),repmat({'DPDPENaltrindole'},1,numExp)]';
latency_export.expDPDPE = repmat({'DPDPE'},numExp*3,1);

export_table = table(latency_export.spike_latency_Damgo, latency_export.conditionDamgo, latency_export.cellDamgo, latency_export.expDamgo, latency_export.spike_latency_DPDPE, latency_export.conditionDPDPE, latency_export.cellDPDPE, latency_export.expDPDPE, 'VariableNames', {'spike_latency_DAMGO','conditionDAMGO', 'cellDAMGO', 'expDAMGO', 'spike_latency_DPDPE','conditionDPDPE', 'cellDPDPE', 'expDPDPE'});
writetable(export_table, 'data/data_cellAttached_latencies.csv')

% Use the below script to perform a Skillings-Mack test with post-hoc 
% paired Wilcoxon signed-rank test. (uncomment the below in R) Also
% described in analyze_cellAttached_latencies.R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # Plot latency data and perform Test for Significance on spike probability
% # across conditions
% # Requires installation of packages: 'Skillings.Mack', 'ggplot2', 'cowplot'
% library(Skillings.Mack)
% library(ggplot2)
% library(cowplot)
% df = read.csv('data/data_cellAttached_latencies.csv')
% 
% # subset for all DAMGO experiments with a baseline value
% df_DAMGO = df[rep(!is.nan(df$spike_latency_DAMGO[grep('baseDAMGO', df$conditionDAMGO)]),3),c('spike_latency_DAMGO', 'conditionDAMGO','cellDAMGO','expDAMGO')]
% print('DAMGO AP latency Skillings-Mack analysis')
% Ski.Mack(df_DAMGO$spike_latency_DAMGO,df_DAMGO$conditionDAMGO,df_DAMGO$cellDAMGO)
% print('Paired Wilcoxon signed rank test baseline vs Damgo')
% # define here the two groups to compare
% first_group = df_DAMGO$spike_latency_DAMGO[df_DAMGO$conditionDAMGO=='baseDAMGO']
% second_group = df_DAMGO$spike_latency_DAMGO[df_DAMGO$conditionDAMGO=='DAMGO']
% # Below is required to calculate n and statistic in a sequence independent way
% # (first vs second group and second vs first group)
% wilcox_antero = wilcox.test(first_group, second_group, paired = TRUE)
% wilcox_retro = wilcox.test(second_group, first_group, paired = TRUE)
% statistic = min(c(wilcox_antero$statistic, wilcox_retro$statistic))
% diff = c(first_group-second_group)
% n = length(na.omit(diff[diff != 0]))
% print(paste('W(', n, ') = ', statistic, '. p = ', wilcox_antero$p.value))
% 
% print('Paired Wilcoxon signed rank test Damgo vs CTAP')
% # define here the two groups to compare
% first_group = df_DAMGO$spike_latency_DAMGO[df_DAMGO$conditionDAMGO=='DAMGO']
% second_group = df_DAMGO$spike_latency_DAMGO[df_DAMGO$conditionDAMGO=='DAMGOCTAP']
% # Below is required to calculate n and statistic in a sequence independent way
% # (first vs second group and second vs first group)
% wilcox_antero = wilcox.test(first_group, second_group, paired = TRUE)
% wilcox_retro = wilcox.test(second_group, first_group, paired = TRUE)
% statistic = min(c(wilcox_antero$statistic, wilcox_retro$statistic))
% diff = c(first_group-second_group)
% n = length(na.omit(diff[diff != 0]))
% print(paste('W(', n, ') = ', statistic, '. p = ', wilcox_antero$p.value))
% 
% # subset for all DPDPE experiments with a baseline value
% df_DPDPE = df[rep(!is.nan(df$spike_latency_DPDPE[grep('baseDPDPE', df$conditionDPDPE)]),3),c('spike_latency_DPDPE', 'conditionDPDPE','cellDPDPE','expDPDPE')]
% print('DPDPE AP latency Skillings-Mack analysis')
% Ski.Mack(df_DPDPE$spike_latency_DPDPE,df_DPDPE$conditionDPDPE,df_DPDPE$cellDPDPE)
% print('Paired Wilcoxon signed rank test baseline vs DPDPE')
% # define here the two groups to compare
% first_group = df_DPDPE$spike_latency_DPDPE[df_DPDPE$conditionDPDPE=='baseDPDPE']
% second_group = df_DPDPE$spike_latency_DPDPE[df_DPDPE$conditionDPDPE=='DPDPE']
% # Below is required to calculate n and statistic in a sequence independent way
% # (first vs second group and second vs first group)
% wilcox_antero = wilcox.test(first_group, second_group, paired = TRUE)
% wilcox_retro = wilcox.test(second_group, first_group, paired = TRUE)
% statistic = min(c(wilcox_antero$statistic, wilcox_retro$statistic))
% diff = c(first_group-second_group)
% n = length(na.omit(diff[diff != 0]))
% print(paste('W(', n, ') = ', statistic, '. p = ', wilcox_antero$p.value))
% 
% print('Paired Wilcoxon signed rank test DPDPE vs CTAP')
% # define here the two groups to compare
% first_group = df_DPDPE$spike_latency_DPDPE[df_DPDPE$conditionDPDPE=='DPDPE']
% second_group = df_DPDPE$spike_latency_DPDPE[df_DPDPE$conditionDPDPE=='DPDPENaltrindole']
% # Below is required to calculate n and statistic in a sequence independent way
% # (first vs second group and second vs first group)
% wilcox_antero = wilcox.test(first_group, second_group, paired = TRUE)
% wilcox_retro = wilcox.test(second_group, first_group, paired = TRUE)
% statistic = min(c(wilcox_antero$statistic, wilcox_retro$statistic))
% diff = c(first_group-second_group)
% n = length(na.omit(diff[diff != 0]))
% print(paste('W(', n, ') = ', statistic, '. p = ', wilcox_antero$p.value))
% 
% #for plotting combine condition and spike_latency columns.
% spike_latencies = c(df_DAMGO$spike_latency_DAMGO,df_DPDPE$spike_latency_DPDPE)
% conditions = factor(c(as.character(df_DAMGO$conditionDAMGO),as.character(df_DPDPE$conditionDPDPE)), levels=c('baseDAMGO', 'DAMGO','DAMGOCTAP','baseDPDPE','DPDPE','DPDPENaltrindole'))
% df_plot = data.frame('condition'=conditions, 'spike_latency'=spike_latencies)
% ggplot(data=df_plot, aes(x=condition, y=spike_latencies))+
%     geom_boxplot()+
%     coord_cartesian(ylim=c(0, 20))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # Till here copy the above to R for stats. and uncomment the above in R


% Rebuttal figure 3 analysis to calculate the half-width at half-max of each identified spike
for experiment = 1:length(dataStruc)
    % hwhm was calculated as unit: datapoints so per experiment calculate
    % correction factor to go from datapoints to ms
    time_factor = dataStruc(experiment).samplerate / 1000;
    [spikes, trials] = size(dataStruc(experiment).spikeLocs);
    %get the baseline half-width at half-max (hwhm) for the APs
    startBaseline = uint64(dataStruc(experiment).startEpochs(dataStruc(experiment).epochLoc(1)));
    endBaseline = uint64(dataStruc(experiment).endEpochs(dataStruc(experiment).epochLoc(1)));
    dataStruc(experiment).hwhm_baseline = dataStruc(experiment).spike_hwhm(~isnan(dataStruc(experiment).spike_hwhm(:,startBaseline:endBaseline))) / time_factor; % corrected by time_factor
    dataStruc(experiment).hwhm_baseline_median = median(dataStruc(experiment).hwhm_baseline);
    hwhm_baseline_dist(experiment) = dataStruc(experiment).hwhm_baseline_median;
    % same for agonists
    startAgonist = uint64(dataStruc(experiment).startEpochs(dataStruc(experiment).epochLoc(2)));
    endAgonist = uint64(dataStruc(experiment).endEpochs(dataStruc(experiment).epochLoc(2)));
    dataStruc(experiment).hwhm_agonist = dataStruc(experiment).spike_hwhm(~isnan(dataStruc(experiment).spike_hwhm(:,startAgonist:endAgonist))) / time_factor; % corrected by time_factor
    dataStruc(experiment).hwhm_agonist_median = nanmedian(dataStruc(experiment).hwhm_agonist);
    hwhm_agonist_dist(experiment) = dataStruc(experiment).hwhm_agonist_median;
end

% Rebuttal figure 3: plot the action potential widths for baseline and
% agonist (DAMGO and DPDPE combined). And perform Mann-Whitney U test.
figure
plot(repmat([1;2],1,14),[hwhm_baseline_dist;hwhm_agonist_dist], 'o-','markersize', 10, 'LineWidth',2, 'color',[0.7 0.7 0.7], 'markerfacecolor',[0.7 0.7 0.7])
hold on
errorbar([1,2], [nanmean(hwhm_baseline_dist), nanmean(hwhm_agonist_dist)], [nanstd(hwhm_baseline_dist)/sqrt(numel(hwhm_baseline_dist)), nanstd(hwhm_agonist_dist)/sqrt(numel(hwhm_agonist_dist))],'k-','LineWidth',2)
xlim([0.5, 2.5])
ylim([0, 5])
[P,H,STATs]=ranksum(hwhm_baseline_dist, hwhm_agonist_dist)