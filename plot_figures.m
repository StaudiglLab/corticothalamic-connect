%% Plot Figures
% This script plots all figures found in the manuscript "Rhythmic 
% interactions between the mediodorsal thalamus and prefrontal cortex 
% precede human visual perception". All necessary data can be found within
% the subfolder "source_data".
%
% Data exceptions: 
%
%   Fig. 1a: This figure is a schematic drawing that contains no data
%
%   Fig. 1b: This figure contains data about electrode implanation
%            trajectories that would compromise the anonymity of the patient
%
%   Fig. 3a: This figure is a schematic drawing that contains no data
%
%   Supp. Fig. 1: This figure contains data about electrode implanation
%            trajectories that would compromise the anonymity of the patient
%
% Following plotting, figures were exported to Microsoft Powerpoint for
% resizing/positioning, and then to GIMP to tidy the complete figures.
% Consequently, some dimensions/colours/fonts/labels may vary.

%% Prepare Workspace
% clear workspace
clearvars
close all
clc

% define data directory
dir_repo = 'C:/Users/ra34fod/github/corticothalamic-connect/';

% add toolboxes
addpath(genpath(sprintf('%ssubfunctions',dir_repo)))

%% Figure 1.
% --- Figure 1a --- %
% This is a schematic drawing that contains no data for plotting


% --- Figure 1b --- %
% Implanation trajectory data cannot be shared


% --- Figure 1c --- %
% load data
tbl = readtable([dir_repo,'/source_data/fig_1.xlsx'],'sheet','figure_1c');

% convert longform to tfr
[vals,freqs,times] = longform_to_tfr(tbl.md_pbi,tbl.freq,tbl.time);

% plot
cmap = flipud(brewermap(11,'RdGy'));
figure('position',[100 100 500 300]); imagesc(times,freqs,vals); hold on; axis xy
caxis([-0.004 0.004]); 
xline(0,'k--'); colorbar(); colormap(cmap);
title('iEEG: Mediodorsal Thalamus');
set(gca,'ytick',5:5:20);
xlabel('Time (s)');
ylabel('Freq. (Hz)')


% --- Figure 1d --- %
% load data
tbl = readtable([dir_repo,'/source_data/fig_1.xlsx'],'sheet','figure_1d');

% plot
figure('position',[100 100 600 220]);
for pp = 1 : 6
    subplot(3,2,pp); hold on
    hits = tbl.(sprintf('vals_pp%d_hits',pp));
    misses = tbl.(sprintf('vals_pp%d_misses',pp));
    times = tbl.times;
    plot(times,misses,'color',[175 171 171]/255,'linewidth',1.5);
    plot(times,hits,'color',[208 84 70]/255,'linewidth',1.5);
    xline(0,'k--'); xlim([-0.8 0.8])
    set(gca,'tickdir','out','xtick',[-0.8 0 0.8],'ytick',[-1 1])
    title(sprintf('P%d',pp))
end


% --- Figure 1e --- %
% load data
tbl = readtable([dir_repo,'/source_data/fig_1.xlsx'],'sheet','figure_1e');

% plot
cmap = [208 84 70]./255;
figure('position',[100 100 1500 200]);
for pp = 1:6
    subplot(1,6,pp); hold on
    vals_perm = tbl.(sprintf('vals_pp%d',pp))(tbl.obs_boolean==0);
    vals_obs = tbl.(sprintf('vals_pp%d',pp))(tbl.obs_boolean==1);
    peak_freq = tbl.(sprintf('peakfreq_pp%d',pp))(tbl.obs_boolean==1);
    histogram(vals_perm,-0.003:0.0005:0.006,'normalization','probability','facecolor',cmap);
    xline(vals_obs);
    ylim([0 0.4])
    set(gca,'xtick',-0.003:0.003:0.006,'ytick',[0 0.2 0.4])
    title(sprintf('P%d - Peak: %dHz',pp,peak_freq))
    if pp==1; ylabel('Probability Count (%)'); end
end

% --- Figure 1f --- %
% load data
tbl = readtable([dir_repo,'/source_data/fig_1.xlsx'],'sheet','figure_1f');

% convert longform to tfr
[vals,freqs,times] = longform_to_tfr(tbl.mpfc_pbi,tbl.freq,tbl.time);

% plot
cmap = flipud(brewermap(11,'RdGy'));
figure('position',[100 100 500 300]); imagesc(times,freqs,vals); hold on; axis xy
caxis([-0.004 0.004]); 
xline(0,'k--'); colorbar(); colormap(cmap);
title('Patient MEG: PFC');
set(gca,'ytick',4:4:20);
xlabel('Time (s)');
ylabel('Freq. (Hz)')


% --- Figure 1g --- %
% load data
tbl = readtable([dir_repo,'/source_data/fig_1.xlsx'],'sheet','figure_1g');

% convert longform to tfr
[vals,freqs,times] = longform_to_tfr(tbl.control_pbi,tbl.freq,tbl.time);

% plot
cmap = flipud(brewermap(11,'RdGy'));
figure('position',[100 100 500 300]); imagesc(times,freqs,vals); hold on; axis xy
caxis([-0.001 0.001]); 
xline(0,'k--'); colorbar(); colormap(cmap);
title('Control MEG: PFC');
set(gca,'ytick',4:4:20);
xlabel('Time (s)');
ylabel('Freq. (Hz)')


%% Figure 2.
% --- Figure 2a --- %
% load data
tbl = readtable([dir_repo,'/source_data/fig_2.xlsx'],'sheet','fig_2a');

% convert longform to tfr
[vals,freqs,times] = longform_to_tfr(tbl.all_plv,tbl.freq,tbl.time);

% plot
c1 = brewermap(6,'YlOrBr');
c2 = brewermap(6,'Greys');
cmap = flipud(cat(1,flipud(c1(2:6,:)),c2));
figure('position',[100 100 500 300]); imagesc(times,freqs,vals); hold on; axis xy
caxis([-1.2 1.2]); 
xline(0,'k--'); colorbar(); colormap(cmap);
title('MD-mPFC Inter-Site Phase Clust.');
set(gca,'ytick',5:5:20);
xlabel('Time (s)');
ylabel('Freq. (Hz)')

% --- Figure 2b --- %
% This is a sourcemap plotted using MRIcron

% --- Figure 2c --- %
% load data
tbl = readtable([dir_repo,'/source_data/fig_2.xlsx'],'sheet','fig_2c');

% plot
figure;
theta = reshape(repmat(tbl.angle(1:6)',[2 1]),[],1);
rho = reshape(cat(1,tbl.length(1:6)',zeros(1,6)),[],1);
polarplot(theta,rho,'k-','color',c1(3,:)); hold on
polarplot([tbl.angle(7) tbl.angle(7)],[0 tbl.length(7)],'k-','color',c1(5,:))
thetaticks([0 90 180 270])
rticks([0 0.05 0.1 0.15 0.2 0.25])

% --- Figure 2d --- %
% load data
tbl = readtable([dir_repo,'/source_data/fig_2.xlsx'],'sheet','fig_2d');

% plot
cmap = c1(5,:);
figure('position',[100 100 1500 200]);
for pp = 1:6
    subplot(1,6,pp); hold on
    vals_perm = tbl.(sprintf('vals_pp%d',pp))(tbl.obs_boolean==0);
    vals_obs = tbl.(sprintf('vals_pp%d',pp))(tbl.obs_boolean==1);
    peak_freq = tbl.(sprintf('freq_pp%d',pp))(tbl.obs_boolean==1);
    histogram(vals_perm,0.04:0.005:0.16,'normalization','probability','facecolor',cmap);
    xline(vals_obs);
    xlim([0.04 0.16]); ylim([0 0.5])
    set(gca,'xtick',-0.04:0.04:0.16,'ytick',[0 0.25 0.5])
    title(sprintf('P%d - Peak: %dHz',pp,peak_freq))
    set(gca,'tickdir','out')
    if pp==1; ylabel('Probability Count (%)'); end
end

% --- Figure 2e --- %
% load data
tbl = readtable([dir_repo,'/source_data/fig_2.xlsx'],'sheet','fig_2e');

% plot
figure; hold on
colmap = brewermap(8,'Purples');
h(2)=shadedErrorBar(tbl.freqs(tbl.condition==2),smooth(tbl.vals_avg(tbl.condition==2),3),smooth(tbl.vals_sem(tbl.condition==2),3),'lineProps',{'color',[0.5 0.5 0.5]});
h(1)=shadedErrorBar(tbl.freqs(tbl.condition==1),smooth(tbl.vals_avg(tbl.condition==1),3),smooth(tbl.vals_sem(tbl.condition==1),3),'lineProps',{'color',colmap(end,:)});
xlim([0 100]); yline(0,'k-'); ylim([-1 2]); xlabel('freq. (Hz)'); ylabel('Phase Slope Index (z)');
set(gca,'tickdir','out'); legend([h(1).patch,h(2).patch],{'hits','misses'}); title('Pre-stimulus PSI');

% --- Figure 2f --- %
% load data
tbl = readtable([dir_repo,'/source_data/fig_2.xlsx'],'sheet','fig_2f');

% plot
cmap = [188 168 211]./255;
figure('position',[100 100 750 400]);
for pp = 1:6
    subplot(2,3,pp); hold on
    vals_perm = tbl.(sprintf('vals_pp%d',pp))(tbl.obs_boolean==0);
    vals_obs = tbl.(sprintf('vals_pp%d',pp))(tbl.obs_boolean==1);
    histogram(vals_perm,-0.04:0.004:0.036,'normalization','probability','facecolor',cmap);
    xline(vals_obs);
    ylim([0 0.4])
    set(gca,'xtick',[-0.04 0 0.04],'ytick',[0 0.2 0.4],'yticklabel',{'0','','0.4'},'tickdir','out')
    xlim([-0.04 0.04])
    title(sprintf('P%d',pp))
    if pp==1; ylabel('Probability Count (%)'); end
    if pp==6; xlabel('Phase Slope Index (arb. units)'); end
end


    % plot histogram
%     bins = -0.04:0.004:0.036;
%     %sorted_avg = sort(perm_avg);
%     figure; hold on
%     histogram(perm_avg,'BinEdges',bins,'normalization','probability');
%     xlim([-0.04 0.04]); ylim([0 0.4])
%     %xline(sorted_avg(round(numel(sorted_avg)*0.95)),'k--')
%     xline(obs_avg,'k-'); set(gca,'tickdir','out');
%     title(sprintf('sub-%02.0f - peak freq: %2.1fHz',pp,freqs(peak_idx+6)))

%% Figure 3.
% --- Figure 3a --- %
% This is a schematic drawing that contains no data for plotting

% --- Figure 3b --- %
% load data
tbl = readtable([dir_repo,'/source_data/fig_3.xlsx']);

% plot
cmap = [3 187 131]./255;
figure('position',[100 100 750 400]);
for pp = 1:6
    subplot(2,3,pp); hold on
    vals_perm = tbl.(sprintf('vals_pp%d',pp))(tbl.obs_boolean==0);
    vals_obs = tbl.(sprintf('vals_pp%d',pp))(tbl.obs_boolean==1);
    histogram(vals_perm,-0.1:0.015:0.5,'normalization','probability','facecolor',cmap);
    xline(vals_obs);
    ylim([0 0.5])
    set(gca,'xtick',[0 0.2 0.4],'ytick',[0 0.25 0.5],'yticklabel',{'0','','0.5'},'tickdir','out')
    xlim([-0.1 0.5])
    title(sprintf('P%d',pp))
    if pp==1; ylabel('Probability Count (%)'); end
    if pp==6; xlabel('Magnitude of Indirect Path (arb. units)'); end
end

%% Supp. Fig. 1
% Implanation trajectory data cannot be shared

%% Supp. Fig. 2
% MD / ANT / MD>ANT PBI plots

% load data
tbl = readtable([dir_repo,'/source_data/supp_fig_2.xlsx']);

% define variable names of interest
var_names = {'md_pbi','ant_pbi','md_vs_ant_pbi'};
plot_names = {'iEEG: Mediodorsal Thalamus','iEEG: Anterior Thalamus','iEEG: Mediodorsal > Anterior'};

% cycle through plots
cmap = flipud(brewermap(11,'RdGy'));
figure('position',[100 100 400 700]);
for i = 1 : numel(var_names)

    % convert longform to tfr
    [vals,freqs,times] = longform_to_tfr(tbl.(var_names{i}),tbl.freq,tbl.time);
    
    % plot
    subplot(3,1,i)
    imagesc(times,freqs,vals); hold on; axis xy
    caxis([-0.004 0.004]); 
    xline(0,'k--'); colorbar(); colormap(cmap);
    title(plot_names{i});
    set(gca,'ytick',5:5:20);
    xlabel('Time (s)');
    ylabel('Freq. (Hz)')
end

%% Supp. Fig. 3
% MD plot minus participant 2

% load data
tbl = readtable([dir_repo,'/source_data/supp_fig_3.xlsx']);

% convert longform to tfr
[vals,freqs,times] = longform_to_tfr(tbl.md_pbi,tbl.freq,tbl.time);

% plot
cmap = flipud(brewermap(11,'RdGy'));
figure;
imagesc(times,freqs,vals); hold on; axis xy
caxis([-0.004 0.004]); 
xline(0,'k--'); colorbar(); colormap(cmap);
title('Mediodorsal Thalamus [excl. sub-02]');
set(gca,'ytick',5:5:20);
xlabel('Time (s)');
ylabel('Freq. (Hz)')

%% Supp. Fig. 4
% thalamic spectral power tfrs

% load data
tbl = readtable([dir_repo,'/source_data/supp_fig_4.xlsx']);

% define variable names of interest
var_names = {'thal_raw','thal_diff','mpfc_raw','mpfc_diff'};
plot_names = {'Thalamus: All Trials','Thalamus: Hits > Misses','mPFC: All Trials','mPFC: Hits > Misses'};

% cycle through plots
figure('position',[100 100 800 400]);
for i = 1 : 4

    % convert longform to tfr
    [vals,freqs,times] = longform_to_tfr(tbl.(var_names{i}),tbl.freq,tbl.time);
    
    % plot
    subplot(2,2,i);
    imagesc(times,freqs,vals); axis xy
    colormap(flipud(brewermap(11,'RdBu')))
    caxis([-0.25 0.25]);
    xline(0,'k--'); xlabel('Time (s)');
    ylabel('Freq. (Hz)')
    colorbar();
    title(plot_names{i})
end

%% Supp. Fig. 5
% phase reset analysis

% load data
tbl = readtable([dir_repo,'/source_data/supp_fig_5.xlsx']);

% organise data
group_pr = zeros(6,2,2);
group_pr(:,1,1) = tbl.phaseReset_singleTrialPreStim;
group_pr(:,1,2) = tbl.phaseReset_singleTrialPostStim;
group_pr(:,2,1) = tbl.phaseReset_averagedPreStim;
group_pr(:,2,2) = tbl.phaseReset_averagedPostStim;

% plot
cmap = cat(1,[208 84 70]/255,[175 171 171]/255);
figure; hold on
for i = 1:2
    for j = 1:2
        plot(i+((j*0.2)-0.2),mean(group_pr(:,i,j)),'ko','markeredgecolor',cmap(j,:),'markerfacecolor',cmap(j,:));
        plot([i+((j*0.2)-0.2) i+((j*0.2)-0.2)],[mean(group_pr(:,i,j))-sem(group_pr(:,i,j)) mean(group_pr(:,i,j))+sem(group_pr(:,i,j))],'k-','color',cmap(j,:),'linewidth',2);
    end
end
h = [];
h(1)=plot([1 2],squeeze(mean(group_pr(:,:,1))),'color',cmap(1,:),'linewidth',2);
h(2)=plot([1.2 2.2],squeeze(mean(group_pr(:,:,2))),'color',cmap(2,:),'linewidth',2);
xlim([0.6 2.6]);
ylabel('Baseline Corrected Power (arb. units)')
set(gca,'xtick',[1.1 2.1],'xticklabel',{'Pre. Stim.','Post. Stim.'},'tickdir','out')
legend(h,{'Single Trial Power','Trial-Average Power'})


%% Supp. Fig. 6
% load data
tbl = readtable([dir_repo,'/source_data/supp_fig_6.csv']);

% plot parameters
count = 1;
xabslim = [0.001 0.004 0.001 0.004];
plot_title = {'Controls: Without ERP', 'Patients: Without ERP', 'Controls: With ERP', 'Patients: With ERP'};

% cycle through conditions
for erp = 0 : 1
    for samp = 0 : 1

        % shorten table
        tbl_red = tbl(tbl.hasERP==erp&tbl.isPatient==samp,:);
        
        % convert longform to tfr
        [vals,freqs,times] = longform_to_tfr(tbl_red.vals,tbl_red.freq,tbl_red.time);

        % plot
        cmap = flipud(brewermap(11,'RdGy'));
        figure;
        imagesc(times,freqs,vals); hold on; axis xy
        caxis([-xabslim(count) xabslim(count)]); 
        xline(0,'k--'); colorbar(); colormap(cmap);
        title(plot_title{count});
        set(gca,'ytick',5:5:20);
        xlabel('Time (s)');
        ylabel('Freq. (Hz)')
        count = count + 1;
    end
end

%% Supp. Fig. 7
% missing sensor plot

% load data
tbl = readtable([dir_repo,'/source_data/supp_fig_6.xlsx']);
load([dir_repo,'/source_data/layout.mat']);

% create fake struct
tml = [];
tml.avg = tbl.participant_count;
tml.time = 1;
tml.dimord = 'chan_time';
tml.label = tbl.MEG_Label;

% plot result
cmap = brewermap(7,'OrRd');
cfg = [];
cfg.layout = lay;
cfg.zlim = [-0.1 5.9];
cfg.style = 'fill';
cfg.gridscale = 400;
cfg.contournum = 6;
cfg.colormap = cmap;
cfg.interpolation = 'nearest';
cfg.comment = 'no';
ft_topoplotER(cfg,tml);


%% Supp. Fig. 8
% plots produced in MRIcron

%% Supp. Fig. 9
% additional ISPC plots

% --- supp. fig 8a --- %
% load data
tbl = readtable([dir_repo,'/source_data/supp_fig_8.xlsx'],'sheet','supp_fig_8a');

% convert longform to tfr
[vals,freqs,times] = longform_to_tfr(tbl.all_plv,tbl.freq,tbl.time);

% plot
c1 = brewermap(6,'YlOrBr');
c2 = brewermap(6,'Greys');
cmap = flipud(cat(1,flipud(c1(2:6,:)),c2));
figure('position',[100 100 500 300]); imagesc(times,freqs,vals); hold on; axis xy
caxis([-1.2 1.2]); 
xline(0,'k--'); colorbar(); colormap(cmap);
title('MD-mPFC Inter-Site Phase Clust.');
set(gca,'ytick',5:5:20);
xlabel('Time (s)');
ylabel('Freq. (Hz)')


% --- supp. fig 8b --- %
% load data
tbl = readtable([dir_repo,'/source_data/supp_fig_8.xlsx'],'sheet','supp_fig_8b');

% convert longform to tfr
[vals,freqs,times] = longform_to_tfr(tbl.beh_plv,tbl.freq,tbl.time);

% plot
c1 = brewermap(6,'YlOrBr');
c2 = brewermap(6,'Greys');
cmap = flipud(cat(1,flipud(c1(2:6,:)),c2));
figure('position',[100 100 500 300]); imagesc(times,freqs,vals); hold on; axis xy
caxis([-0.04 0.04]); 
xline(0,'k--'); colorbar(); colormap(cmap);
title('MD-mPFC Inter-Site Phase Clust.');
set(gca,'ytick',5:5:20);
xlabel('Time (s)');
ylabel('Freq. (Hz)')

%% Supp. Fig. 10
% plots produced in MRIcron

%% Supp. Fig. 11
% partial correlation histograms

% load data
tbl = readtable([dir_repo,'/source_data/supp_fig_10.xlsx']);

% plot
cmap = [3 187 131]./255;
figure('position',[100 100 750 400]);
for pp = 1:6
    subplot(2,3,pp); hold on
    vals_perm = tbl.(sprintf('vals_pp%d',pp))(tbl.obs_boolean==0);
    vals_obs = tbl.(sprintf('vals_pp%d',pp))(tbl.obs_boolean==1);
    histogram(vals_perm,-0.006:0.0003:0.006,'normalization','probability','facecolor',cmap);
    xline(vals_obs);
    ylim([0 0.5])
    set(gca,'xtick',-0.006:0.003:0.006,'ytick',[0 0.5])
    title(sprintf('P%d',pp))
    if pp==1; ylabel('Probability Count (%)'); end
    if pp==6; xlabel('Raw > Partial Correlation (arb. units)'); end
end

%% Supp. Fig 12
% thalamus ERPs

% load data
tbl = readtable([dir_repo,'/source_data/supp_fig_11.xlsx']);

figure('position',[100 100 800 800]); 

% plot MD
for patient = 1:6
    subplot(4,3,patient)
    plot(tbl.time,tbl.(sprintf('md_pp%d',patient)),'b-','linewidth',1);
    xlim([-0.25 1]); xline(0,'k--');
    xlabel('Time (s)'); ylabel('Amp. (a.u.)')
    set(gca,'box','off','tickdir','out','fontsize',8)
    
    if patient < 6
        subplot(4,3,patient+6)
        plot(tbl.time,tbl.(sprintf('ant_pp%d',patient)),'r-','linewidth',1);
        xlim([-0.25 1]); xline(0,'k--');
        xlabel('Time (s)'); ylabel('Amp. (a.u.)')
        set(gca,'box','off','tickdir','out','fontsize',8)
    end
end
