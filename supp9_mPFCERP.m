%% Prepare Workspace
% this script runs additional analyses required for the plots of figure 1

% clear workspace
clearvars
close all
clc

% define data directory
dir_data = 'F:/meg-ieeg_magdeburg/';
dir_repo = 'C:/Users/ra34fod/github/corticothalamic-connect/';

% add toolboxes
addpath(genpath(sprintf('%ssubfunctions',dir_repo)))

% define participants
npp = 6;

%% Estimate ERP
% get sourcemodel
[atlas,sourcemodel] = get_sourcemodel('mpfc');

% cycle through participants
group_pbi = []; trl_count = []; group_erp = [];
for pp = 1 : 6
                
    % define dataset
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,pp);
    filename = sprintf('sub-%02.0f_ephys-source.mat',pp);
    
    % load data
    load([directory,filename],'ephys')

    % use all channels
    label = ephys.label(atlas.mpfc_mask(sourcemodel.inside)>1);

    % select trials
    ephys = ft_selectdata(struct('channel',{label}),ephys);
    
    % get trialinfo
    trialinfo = ephys.trialinfo(:,3);
      
    % subtract erp
    tml = ft_timelockanalysis([],ephys);
    group_erp(pp,:,1) = mean(tml.avg);
    tml = ft_timelockanalysis(struct('trials',find(trialinfo==1)),ephys);
    tml = ft_timelockbaseline(struct('baseline',[-0.8 0]),tml);
    group_erp(pp,:,2) = mean(tml.avg);
    tml = ft_timelockanalysis(struct('trials',find(trialinfo==0)),ephys);
    tml = ft_timelockbaseline(struct('baseline',[-0.8 0]),tml);
    group_erp(pp,:,3) = mean(tml.avg);
end

% cycle through participants
hlth_pbi = []; hlth_count = []; hlth_erp = [];
for hlth = 1 : 12
                    
    % load data
    directory = sprintf('%s/derivatives/hlth-%02.0f/ephys/',dir_data,hlth);
    filename = sprintf('%s/hlth-%02.0f_ephys-sourceMAG_5Hz.mat',directory,hlth);
    load(filename,'ephys')

    % use all channels
    label = ephys.label(atlas.mpfc_mask(sourcemodel.inside)>1);

    % select trials
    ephys = ft_selectdata(struct('channel',{label}),ephys);
    
    % filter
    %ephys = ft_preprocessing(struct('hpfilter','yes','hpfreq',5),ephys);
    
    % get trialinfo
    trialinfo = ephys.trialinfo(:,4);
    hlth_count(hlth,:) = [sum(trialinfo==1) sum(trialinfo==0)];
     
    % subtract erp
    tml = ft_timelockanalysis([],ephys);
    hlth_erp(hlth,:,1) = mean(tml.avg);
    tml = ft_timelockanalysis(struct('trials',find(trialinfo==1)),ephys);
    tml = ft_timelockbaseline(struct('baseline',[-0.8 0]),tml);
    hlth_erp(hlth,:,2) = mean(tml.avg);
    tml = ft_timelockanalysis(struct('trials',find(trialinfo==0)),ephys);
    tml = ft_timelockbaseline(struct('baseline',[-0.8 0]),tml);
    hlth_erp(hlth,:,3) = mean(tml.avg);
end

% update user
fprintf('\ndone...\n\n')

% define data as fieldtrip structure
tml1 = struct('cfg',[],'dimord','subj_chan_time','time',linspace(-2,2,size(hlth_erp,2)),...
                'label',{{'dummy'}},'avg',permute(hlth_erp(:,:,2),[1 3 2]));
tml2 = struct('cfg',[],'dimord','subj_chan_time','time',linspace(-2,2,size(hlth_erp,2)),...
                'label',{{'dummy'}},'avg',permute(hlth_erp(:,:,3),[1 3 2]));
tml3 = struct('cfg',[],'dimord','subj_chan_time','time',linspace(-2,2,size(group_erp,2)),...
                'label',{{'dummy'}},'avg',permute(group_erp(:,:,2),[1 3 2]));
tml4 = struct('cfg',[],'dimord','subj_chan_time','time',linspace(-2,2,size(group_erp,2)),...
                'label',{{'dummy'}},'avg',permute(group_erp(:,:,3),[1 3 2]));
           
% set random seed
rng(1)

% define general stat config structure
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 2000;
cfg.avgoverchan         = 'yes';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.parameter           = 'avg';
cfg.design              = zeros(2,size(tml1.avg,1)*2);
cfg.design(1,:)         = repmat(1:size(tml1.avg,1),[1 2]);
cfg.design(2,:)         = [ones(1,size(tml1.avg,1)),ones(1,size(tml1.avg,1))+1];
cfg.latency             = [0 0.8];
cfg.statistic           = 'ft_statfun_depsamplesT';  
stat_hlth               = ft_timelockstatistics(cfg,tml1,tml2);
cfg.design              = zeros(2,size(tml3.avg,1)*2);
cfg.design(1,:)         = repmat(1:size(tml3.avg,1),[1 2]);
cfg.design(2,:)         = [ones(1,size(tml3.avg,1)),ones(1,size(tml3.avg,1))+1];
stat_pat                = ft_timelockstatistics(cfg,tml3,tml4);

% get bayes at peak
addpath('C:\Users\ra34fod\github\bayesFactor\')
pidx_hlth = find(max(stat_hlth.stat(:))==stat_hlth.stat(:));
stat_hlth.bayes = bf.ttest(squeeze(tml2.avg(:,pidx_hlth)-tml1.avg(:,pidx_hlth)));
pidx_pat = find(max(stat_pat.stat(:))==stat_pat.stat(:));
stat_pat.bayes = bf.ttest(squeeze(tml4.avg(:,pidx_pat)-tml3.avg(:,pidx_pat)));

% get effect size at peak
stat_hlth.dz = abs(stat_hlth.stat ./ sqrt(size(tml1.avg,1)));
stat_pat.dz = abs(stat_pat.stat ./ sqrt(size(tml3.avg,1)));

% report result
fprintf('\n--- Statistics: Patients ---\nt(%d): %3.3f\np: %3.3f\nBF: %3.3f\nd: %3.3f\n\n',size(tml3.avg,1)-1,stat_pat.stat(pidx_pat),stat_pat.negclusters(1).prob,stat_pat.bayes,stat_pat.dz(pidx_pat))
fprintf('\n--- Statistics: Healthy ---\nt(%d): %3.3f\np: %3.3f\nBF: %3.3f\nd: %3.3f\n\n',size(tml1.avg,1)-1,stat_hlth.stat(pidx_hlth),stat_hlth.negclusters(1).prob,stat_hlth.bayes,stat_hlth.dz(pidx_hlth))

%% Plot 
cmap = flipud(brewermap(11,'RdGy'));
            
figure;
plot(linspace(-2,2,size(hlth_erp,2)),smooth(mean(hlth_erp(:,:,3),1)),'color',cmap(4,:),'linewidth',2); hold on
plot(linspace(-2,2,size(hlth_erp,2)),smooth(mean(hlth_erp(:,:,2),1)),'color',cmap(9,:),'linewidth',2); hold on
xlim([-0.8 0.8])
xlabel('time (s)')
ylabel('amp. (uV)')
xline(0,'k--')
yline(0,'k-')
set(gca,'box','off','tickdir','out')


figure;
plot(linspace(-2,2,size(hlth_erp,2)),smooth(mean(group_erp(:,:,3),1)),'color',cmap(4,:),'linewidth',2); hold on
plot(linspace(-2,2,size(hlth_erp,2)),smooth(mean(group_erp(:,:,2),1)),'color',cmap(9,:),'linewidth',2); hold on
xlim([-0.8 0.8])
xlabel('time (s)')
ylabel('amp. (uV)')
xline(0,'k--')
yline(0,'k-')
set(gca,'box','off','tickdir','out')
