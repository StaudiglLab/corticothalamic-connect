%% Phase Reset
% clear workspace
clearvars
close all
clc

% define data directory
dir_data = 'F:/meg-ieeg_data_v2/';
dir_repo = 'C:/Users/ra34fod/github/corticothalamic-connect/';

% add toolboxes
addpath(genpath(sprintf('%ssubfunctions',dir_repo)))

% define participants
npp = 6;

%% Phase Reset: All Trials
% cycle through participants
for pp = 1 : 6
                
    % define dataset
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,pp);
    filename = sprintf('sub-%02.0f_ephys-bipolar.mat',pp);
    
    % load data
    load([directory,filename],'ephys')
    
    % record trialinfo
    trialinfo = ephys.trialinfo(:,3);
    
    % get averaged signal
    ephys_avg = ephys;
    ephys_avg.time = ephys_avg.time(1);
    signal = reshape(cell2mat(ephys_avg.trial),[numel(ephys.label) numel(ephys.time{1}) numel(ephys.time)]);
    ephys_avg.trial = {mean(signal,3)};
    
    % get phase of ROI channels
    cfg         = [];
    cfg.output  = 'pow';
    cfg.method  = 'wavelet';
    cfg.width   = 6;
    cfg.toi     = -0.2:0.025:0.2;
    cfg.foi     = 6:9;
    cfg.pad     = 'nextpow2';
    cfg.channel = 'MD';
    freq        = ft_freqanalysis(cfg,ephys); clear ephys
    freq_avg    = ft_freqanalysis(cfg,ephys_avg); clear ephys_avg
    
    % get mean pow
    mean_pow = (mean(mean(freq.powspctrm(:)))+mean(mean(freq_avg.powspctrm(:)))) / 2;
    
    % extract metric
    group_pr(pp,1,1) = mean(mean(freq.powspctrm(:,:,freq.time<0)))-mean_pow;
    group_pr(pp,1,2) = mean(mean(freq.powspctrm(:,:,freq.time>0)))-mean_pow;
    group_pr(pp,2,1) = mean(mean(freq_avg.powspctrm(:,:,freq_avg.time<0)))-mean_pow;
    group_pr(pp,2,2) = mean(mean(freq_avg.powspctrm(:,:,freq_avg.time>0)))-mean_pow;
end

% write xls data
tbl = table(group_pr(:,1,1),group_pr(:,1,2),group_pr(:,2,1),group_pr(:,2,1),...
    'VariableNames',{'phaseReset_singleTrialPreStim','phaseReset_singleTrialPostStim',...
                     'phaseReset_averagedPreStim','phaseReset_averagedPostStim'});
writetable(tbl,[dir_repo,'/source_data/supp_fig_5.xlsx'],'sheet','supp_fig_5')

    