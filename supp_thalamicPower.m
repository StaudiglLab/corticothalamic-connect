%% Prepare Workspace
% clear workspace
clearvars
close all
clc

% define data directory
dir_data = 'C:/Users/ra34fod/Documents/meg-ieeg_data_v2/';
dir_repo = 'C:/Users/ra34fod/github/corticothalamic-connect/';

% add toolboxes
addpath(genpath(sprintf('%ssubfunctions',dir_repo)))

% define participants
npp = 6;

%% Compute Phase Bifurcation  
% predefine group output
group_freq = cell(6,1);

% cycle through files
for patient = 1 : 6
                               
    % define patient directory
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,patient);
    filename = sprintf('%s/sub-%02.0f_ephys-bipolar.mat',directory,patient);
    
    % load data
    load(filename,'ephys')

    % select MD
    ephys = ft_selectdata(struct('channel','MD'),ephys);
    
    % get phase using variable width wavelets
    [~,trialinfo,label,freqs,times,power] = get_tfr(ephys,ephys.label);

    % get hit > miss power difference
    power = permute(mean(power(trialinfo==1,:,:,:)) - mean(power(trialinfo==0,:,:,:)),[2 3 4 1]);
    
    % create frequency structure
    group_freq{patient} = struct('label',{{'MD'}},'freq',freqs,'time',times,...
                  'dimord','chan_freq_time','cfg',[],'powspctrm',power);
end

%% Run Statistics
% get grand average of participants
cfg                 = [];
cfg.keepindividual  = 'yes';
freq                = ft_freqgrandaverage(cfg,group_freq{:});

% set random seed
rng(1)

% define null hyp
null_hyp     = freq;
null_hyp.powspctrm = zeros(size(null_hyp.powspctrm));

% define stat design
design      = zeros(2,size(freq.powspctrm,1)*2);
design(1,:) = repmat(1:size(freq.powspctrm,1),[1 2]);
design(2,:) = [ones(1,size(freq.powspctrm,1)),ones(1,size(freq.powspctrm,1))+1];

% define general stat config structure
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 'all';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.design              = design;
cfg.statistic           = 'ft_statfun_depsamplesT';  
cfg.latency             = [-0.8 0];
cfg.frequency           = [6 14];
stat                    = ft_freqstatistics(cfg,freq,null_hyp);

% get bayes at peak
addpath('C:\Users\ra34fod\github\bayesFactor\')
pidx = find(max(stat.stat(:))==stat.stat(:));
[~,y,z] = ind2sub(size(stat.stat),pidx);
tidx = stat.time(z); fidx = stat.freq(y);
stat.bayes = bf.ttest(squeeze(freq.powspctrm(:,:,knnsearch(freq.freq',fidx),knnsearch(freq.time',tidx))));

% get effect size at peak
stat.dz = abs(stat.stat ./ sqrt(size(freq.powspctrm,1)));

% report result
fprintf('\n--- Statistics ---\nt(%d): %3.3f\np: %3.3f\nBF: %3.3f\nd: %3.3f\n\n',size(freq.powspctrm,1)-1,stat.stat(pidx),stat.posclusters(1).prob,stat.bayes,stat.dz(pidx))

%% Visualise Effect
% predefine group output
group_freq = cell(6,1);

% cycle through files
for patient = 1 : 6
                               
    % define patient directory
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,patient);
    filename = sprintf('%s/sub-%02.0f_ephys-bipolar.mat',directory,patient);
    
    % load data
    load(filename,'ephys')

    % select MD
    ephys = ft_selectdata(struct('channel','MD'),ephys);
    
    % record trialinfo
    trialinfo = ephys.trialinfo(:,3);

    % get phase of ROI channels
    cfg         = [];
    cfg.output  = 'fourier';
    cfg.method  = 'wavelet';
    cfg.width   = 6;
    cfg.toi     = -0.8:0.025:0.8;
    cfg.foi     = 5:0.5:20;
    cfg.pad     = 'nextpow2';
    freq        = ft_freqanalysis(cfg,ephys); clear ephys
    freqs       = freq.freq;
    times       = freq.time;
    power       = abs(freq.fourierspctrm);

    % get hit > miss power difference
    power = permute(mean(power(trialinfo==1,:,:,:)) - mean(power(trialinfo==0,:,:,:)),[2 3 4 1]);
    
    % create frequency structure
    group_freq{patient} = struct('label',{{'MD'}},'freq',freqs,'time',times,...
                  'dimord','chan_freq_time','cfg',[],'powspctrm',power);
end

% get grand average of participants
cfg                 = [];
cfg.keepindividual  = 'yes';
freq                = ft_freqgrandaverage(cfg,group_freq{:});

% get key variables
vals = interp2(squeeze(mean(freq.powspctrm)),4);
freqs = linspace(freq.freq(1),freq.freq(end),size(vals,1));
times = linspace(freq.time(1),freq.time(end),size(vals,2));
imagesc(times,freqs,vals); axis xy
colormap(flipud(brewermap(11,'RdBu')))
maxabs = max(abs(vals(:)));
caxis([-maxabs maxabs])
xline(0,'k--')


