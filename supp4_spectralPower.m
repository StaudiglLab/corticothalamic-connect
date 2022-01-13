%% Prepare Workspace
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

%% Compute Phase Bifurcation  
% predefine group output
group_freq = cell(6,1);

% cycle through files
for patient = 1 : 6
          
    % --- thalamus --- %
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
    group_freq{patient,1} = struct('label',{label},'freq',freqs,'time',times,...
                  'dimord','chan_freq_time','cfg',[],'powspctrm',power);
              
    % --- mpfc --- %
    % load source dataset
    filename = sprintf('sub-%02.0f_ephys-sourceMD.mat',patient);
    load([directory,filename],'source')

    % get labels inside source
    [atlas,sourcemodel] = get_sourcemodel();
    source = ft_selectdata(struct('channel',{cat(1,source.label(atlas.inside(sourcemodel.inside(:))),'MD')}),source); 

    % get mpfc mask
    load('F:\meg-ieeg_data_v2\derivatives\group\stat_pli-all_sourceMD.mat','stat_mask'); %#ok<*LOAD>
    mpfc_mask = stat_mask.stat>5;
    mpfc_mask = mpfc_mask(stat_mask.inside);                        
    
    % select MD
    ephys = ft_selectdata(struct('channel',{source.label(mpfc_mask)}),source);
    clear source
    
    % get phase using variable width wavelets
    [~,trialinfo,label,freqs,times,power] = get_tfr(ephys,ephys.label);

    % get hit > miss power difference
    power = permute(mean(power(trialinfo==1,:,:,:)) - mean(power(trialinfo==0,:,:,:)),[2 3 4 1]);
    
    % create frequency structure
    group_freq{patient,2} = struct('label',{label},'freq',freqs,'time',times,...
                  'dimord','chan_freq_time','cfg',[],'powspctrm',power);
end

%% Run Thalamus Statistics
% get grand average of participants
cfg                 = [];
cfg.keepindividual  = 'yes';
freq                = ft_freqgrandaverage(cfg,group_freq{:,1});

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

%% Run mPFC Statistics
% average over channels
group_freq(:,2) = cellfun(@(x) {ft_selectdata(struct('avgoverchan','yes'),x)}, group_freq(:,2));

% get grand average of participants
cfg                 = [];
cfg.keepindividual  = 'yes';
freq                = ft_freqgrandaverage(cfg,group_freq{:,2});

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
pidx = find(max(abs(stat.stat(:)))==abs(stat.stat(:)));
[~,y,z] = ind2sub(size(stat.stat),pidx);
tidx = stat.time(z); fidx = stat.freq(y);
stat.bayes = bf.ttest(squeeze(freq.powspctrm(:,:,knnsearch(freq.freq',fidx),knnsearch(freq.time',tidx))));

% get effect size at peak
stat.dz = abs(stat.stat ./ sqrt(size(freq.powspctrm,1)));

% report result
fprintf('\n--- Statistics ---\nt(%d): %3.3f\np: %3.3f\nBF: %3.3f\nd: %3.3f\n\n',size(freq.powspctrm,1)-1,stat.stat(pidx),stat.negclusters(1).prob,stat.bayes,stat.dz(pidx))

%% Compute High-Res. for Visualisation 
% predefine group output
group_freq = cell(6,2);

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

    % get power of ROI channels
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

    % create frequency structure
    power_corr = (power - repmat(mean(power,4),[1 1 1 size(power,4)])) ./ repmat(std(power,[],4),[1 1 1 size(power,4)]);
    group_freq{patient,1} = struct('label',{{'MD'}},'freq',freqs,'time',times,...
                  'dimord','chan_freq_time','cfg',[],'powspctrm',permute(mean(power_corr,1),[2 3 4 1]));

    % get hit > miss power difference
    power = permute(mean(power_corr(trialinfo==1,:,:,:)) - mean(power_corr(trialinfo==0,:,:,:)),[2 3 4 1]);
    
    % create frequency structure
    group_freq{patient,2} = struct('label',{{'MD'}},'freq',freqs,'time',times,...
                  'dimord','chan_freq_time','cfg',[],'powspctrm',power);
end

% get grand average of participants
cfg                 = [];
cfg.keepindividual  = 'yes';
freq_raw            = ft_freqgrandaverage(cfg,group_freq{:,1});
freq_diff           = ft_freqgrandaverage(cfg,group_freq{:,2});

% plot raw
figure;
subplot(2,1,1);
vals = interp2(squeeze(mean(freq_raw.powspctrm)),4);
freqs = linspace(freq_raw.freq(1),freq_raw.freq(end),size(vals,1));
times = linspace(freq_raw.time(1),freq_raw.time(end),size(vals,2));
imagesc(times,freqs,vals); axis xy
colormap(flipud(brewermap(11,'RdBu')))
maxabs = round(max(abs(vals(:))),7);
caxis([-maxabs maxabs])
colorbar();
xline(0,'k--')

% convert diff power matrix to long form
xls_data{1} = tfr_to_longform(vals,freqs,times);

% plot difference
subplot(2,1,2);
vals = interp2(squeeze(mean(freq_diff.powspctrm)),4);
freqs = linspace(freq_diff.freq(1),freq_diff.freq(end),size(vals,1));
times = linspace(freq_diff.time(1),freq_diff.time(end),size(vals,2));
imagesc(times,freqs,vals); axis xy
colormap(flipud(brewermap(11,'RdBu')))
maxabs = round(max(abs(vals(:))),7);
caxis([-maxabs maxabs])
colorbar();
xline(0,'k--')

% convert diff power matrix to long form
xls_data{2} = tfr_to_longform(vals,freqs,times);

% --- mPFC Power for Visualisation --- %
% predefine group output
group_freq = cell(6,2);

% cycle through files
for patient = 1 : 6
       
    % define source dataset
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,patient);
    filename = sprintf('sub-%02.0f_ephys-sourceMD.mat',patient);
    load([directory,filename],'source')

    % get labels inside source
    [atlas,sourcemodel] = get_sourcemodel();
    source = ft_selectdata(struct('channel',{cat(1,source.label(atlas.inside(sourcemodel.inside(:))),'MD')}),source); 

    % get mpfc mask
    load('F:\meg-ieeg_data_v2\derivatives\group\stat_pli-all_sourceMD.mat','stat_mask'); %#ok<*LOAD>
    mpfc_mask = stat_mask.stat>5;
    mpfc_mask = mpfc_mask(stat_mask.inside);                        
    
    % select MD
    ephys = ft_selectdata(struct('channel',{source.label(mpfc_mask)}),source);
    clear source
    
    % record trialinfo
    trialinfo = ephys.trialinfo(:,3);

    % get power of ROI channels
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

    % create frequency structure
    power_corr = (power - repmat(mean(power,4),[1 1 1 size(power,4)])) ./ repmat(std(power,[],4),[1 1 1 size(power,4)]);
    group_freq{patient,1} = struct('label',{{'MD'}},'freq',freqs,'time',times,...
                  'dimord','chan_freq_time','cfg',[],'powspctrm',permute(mean(power_corr,1),[2 3 4 1]));

    % get hit > miss power difference
    power = permute(mean(power_corr(trialinfo==1,:,:,:)) - mean(power_corr(trialinfo==0,:,:,:)),[2 3 4 1]);
    
    % create frequency structure
    group_freq{patient,2} = struct('label',{{'MD'}},'freq',freqs,'time',times,...
                  'dimord','chan_freq_time','cfg',[],'powspctrm',power);
end

% get grand average of participants
cfg                 = [];
cfg.keepindividual  = 'yes';
freq_raw            = ft_freqgrandaverage(cfg,group_freq{:,1});
freq_diff           = ft_freqgrandaverage(cfg,group_freq{:,2});

% plot raw
figure;
vals = interp2(squeeze(mean(freq_raw.powspctrm)),4);
freqs = linspace(freq_raw.freq(1),freq_raw.freq(end),size(vals,1));
times = linspace(freq_raw.time(1),freq_raw.time(end),size(vals,2));
imagesc(times,freqs,vals); axis xy
colormap(flipud(brewermap(11,'RdBu')))
maxabs = round(max(abs(vals(:))),7);
caxis([-maxabs maxabs])
xline(0,'k--')

% convert raw power matrix to long form
xls_data{3} = tfr_to_longform(vals,freqs,times);

% plot difference
figure;
vals = interp2(squeeze(mean(freq_diff.powspctrm)),4);
freqs = linspace(freq_diff.freq(1),freq_diff.freq(end),size(vals,1));
times = linspace(freq_diff.time(1),freq_diff.time(end),size(vals,2));
imagesc(times,freqs,vals); axis xy
colormap(flipud(brewermap(11,'RdBu')))
maxabs = round(max(abs(vals(:))),7);
caxis([-maxabs maxabs])
xline(0,'k--')

% convert diff power matrix to long form
xls_data{4} = tfr_to_longform(vals,freqs,times);

% write xls data
tbl = table(xls_data{1}(:,1),xls_data{2}(:,1),xls_data{3}(:,1),xls_data{4}(:,1),xls_data{1}(:,2),xls_data{1}(:,3),...
    'VariableNames',{'thal_raw','thal_diff','mpfc_raw','mpfc_diff','freq','time'});
writetable(tbl,[dir_repo,'/source_data/supp_fig_4.xlsx'])


