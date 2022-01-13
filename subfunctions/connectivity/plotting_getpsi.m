
function [psi_hits,psi_misses,psi_all,freq_hits,freq_misses,freq] = plotting_getpsi(source,randomise,freq_hits,freq_misses,freq)
        
% compute power spectrm if not included as input
if nargin < 3

    % get fourier spectrum for hits
    cfg             = [];
    cfg.output      = 'fourier';
    cfg.channel     = {'S*','MD'};
    cfg.method      = 'mtmfft';
    cfg.taper       = 'hanning';
    cfg.keeptrials  = 'yes';
    cfg.foilim      = [0 100];
    cfg.pad         = 'nextpow2';
    freq            = ft_freqanalysis(cfg,source);

    % fix freq precision error
    freq.freq = linspace(freq.freq(1),freq.freq(end),numel(freq.freq));

    % split data
    cfg             = [];
    cfg.trials      = freq.trialinfo(:,3) == 1;
    freq_hits       = ft_selectdata(cfg,freq);
    cfg.trials      = freq.trialinfo(:,3) == 0;
    freq_misses     = ft_selectdata(cfg,freq);
end

% get channel combinations
channelcmb = cat(2,source.label(1:end-1),repmat({'MD'},[numel(source.label)-1 1]));
channelcmb = cat(1,channelcmb);

% if randomised order requested
if randomise
    
    % get MD position
    md_idx = ismember(freq_hits.label,'MD');
    
    % randomise hits and misses
    freq.fourierspctrm(:,md_idx,:) = freq.fourierspctrm(randperm(size(freq.fourierspctrm,1)),md_idx,:);
    freq_hits.fourierspctrm(:,md_idx,:) = freq_hits.fourierspctrm(randperm(size(freq_hits.fourierspctrm,1)),md_idx,:);
    freq_misses.fourierspctrm(:,md_idx,:) = freq_misses.fourierspctrm(randperm(size(freq_misses.fourierspctrm,1)),md_idx,:);
end

% run psi
cfg                     = [];
cfg.method              = 'psi';
cfg.channelcmb          = channelcmb;
cfg.bandwidth           = 5;
data_hits               = ft_connectivityanalysis(cfg,freq_hits);
psi_hits                = data_hits.psispctrm;

% run psi for misses [if requested, else save comp. time]
if nargout >= 2
    data_misses         = ft_connectivityanalysis(cfg,freq_misses);
    psi_misses          = data_misses.psispctrm;
    data_all         = ft_connectivityanalysis(cfg,freq);
    psi_all          = data_all.psispctrm;
end