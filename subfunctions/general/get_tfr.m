function [phase,trialinfo,label,freqs,times,power] = get_tfr(ephys,channels)
    
% if channels is not specified
if ~exist('channels','var'); channels = ephys.label; end

% record trialinfo
trialinfo = ephys.trialinfo(:,3);

% get phase of ROI channels
cfg         = [];
cfg.output  = 'fourier';
cfg.method  = 'wavelet';
cfg.width   = 6;
cfg.toi     = -0.8:0.05:0.8;
cfg.foi     = 4:20;
cfg.pad     = 'nextpow2';
cfg.channel = channels;
freq        = ft_freqanalysis(cfg,ephys); clear ephys
freqs       = freq.freq;
times       = freq.time;
label       = freq.label;

% extract key info
phase = single(angle(freq.fourierspctrm));
if nargout == 6; power = single(abs(freq.fourierspctrm)); end
clear freq

