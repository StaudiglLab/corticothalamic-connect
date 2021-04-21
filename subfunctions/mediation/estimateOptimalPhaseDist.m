function [phase_lag,cidx] = estimateOptimalPhaseDist(ephys,chan)

if nargin == 1; chan = 'all'; end

% get phase of ROI channels
cfg             = [];
cfg.keeptrials  = 'yes';
cfg.output      = 'fourier';
cfg.method      = 'wavelet';
cfg.pad         = 'nextpow2';
cfg.width       = 6;
cfg.toi         = -0.8:0.05:0;
cfg.foi         = 6:1:14;
cfg.channel     = chan;
freq            = ft_freqanalysis(cfg,ephys);
phase           = angle(freq.fourierspctrm);
trlinfo         = freq.trialinfo(:,3);

% --- get peak PBI band --- % 
% split phase data
phase_hits = phase(trlinfo==1,:,:,:);
phase_misses = phase(trlinfo==0,:,:,:);

% get trials
A = phase_hits;
B = phase_misses;
X = cat(1,A,B);

% get PLV
Ap = permute(abs(mean(exp(1i.*A))),[2 3 4 1]);
Bp = permute(abs(mean(exp(1i.*B))),[2 3 4 1]);
Xp = permute(abs(mean(exp(1i.*X))),[2 3 4 1]);

% get phase bifurcation
pbf = (Ap - Xp) .* (Bp - Xp);

% get peak pbf
for chan = 1 : size(pbf,1)
    
    % get peak pbf
    freq_pbf = squeeze(mean(mean(pbf(chan,:,:),3),1));
    [~,fidx] = max(freq_pbf);
    time_pbf = squeeze(mean(pbf(chan,fidx,:),1));
    [~,tidx] = max(time_pbf);
    %if size(pbf,1)>1; fidx = knnsearch(freq.freq',8);
    %else; fidx = knnsearch(freq.freq',13);
    %end
    
    % restrict phase to peak 
    phase_peak(:,chan) = phase(:,chan,fidx,tidx);
    pbf_peak(chan) = mean(pbf(chan,fidx,tidx),3);
end

% get channel peak
[~,cidx] = max(pbf_peak);

% ------ % 

% restrict to peak phase
phase = phase_peak;

% get optimal phase for hits
A = phase(trlinfo==1,:,:,:);
X = repmat(angle(mean(exp(1i.*A))),[size(phase,1) 1 1 1]);

% get resultant vector
Y = cat(ndims(phase)+1,phase,X);
phase_lag = abs(mean(exp(1i.*Y),ndims(Y)));
