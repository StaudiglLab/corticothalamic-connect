function run_psi_estimation(directory,pp)
    
% load data
load(sprintf('%s/derivatives/sub-%02.0f/ephys/sub-%02.0f_ephys-sourceMD.mat',directory,pp,pp),'source')  

% define window of interest
time_bins = [-0.8 0];

% define data dimensions
nperms = 200;
nfreqs = 103;

% predefine permutation output
perm_hits   = zeros(size(source.label,1)-1,nfreqs,nperms);
perm_misses = zeros(size(source.label,1)-1,nfreqs,nperms);
perm_all    = zeros(size(source.label,1)-1,nfreqs,nperms);

% select time
source = ft_selectdata(struct('latency',time_bins),source);  

% calculate freq-domain psi for hits and misses
[psi_hits,psi_misses,psi_all,freq_hits,freq_misses,freq_all] = internal_getpsi(source,false);

% cycle through permutations
for perm = 1 : nperms

    % calculate freq-domain psi for hits and misses
    [perm_hits(:,:,perm),perm_misses(:,:,perm),perm_all(:,:,perm)] = internal_getpsi(source,true,freq_hits,freq_misses,freq_all);
    
    % update user
    fprintf('\n\nsub-%02.0f: perm %d of %d complete...\n',pp,perm,nperms)
end

% z-transform result
psi_z_hits = (psi_hits - mean(perm_hits,3)) ./ std(perm_hits,[],3);
psi_z_misses = (psi_misses - mean(perm_misses,3)) ./ std(perm_misses,[],3);
psi_z_all = (psi_all - mean(perm_all,3)) ./ std(perm_all,[],3);

% create data structures
labels = source.label(1:end-1);
freqs = linspace(0,100,size(psi_hits,2)); % define frequencies resolved
data_hits = internal_createDataStruct(psi_hits,labels,freqs);
data_misses = internal_createDataStruct(psi_misses,labels,freqs);
data_all = internal_createDataStruct(psi_all,labels,freqs);

% if psi z computed
data_hits.psi_z = psi_z_hits;
data_misses.psi_z = psi_z_misses;
data_all.psi_z = psi_z_all;
% 
% % plot histogram
% [~,peak_idx] = max(mean(psi_z_hits(:,freqs<=14,1)));
% bins = -0.05:0.0025:0.0475;
% perm_avg = squeeze(mean(mean(perm_hits(:,peak_idx,1,:),1),2));
% obs_avg = squeeze(mean(mean(psi_hits(:,peak_idx,1,:),1),2));
% sorted_avg = sort(perm_avg);
% figure; hold on
% histogram(perm_avg,'BinEdges',bins,'normalization','probability');
% xlim([-0.05 0.05]); ylim([0 0.2])
% %xline(sorted_avg(round(numel(sorted_avg)*0.95)),'k--')
% xline(obs_avg,'k-'); set(gca,'tickdir','out');
% title(sprintf('sub-%02.0f - peak freq: %2.1fHz',pp,freqs(peak_idx)))

% save data
save(sprintf('%s/derivatives/sub-%02.0f/ephys/sub-%02.0f_ephys-psiMD_contrast-timeResolved.mat',directory,pp,pp),'data_hits','data_misses','data_all');

end

function [psi_hits,psi_misses,psi_all,freq_hits,freq_misses,freq] = internal_getpsi(source,randomise,freq_hits,freq_misses,freq)
      
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

end

function data = internal_createDataStruct(val,label,freq)
    
data = struct('label',{label},'freq',freq,'time',1,'dimord','chan_freq_time',...
              'cfg',[],'psi',val);
              %'cfg',[],'gc_thal2Source',val(2:2:end,:,:),'gc_source2Thal',val(1:2:end,:,:));
    
end
