function run_phaseBifurcation_contrastStatistics(directory)
    
% define number of patients and controls
npp = 5;

% predefine cell for group data
group_data = {};
    
% cycle through participants 
for pp = 1 : npp
    
    % load data
    pp_directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',directory,pp);
    filename = sprintf('sub-%02.0f_phaseBif-MD.mat',pp);
    load([pp_directory,filename],'freq')

    % add data to group
    group_data{pp,1} = freq;
    clear freq
    
    % load data
    filename = sprintf('sub-%02.0f_phaseBif-ANT.mat',pp);
    load([pp_directory,filename],'freq')

    % add data to group
    group_data{pp,2} = freq;
    clear freq
end
  
% get grand average of participants
cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter = {'pbi','pbi_raw','pbi_dist'};
freq1 = ft_freqgrandaverage(cfg,group_data{:,1});
freq2 = ft_freqgrandaverage(cfg,group_data{:,2});

% create contrast
freq = freq1; 
freq.pbi = freq1.pbi - freq2.pbi;
 
% set random seed
rng(1)

% define null hyp
null_hyp     = freq;
null_hyp.pbi = zeros(size(null_hyp.pbi));

% define stat design
design      = zeros(2,size(freq.pbi,1)*2);
design(1,:) = repmat(1:size(freq.pbi,1),[1 2]);
design(2,:) = [ones(1,size(freq.pbi,1)),ones(1,size(freq.pbi,1))+1];

% load MD cluster
load(sprintf('%s/derivatives/group/stat_pbi-MD.mat',directory),'stat');

% define general stat config structure
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 'all';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.parameter           = 'pbi';
cfg.design              = design;
cfg.statistic           = 'ft_statfun_depsamplesT';  
cfg.tail                = 1;
% cfg.avgoverfreq         = 'yes';
% cfg.avgovertime         = 'yes';
cfg.latency             = [min(stat.time(any(any(stat.posclusterslabelmat==1,2),1))) max(stat.time(any(any(stat.posclusterslabelmat==1,2),1)))];
cfg.frequency           = [min(stat.freq(any(any(stat.posclusterslabelmat==1,3),1))) max(stat.freq(any(any(stat.posclusterslabelmat==1,3),1)))];
stat                    = ft_freqstatistics(cfg,freq,null_hyp);

% get bayes at peak
addpath('C:\Users\ra34fod\github\bayesFactor\')
pidx = find(max(stat.stat(:))==stat.stat(:));
[~,y,z] = ind2sub(size(stat.stat),pidx);
tidx = stat.time(z); fidx = stat.freq(y);
stat.bayes = bf.ttest(squeeze(freq.pbi(:,:,knnsearch(freq.freq',fidx),knnsearch(freq.time',tidx))));

% get effect size at peak
stat.dz = abs(stat.stat ./ sqrt(size(freq.pbi,1)));

% report result
fprintf('\n--- Statistics: "MD > ANT" ---\nt(%d): %3.3f\np: %3.3f\nBF: %3.3f\nd: %3.3f\n\n',size(freq.pbi,1)-1,stat.stat(pidx),stat.posclusters(1).prob,stat.bayes,stat.dz(pidx))

% save outcome
save(sprintf('%s/derivatives/group/stat_pbi-thalContrast.mat',directory),'stat','freq','-v7.3');

end

