function run_connectivity_statistics(directory,npp,roi_str)
   
%% Prepare Data
% predefine cell for group data
group_data = cell(npp,1);
    
% cycle through participants 
for pp = 1 : npp
    
    % load data
    pp_directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',directory,pp);
    filename = sprintf('sub-%02.0f_pliconnect-%s_PLV.mat',pp,roi_str);
    load([pp_directory,filename],'data')

    % add data to group
    group_data{pp} = data;
    clear data
end

% get grand average of participants
cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter = {'all','beh'}; 
grand_freq = ft_freqgrandaverage(cfg,group_data{:});

%% convert to source data
% get sourcemodel
[atlas,sourcemodel] = get_sourcemodel();

% create source structure
source = struct('dim',sourcemodel.dim,'inside',atlas.inside(:),'pos',sourcemodel.pos.*10,...
                'label',{grand_freq.label},'all',zeros(size(grand_freq.all,1),numel(atlas.inside)),...
                'beh',zeros(size(grand_freq.all,1),numel(atlas.inside)),...
                'unit','mm','dimord','rpt_pos','transform',atlas.transform);

% add in measure
source.all(:,source.inside) = mean(mean(grand_freq.all(:,:,grand_freq.freq<=14,grand_freq.time<0),4),3);
source.beh(:,source.inside) = mean(mean(grand_freq.beh(:,:,grand_freq.freq<=14,grand_freq.time<0),4),3);

% rename source to fit code
freq = source;

%% Run Stats (All > Chance)
% set random seed
rng(1)

% define null hyp
null_hyp     = freq;
null_hyp.all = zeros(size(null_hyp.all));

% define stat design
design      = zeros(2,size(freq.all,1)*2);
design(1,:) = repmat(1:size(freq.all,1),[1 2]);
design(2,:) = [ones(1,size(freq.all,1)),ones(1,size(freq.all,1))+1];

% define general stat config structure
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 'all';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.parameter           = 'all';
cfg.design              = design;
cfg.statistic           = 'ft_statfun_depsamplesT';  
cfg.tail                = 1;
stat                    = ft_sourcestatistics(cfg,freq,null_hyp);

% get bayes at peak
addpath('C:\Users\ra34fod\github\bayesFactor\')
pidx = find(max(stat.stat(:))==stat.stat(:));
stat.bayes = bf.ttest(squeeze(freq.all(:,pidx)));

% get effect size at peak
stat.dz = abs(stat.stat ./ sqrt(size(freq.all,1)));

% report result
fprintf('\n--- Statistics: "%s" for "all trials" ---\nt(%d): %3.3f\np: %3.3f\nBF: %3.3f\nd: %3.3f\n\n',roi_str,size(freq.all,1)-1,stat.stat(pidx),stat.posclusters(1).prob,stat.bayes,stat.dz(pidx))

%% Run Stats (Hits > Misses)
% restrict analysis to cluster
freq.beh(:,stat.posclusterslabelmat~=1) = 0;

% set random seed
rng(1)

% define null hyp
null_hyp     = freq;
null_hyp.beh = zeros(size(null_hyp.all));

% define stat design
design      = zeros(2,size(freq.beh,1)*2);
design(1,:) = repmat(1:size(freq.all,1),[1 2]);
design(2,:) = [ones(1,size(freq.all,1)),ones(1,size(freq.all,1))+1];

% define general stat config structure
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 'all';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.parameter           = 'beh';
cfg.design              = design;
cfg.statistic           = 'ft_statfun_depsamplesT';  
cfg.tail                = 1;
stat_beh                = ft_sourcestatistics(cfg,freq,null_hyp);

% get bayes at peak
addpath('C:\Users\ra34fod\github\bayesFactor\')
pidx = find(max(stat_beh.stat(:))==stat_beh.stat(:));
stat_beh.bayes = bf.ttest(squeeze(freq.beh(:,pidx)));

% get effect size at peak
stat_beh.dz = abs(stat_beh.stat ./ sqrt(size(freq.all,1)));

% report result
fprintf('\n--- Statistics: "%s" for "hits > misses" ---\nt(%d): %3.3f\np: %3.3f\nBF: %3.3f\nd: %3.3f\n\n',roi_str,size(freq.all,1)-1,stat_beh.stat(pidx),stat_beh.posclusters(1).prob,stat_beh.bayes,stat_beh.dz(pidx))

%% Get 8Hz Mask
% get sourcemodel
[atlas,sourcemodel] = get_sourcemodel();

% create source structure
source = struct('dim',sourcemodel.dim,'inside',atlas.inside(:),'pos',sourcemodel.pos.*10,...
                'label',{grand_freq.label},'all',zeros(size(grand_freq.all,1),numel(atlas.inside)),...
                'beh',zeros(size(grand_freq.all,1),numel(atlas.inside)),...
                'unit','mm','dimord','rpt_pos','transform',atlas.transform);

% add in measure
source.all(:,source.inside) = mean(mean(grand_freq.all(:,:,knnsearch(grand_freq.freq',8),grand_freq.time<0),4),3);
source.beh(:,source.inside) = mean(mean(grand_freq.beh(:,:,knnsearch(grand_freq.freq',8),grand_freq.time<0),4),3);

% rename source to fit code
freq = source;

% set random seed
rng(1)

% define null hyp
null_hyp     = freq;
null_hyp.all = zeros(size(null_hyp.all));

% define stat design
design      = zeros(2,size(freq.all,1)*2);
design(1,:) = repmat(1:size(freq.all,1),[1 2]);
design(2,:) = [ones(1,size(freq.all,1)),ones(1,size(freq.all,1))+1];

% define general stat config structure
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 'all';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.parameter           = 'all';
cfg.design              = design;
cfg.statistic           = 'ft_statfun_depsamplesT';  
cfg.tail                = 1;
stat_mask               = ft_sourcestatistics(cfg,freq,null_hyp);

%% Save Results
save(sprintf('%s/derivatives/group/stat_pli-all_sourceMD.mat',directory),'grand_freq','stat','stat_beh','stat_mask'); %#ok<*LOAD>

end
