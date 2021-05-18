function run_psi_statistics(directory,npp)
    
% predefine cell for group data
group_data = cell(npp,1);
    
% cycle through participants 
for pp = 1 : npp
        
    % load data
    load(sprintf('%s/derivatives/sub-%02.0f/ephys/sub-%02.0f_ephys-psiMD_contrast.mat',directory,pp,pp),'data_hits','data_misses','data_all');
    
    % add data to group
    group_data{pp,1} = data_hits;
    group_data{pp,2} = data_misses;
    group_data{pp,3} = data_all;
    clear data
end

% get grand average of participants
cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter = {'psi','psi_z'};
grand_hits = ft_freqgrandaverage(cfg,group_data{:,1});
grand_misses = ft_freqgrandaverage(cfg,group_data{:,2});
grand_all = ft_freqgrandaverage(cfg,group_data{:,3});

%% convert to source data
% get sourcemodel
[atlas,sourcemodel,partial] = get_sourcemodel();

% create source structure
source = struct('dim',sourcemodel.dim,'inside',atlas.inside(:),'pos',sourcemodel.pos.*10,...
                'label',{grand_all.label},'all',zeros(size(grand_all.psi_z,1),numel(atlas.inside)),...
                'diff',zeros(size(grand_all.psi_z,1),numel(atlas.inside)),...
                'unit','mm','dimord','rpt_pos','transform',atlas.transform);

% add in measure
source.all(:,source.inside) = mean(mean(grand_all.psi_z(:,:,grand_all.freq>=6&grand_all.freq<=14),4),3);
source.diff(:,source.inside) = mean(mean(grand_hits.psi_z(:,:,grand_all.freq>=6&grand_all.freq<=14),4),3) - mean(mean(grand_misses.psi_z(:,:,grand_all.freq>=6&grand_all.freq<=14),4),3);

% load in partial mask
source.inside = partial.inside(:);
source.all(:,source.inside~=1) = NaN;
source.diff(:,source.inside~=1) = NaN;

% rename source to fit code
freq = source;

%% Run Stats (All > Chance)
% set random seed
rng(1)

% define null hyp
null_hyp     = freq;
null_hyp.all = zeros(size(null_hyp.all));
null_hyp.diff = zeros(size(null_hyp.diff));

% define stat design
design      = zeros(2,size(freq.all,1)*2);
design(1,:) = repmat(1:size(freq.all,1),[1 2]);
design(2,:) = [ones(1,size(freq.all,1)),ones(1,size(freq.all,1))+1];

% define general stat config structure
stat                    = {};
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 'all';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.design              = design;
cfg.statistic           = 'ft_statfun_depsamplesT';  
cfg.parameter           = 'all';
stat{1}                 = ft_sourcestatistics(cfg,freq,null_hyp);
cfg.parameter           = 'diff';
stat{2}                 = ft_sourcestatistics(cfg,freq,null_hyp);

% get bayes at peak
cond_str = {'all','h>m'};
for i = 1 : 2
    addpath('C:\Users\ra34fod\github\bayesFactor\')
    pidx = find(max(stat{i}.stat(:))==stat{i}.stat(:));
    stat{i}.bayes = bf.ttest(squeeze(freq.diff(:,pidx)));

    % get effect size at peak
    stat{i}.dz = abs(stat{i}.stat ./ sqrt(size(freq.all,1)));

    % report result
    fprintf('\n--- Statistics: "%s" for "all trials" ---\nt(%d): %3.3f\np: %3.3f\nBF: %3.3f\nd: %3.3f\n\n',cond_str{i},size(freq.all,1)-1,stat{i}.stat(pidx),stat{i}.posclusters(1).prob,stat{i}.bayes,stat{i}.dz(pidx))
end

%% Save
% save outcome
save(sprintf('%s/derivatives/group/stat_psi-MD.mat',directory),'stat','grand_hits','grand_misses','grand_all');

end

