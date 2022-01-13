%% Supplementary Mediation Analyses
% this script runs a partial correlation analysis to cross-check the
% mediation result presented in figure 3

% clear workspace
clearvars
close all
clc

% define data directory
dir_data = 'F:/meg-ieeg_data_v2/';
dir_repo = 'C:/Users/ra34fod/github/corticothalamic-connect/';

% add toolboxes
addpath(genpath(sprintf('%ssubfunctions',dir_repo)))
addpath('C:/Users/ra34fod/github/export_fig/')
addpath('C:/Users/ra34fod/github/circstat-matlab/')

% remove interfering fieldtrip function "hann()"
rmpath('C:\Users\ra34fod\github\fieldtrip\external\signal\')

% define participants
npp = 6;

%% Partial Correlation Analyses
% cycle through participants
vals = []; obs = [];
rng(1)
for pp = 1 : 6
                        
    % define source dataset
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,pp);
    filename = sprintf('sub-%02.0f_ephys-sourceMD.mat',pp);
    load([directory,filename],'source')

    % get labels inside source
    [atlas,sourcemodel] = get_sourcemodel();
    source = ft_selectdata(struct('channel',{cat(1,source.label(atlas.inside(sourcemodel.inside(:))),'MD')}),source); 

    % get mpfc mask
    load('F:\meg-ieeg_data_v2\derivatives\group\stat_pli-all_sourceMD.mat','stat_mask'); %#ok<*LOAD>
    mpfc_mask = stat_mask.stat>2;
    mpfc_mask = mpfc_mask(stat_mask.inside);
    
    %% Prepare Data
    % estimate distance to optimal phase
    MD_PBI = estimateOptimalPhaseDist(source,'MD');
    mPFC_PBI = estimateOptimalPhaseDist(source,source.label(mpfc_mask));%
    
    % get behaviour
    percept_performance = source.trialinfo(:,3)+1;
    
    % drop NaNs
    MD_PBI(isnan(percept_performance),:) = [];
    mPFC_PBI(isnan(percept_performance),:) = [];
    percept_performance(isnan(percept_performance),:) = [];
    
    %% Run Partial Corrs
    % cycle through source channels
    r = zeros(size(mPFC_PBI,2),2); rp = r;
    obs_d = zeros(size(mPFC_PBI,2),1);
    for chan = 1 : size(mPFC_PBI,2)
    
        % define matrices
        y = percept_performance;
        X = cat(2,mPFC_PBI(:,chan),MD_PBI);

        % get correlations
        r(chan,1) = corr(mPFC_PBI(:,chan),percept_performance,'type','Spearman');
        r(chan,2) = partialcorr(mPFC_PBI(:,chan),percept_performance,MD_PBI,'type','Spearman');
        
        % run permutations
        for p = 1 : 1000

            % get shuffled corr
            IVr = mPFC_PBI(randperm(size(mPFC_PBI,1)),chan);
            rp(p,1) = corr(IVr,percept_performance,'type','Spearman');
            rp(p,2) = partialcorr(IVr,percept_performance,MD_PBI,'type','Spearman');
        end

        % get z-transformed diff.
        d = r(chan,1) - r(chan,2);
        dp = rp(:,1) - rp(:,2);
        z(chan,1) = (d - mean(dp)) ./ std(dp);
        
        % store variables for plotting
        phist(chan,:) = dp;
        obs_d(chan,1) = d;
        
        % update user
        fprintf('channel %d of %d complete...\n',chan,size(mPFC_PBI,2))
    end   
    
    % get vals        
    perm_avg = mean(phist,1);
    obs_avg = mean(obs_d,1);

    % store data    
    vals(pp,:) = cat(1,obs_avg',perm_avg');
    obs(pp,:) = cat(1,1,zeros(numel(perm_avg),1));
    
    % plot histogram
    figure; hold on
   % bins = -0.006:0.0005:0.006;
    histogram(perm_avg,20,'normalization','probability');%,'BinEdges',bins);
 %   xlim([-8 8].*10.^-3); ylim([0 0.4])
    xline(obs_avg,'k-'); set(gca,'tickdir','out','box','off');
    
    % store group z
    group_z(pp,1) = mean(z);
    fprintf('sub-%02.0f complete...\n',pp);
end
    
fprintf('\ndone');

% create dummy structure
freq = struct('label',{{'z'}},'time',1,'freq',1,...
                'powspctrm',nan(6,1),'cfg',[],'dimord','subj_chan_freq_time');

% add data
freq.powspctrm = group_z;

% set random seed
rng(1)

% define null hyp
null_hyp    = freq;
null_hyp.powspctrm  = zeros(size(null_hyp.powspctrm));

% define stat design
design      = zeros(2,size(freq.powspctrm,1)*2);
design(1,:) = repmat(1:size(freq.powspctrm,1),[1 2]);
design(2,:) = [ones(1,size(freq.powspctrm,1)),ones(1,size(freq.powspctrm,1))+1];

% define general stat config structure
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'none';
cfg.numrandomization    = 'all';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.parameter           = 'powspctrm';
cfg.design              = design;
cfg.statistic           = 'ft_statfun_depsamplesT';  
cfg.tail                = 1;
stat                    = ft_freqstatistics(cfg,freq,null_hyp);

% get bayes at peak
stat.bayes = bf.ttest(squeeze(freq.powspctrm),'tail','right');

% save data
tbl = array2table(cat(1,vals,obs(1,:))','VariableNames',...
             {'vals_pp1','vals_pp2','vals_pp3','vals_pp4','vals_pp5','vals_pp6','obs_boolean'});
writetable(tbl,[dir_repo,'/source_data/supp_fig_11.xlsx'])
