%% GLM Analyses
% clear workspace
clearvars
close all
clc

% define data directory
dir_data = 'C:/Users/ra34fod/Documents/meg-ieeg_data_v2/';
dir_repo = 'C:/Users/ra34fod/github/corticothalamic-connect/';

% add toolboxes
addpath(genpath(sprintf('%ssubfunctions',dir_repo)))
addpath('C:/Users/ra34fod/github/export_fig/')
addpath('C:/Users/ra34fod/github/circstat-matlab/')

% remove interfering fieldtrip function "hann()"
rmpath('C:\Users\ra34fod\github\fieldtrip\external\signal\')

% define participants
npp = 6;

%% Run Mediation 
% cycle through participants
for pp = 1 : 6
                        
    % define source dataset
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,pp);
    filename = sprintf('sub-%02.0f_ephys-sourceMD.mat',pp);
    load([directory,filename],'source')

    % get labels inside source
    [atlas,sourcemodel] = get_sourcemodel();
    source = ft_selectdata(struct('channel',{cat(1,source.label(atlas.inside(sourcemodel.inside(:))),'MD')}),source); 

    % get mpfc mask
    load('C:\Users\ra34fod\Documents\meg-ieeg_data_v2\derivatives\group\stat_pli-all_sourceMD.mat','stat_mask'); %#ok<*LOAD>
    mpfc_mask = stat_mask.stat>5;
    mpfc_mask = mpfc_mask(stat_mask.inside);
    
    %% Prepare Data
    % estimate distance to optimal phase
    MD_PBI = estimateOptimalPhaseDist(source,'MD');
    mPFC_PBI = estimateOptimalPhaseDist(source,source.label(mpfc_mask));
    
    % get behaviour
    percept_performance = source.trialinfo(:,3)+1;
    
    % drop NaNs
    MD_PBI(isnan(percept_performance),:) = [];
    mPFC_PBI(isnan(percept_performance),:) = [];
    percept_performance(isnan(percept_performance),:) = [];
    
    %% GLM: Does Thalamic Bifurcation Mediate Performance
    % cycle through source channels
    for chan = 1 : size(mPFC_PBI,2)
    
        % define matrices
        y = percept_performance;
        X = cat(2,mPFC_PBI(:,chan),MD_PBI);

        % run mediation analysis
        glm_beh(chan,1) = runMediationAnalysis(X,y);

        % run permutations
        for p = 1 : 500

            % shuffle thalamus data
            Xi = X;
            Xi(:,1) = X(randperm(numel(y)),1);
            Xi(:,2) = X(randperm(numel(y)),2);

            % run change prediction
            stat = runMediationAnalysis(Xi,y);

            % add to group
            perm_a(p,1) = stat.a;
            perm_b(p,1) = stat.b;
            perm_ab(p,1) = stat.ab;
            perm_d(p,1) = stat.d;
            perm_cdash(p,1) = stat.cdash;
        end

        % z-transform
        glm_z(pp,chan,1).a = (glm_beh(chan,1).a - mean(perm_a(:,1))) ./ std(perm_a(:,1));
        glm_z(pp,chan,1).b = (glm_beh(chan,1).b - mean(perm_b(:,1))) ./ std(perm_b(:,1));
        glm_z(pp,chan,1).ab = (glm_beh(chan,1).ab - mean(perm_ab(:,1))) ./ std(perm_ab(:,1));
        glm_z(pp,chan,1).d = (glm_beh(chan,1).d - mean(perm_d(:,1))) ./ std(perm_d(:,1));
        glm_z(pp,chan,1).cdash = (glm_beh(chan,1).cdash - mean(perm_cdash(:,1))) ./ std(perm_cdash(:,1));

        % store mean permuted M for channel
        permutedVals(chan,:,1) = perm_a(:,1);
        permutedVals(chan,:,2) = perm_b(:,1);
        permutedVals(chan,:,3) = perm_ab(:,1);
        permutedVals(chan,:,4) = perm_d(:,1);
        permutedVals(chan,:,5) = perm_cdash(:,1);
        
        % update user
        fprintf('channel %d of %d complete...\n',chan,size(mPFC_PBI,2))
    end   
    
    %% Plot % grab obs beh
    obs_beh = cat(1,[glm_beh(:,1).ab],[glm_beh(:,1).d],[glm_beh(:,1).cdash])';
    label = {'ab','d','c'''};
    %xlims = [-0.2 1; -3 3; 5 15];
    %binsize = [0.05 0.05 0.05];
    
    % cycle through plots
    for i = 1 : numel(label)

        % get vals        
        perm_avg = squeeze(mean(permutedVals(:,:,i))); %squeeze(mean(permutedVals(:,:,i),1)./std(permutedVals(:,:,i),[],1));
        obs_avg = mean(obs_beh(:,i)); %mean(obs_beh(:,i))./std(obs_beh(:,i),[],1);
        
        % plot histogram
        figure; hold on
        %bins = xlims(i,1):binsize(i):xlims(i,2);
        histogram(perm_avg,'normalization','probability');%,'BinEdges',bins);
        title(sprintf('sub-%02.0f - %s',pp,label{i}));
        %xlim(xlims(i,:)); ylim([0 0.4])
        xline(obs_avg,'k-'); set(gca,'tickdir','out','box','off');
    end
    
end
fprintf('\ndone');

%% Do Statistics
% create dummy structure
freq = struct('label',{{'a','b','ab','d','cdash'}},'time',1,'freq',1,...
                'powspctrm',nan(6,5),'cfg',[],'dimord','subj_chan_freq_time');

% add data
freq.powspctrm(:,1) = mean(reshape([glm_z(:,:,1).a],[size(glm_z,1) size(glm_z,2)]),2);
freq.powspctrm(:,2) = mean(reshape([glm_z(:,:,1).b],[size(glm_z,1) size(glm_z,2)]),2);
freq.powspctrm(:,3) = mean(reshape([glm_z(:,:,1).ab],[size(glm_z,1) size(glm_z,2)]),2);
freq.powspctrm(:,4) = mean(reshape([glm_z(:,:,1).d],[size(glm_z,1) size(glm_z,2)]),2);
freq.powspctrm(:,5) = mean(reshape([glm_z(:,:,1).cdash],[size(glm_z,1) size(glm_z,2)]),2);

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
%cfg.tail                = 1;
stat                    = ft_freqstatistics(cfg,freq,null_hyp);

% get bayes at peak
addpath('C:\Users\ra34fod\github\bayesFactor\')
pidx = find(max(stat.stat(:))==stat.stat(:));
for i=1:5; stat.bayes(i) = bf.ttest(squeeze(freq.powspctrm(:,i)),'tail','right'); end

%% Partial Correlation Analyses
% cycle through participants
for pp = 1 : 6
                        
    % define source dataset
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,pp);
    filename = sprintf('sub-%02.0f_ephys-sourceMD.mat',pp);
    load([directory,filename],'source')

    % get labels inside source
    [atlas,sourcemodel] = get_sourcemodel();
    source = ft_selectdata(struct('channel',{cat(1,source.label(atlas.inside(sourcemodel.inside(:))),'MD')}),source); 

    % get mpfc mask
    load('C:\Users\ra34fod\Documents\meg-ieeg_data_v2\derivatives\group\stat_pli-all_sourceMD.mat','stat'); %#ok<*LOAD>
    mpfc_mask = stat.stat>6;
    mpfc_mask = mpfc_mask(stat.inside);
    
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
        r(chan,1) = corr(mPFC_PBI(:,chan),percept_performance);
        r(chan,2) = partialcorr(mPFC_PBI(:,chan),percept_performance,MD_PBI);
        
        % run permutations
        for p = 1 : 1000

            % get shuffled corr
            IVr = mPFC_PBI(randperm(size(mPFC_PBI,1)),chan);
            rp(p,1) = corr(IVr,percept_performance);
            rp(p,2) = partialcorr(IVr,percept_performance,MD_PBI);
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

    % plot histogram
    figure; hold on
    bins = -0.006:0.0005:0.006;
    histogram(perm_avg,'normalization','probability','BinEdges',bins);
    xlim([-6 6].*10.^-3); ylim([0 0.4])
    xline(obs_avg,'k-'); set(gca,'tickdir','out','box','off');
    
    % store group z
    group_z(pp,1) = mean(z);
    fprintf('sub-%02.0f complete...\n',pp);
end
    
fprintf('\ndone');

%% Run Stats
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
%cfg.tail                = 1;
stat                    = ft_freqstatistics(cfg,freq,null_hyp);

% get bayes at peak
stat.bayes = bf.ttest(squeeze(freq.powspctrm),'tail','right');
