
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
    MV = estimateOptimalPhaseDist(source,'MD');
    IV = estimateOptimalPhaseDist(source,source.label(mpfc_mask));%
    
    % get behaviour
    DV = source.trialinfo(:,3)+1;
    
    % drop NaNs
    MV(isnan(DV),:) = [];
    IV(isnan(DV),:) = [];
    DV(isnan(DV),:) = [];
    
    %% GLM: Does Thalamic Bifurcation Mediate Performance
    % cycle through source channels
    for chan = 1 : size(IV,2)
    
        % define matrices
        y = DV;
        X = cat(2,IV(:,chan),MV);

        % run mediation analysis
        stat = runMediationAnalysis(X,y);

        % add to group
        glm_beh(chan,1) = stat{1}; % MD mediates mPFC
        glm_beh(chan,2) = stat{2}; % mPFC mediates MD
        
        % run permutations
        for p = 1 : 1000

            % shuffle thalamus data
            Xi = X;
            Xi(:,2) = X(randperm(numel(y)),2);

            % run change prediction
            stat = runMediationAnalysis(Xi,y);

            % add to group
            for i=1:2
                perm_a(p,i) = stat{i}.a;
                perm_b(p,i) = stat{i}.b;
                perm_ab(p,i) = stat{i}.ab;
                perm_d(p,i) = stat{i}.d;
                perm_cdash(p,i) = stat{i}.cdash;
            end
        end

        % z-transform
        for i = 1:2
            glm_z(pp,chan,i).a = (glm_beh(chan,i).a - mean(perm_a(:,i))) ./ std(perm_a(:,i));
            glm_z(pp,chan,i).b = (glm_beh(chan,i).b - mean(perm_b(:,i))) ./ std(perm_b(:,i));
            glm_z(pp,chan,i).ab = (glm_beh(chan,i).ab - mean(perm_ab(:,i))) ./ std(perm_ab(:,i));
            glm_z(pp,chan,i).d = glm_beh(chan,i).d;%(glm_beh(chan,i).d - mean(perm_d(:,i))) ./ std(perm_d(:,i));
            glm_z(pp,chan,i).cdash = (glm_beh(chan,i).cdash - mean(perm_cdash(:,i))) ./ std(perm_cdash(:,i));
        end

        % store mean permuted M for channel
        permutedVals(chan,:,1) = perm_a(:,1);
        permutedVals(chan,:,2) = perm_b(:,1);
        permutedVals(chan,:,3) = perm_ab(:,1);
        permutedVals(chan,:,4) = perm_d(:,1);
        permutedVals(chan,:,5) = perm_cdash(:,1);
        
        % update user
        fprintf('channel %d of %d complete...\n',chan,size(IV,2))
    end   
    
    %% Plot % grab obs beh
    obs_beh = cat(1,[glm_beh(:,2).a],[glm_beh(:,2).b],[glm_beh(:,2).ab],[glm_beh(:,2).d],[glm_beh(:,2).cdash])';
    label = {'a','b','ab','d','c'''};
    xlims = [-0.3 1; -1 6; -0.2 1; -3 3; 5 15];
    binsize = [0.05 0.05 0.05 0.05 0.05];
    
    % cycle through plots
    for i = 1 : numel(label)

        % get vals        
        perm_avg = squeeze(mean(permutedVals(:,:,i),1)./std(permutedVals(:,:,i),[],1));
        obs_avg = mean(obs_beh(:,i))./std(obs_beh(:,i),[],1);
        
        % plot histogram
        figure; hold on
        bins = xlims(i,1):binsize(i):xlims(i,2);
        histogram(perm_avg,'normalization','probability','BinEdges',bins);
        title(label{i});
        xlim(xlims(i,:)); ylim([0 0.4])
        xline(obs_avg,'k-'); set(gca,'tickdir','out','box','off');
        %title(sprintf('sub-%02.0f - peak freq: %2.1fHz',pp,freqs(peak_idx)))
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
for i=1:5; stat.bayes(i) = bf.ttest(squeeze(freq.powspctrm(:,i))); end
