%% Prepare Workspace
% clear workspace
clearvars
close all
clc

% define data directory
dir_data = 'F:/meg-ieeg_data_v2/';
dir_repo = 'C:\Users\ra34fod\github\corticothalamic-connect/';

% add toolboxes
addpath(genpath(sprintf('%ssubfunctions',dir_repo)))

% define participants
npp = 6;

%% Compute Phase Bifurcation 
% prepare files
files = cell(6,1);

% cycle through patients
for patient = 1 : 6
               
    % define patient directory
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,patient);
    filename = sprintf('%s/sub-%02.0f_ephys-bipolar.mat',directory,patient);
    load(filename)
    
    % filter
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 30;
    ephys = ft_preprocessing(cfg,ephys);
    
    % timelock data
    tml = ft_timelockanalysis([],ephys);
    tml = ft_timelockbaseline(struct('baseline',[-0.25 0]),tml);
    
    % store group MD
    ms = [-1 -1 1 1 1 1];
    group_md(:,patient) = tml.avg(ismember(tml.label,'MD'),:).*ms(patient);
    
    % store group ANT
    if patient ~= 6        
        as = [1 1 1 1 -1];
        group_ant(:,patient) = tml.avg(ismember(tml.label,'ANT'),:).*as(patient);
    end
    
    % get difference
    group_hits{patient} = ft_timelockanalysis(struct('channel','MD','trials',find(ephys.trialinfo(:,3)==1)),ephys);
    group_hits{patient} = ft_timelockbaseline(struct('baseline',[-0.25 0]),group_hits{patient});
    group_misses{patient} = ft_timelockanalysis(struct('channel','MD','trials',find(ephys.trialinfo(:,3)==0)),ephys);
    group_misses{patient} = ft_timelockbaseline(struct('baseline',[-0.25 0]),group_misses{patient});
end
    
% save data
tbl = array2table(cat(2,group_md,group_ant,tml.time'),'VariableNames',...
             {'md_pp1','md_pp2','md_pp3','md_pp4','md_pp5','md_pp6',...
             'ant_pp1','ant_pp2','ant_pp3','ant_pp4','ant_pp5','time'});
writetable(tbl,[dir_repo,'/source_data/supp_fig_12.xlsx'])

% test differences
cfg = struct('keepindividual','yes');
grand_hits = ft_timelockgrandaverage(cfg,group_hits{:});
grand_misses = ft_timelockgrandaverage(cfg,group_misses{:});
        
% set random seed
rng(1)

% define general stat config structure
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 2000;
cfg.avgoverchan         = 'yes';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.parameter           = 'individual';
cfg.design              = zeros(2,size(grand_hits.individual,1)*2);
cfg.design(1,:)         = repmat(1:size(grand_hits.individual,1),[1 2]);
cfg.design(2,:)         = [ones(1,size(grand_hits.individual,1)),ones(1,size(grand_hits.individual,1))+1];
cfg.latency             = [0 0.8];
cfg.statistic           = 'ft_statfun_depsamplesT';  
stat                    = ft_timelockstatistics(cfg,grand_hits,grand_misses);
