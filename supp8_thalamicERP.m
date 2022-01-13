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
end
    
% save data
tbl = array2table(cat(2,group_md,group_ant,tml.time'),'VariableNames',...
             {'md_pp1','md_pp2','md_pp3','md_pp4','md_pp5','md_pp6',...
             'ant_pp1','ant_pp2','ant_pp3','ant_pp4','ant_pp5','time'});
writetable(tbl,[dir_repo,'/source_data/supp_fig_12.xlsx'])
