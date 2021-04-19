%% Prepare Workspace
% clear workspace
clearvars
close all
clc

% define data directory
dir_data = 'C:/Users/ra34fod/Documents/meg-ieeg_data_v2/';
dir_repo = 'C:/Users/ra34fod/github/thalamic_detection/';

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
    
    % plot MD
    ms = [-1 -1 1 1 1 1];
    figure; 
    plot(tml.time,tml.avg(ismember(tml.label,'MD'),:).*ms(patient),'b-','linewidth',1);
    xlim([-0.25 1]); xline(0,'k--');
    xlabel('Time (s)'); ylabel('Amp. (a.u.)')
    set(gca,'box','off','tickdir','out','fontsize',16)
    
    % plot ANT
    as = [1 1 1 1 -1];
    if patient < 6
        figure; 
        plot(tml.time,tml.avg(ismember(tml.label,'ANT'),:).*as(patient),'r-','linewidth',1);
        xlim([-0.25 1]); xline(0,'k--');
        xlabel('Time (s)'); ylabel('Amp. (a.u.)')
        set(gca,'box','off','tickdir','out','fontsize',16)
    end
end
    