%% Prepare Workspace
% clear workspace
clearvars
close all
clc

% define data directory
dir_data = 'C:/Users/ra34fod/Documents/meg-ieeg_data_v2/';
dir_repo = 'C:/Users/ra34fod/github/thalamic_detection/';
cd(dir_repo);

% add toolboxes
addpath(genpath(sprintf('%ssubfunctions',dir_repo)))

% define participants
npp = 6;

%% Calculate Phase-Locking Value between Thalamus and Source Sensors
% cycle through participants
for pp = 1 : npp
                  
    % define dataset directory
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,pp);
    
    % estimate connectivity for MD
    run_connectivity_estimation(directory,pp);    
end

% run main statistics
run_connectivity_statistics(dir_data,6,'sourceMD')

%% Get Phase Slope Index
% cycle through participants
for pp = 1 : npp
               
    % estimate connectivity for MD
    run_psi_estimation(dir_data,pp);    
end

%% Run Statistics
% run MD statistics
run_psi_statistics(dir_data,6)

