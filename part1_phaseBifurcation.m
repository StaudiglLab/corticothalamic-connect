%% Prepare Workspace
% clear workspace
clearvars
close all
clc

% define data directory
dir_data = 'F:/meg-ieeg_data_v2/'; % NOW H:/meg-ieeg_magdeburg_v4
dir_repo = 'C:/Users/ra34fod/github/corticothalamic-connect/';

% add toolboxes
addpath(genpath(sprintf('%ssubfunctions',dir_repo)))

% define participants
npp = 6;

%% Get Data Filenames
% prepare files
files = cell(6,1);

% cycle through patients
for patient = 1 : 6
               
    % define patient directory
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,patient);
    
    % add to list
    files{patient}   = sprintf('%s/sub-%02.0f_ephys-bipolar.mat',directory,patient);
    files{patient+6} = sprintf('%s/sub-%02.0f_ephys-source.mat',directory,patient);
end

% cycle through controls
for hlth = 1 : 12
              
    % define patient directory
    directory = sprintf('%s/derivatives/hlth-%02.0f/ephys/',dir_data,hlth);
    
    % add to list
    files{hlth+12} = sprintf('%s/hlth-%02.0f_ephys-sourceMAG.mat',directory,hlth);   
end
   
%% Compute Phase Bifurcation  
% cycle through files
for i = 1 : numel(files)
                 
    % estimate phase bifurcation for all files specified
    run_phaseBifurcation_estimation(files{i});
    
    % update user
    fprintf('loop %d complete...\n',i)
end

%% Run TFR Stats on Phase Bifurcation in Thalamus
% run phase bifurcation statistics
run_phaseBifurcation_statistics(dir_data,'MD')
run_phaseBifurcation_statistics(dir_data,'ANT')
run_phaseBifurcation_statistics(dir_data,'source')
run_phaseBifurcation_statistics(dir_data,'source','control')

% run MD > ANT contrast
run_phaseBifurcation_contrastStatistics(dir_data)

% run MD > ANT contrast
run_phaseBifurcation_patientControlContrast(dir_data)
