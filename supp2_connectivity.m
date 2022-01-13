%% supp2_connectivity
% clear workspace
clearvars
close all
clc

% define data directory
dir_data = 'F:/meg-ieeg_data_v2/';
dir_repo = 'C:/Users/ra34fod/github/corticothalamic-connect/';
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
      
    % prepare filename   
    filename = sprintf('sub-%02.0f_ephys-sourceMD.mat',pp);
    
    % load data
    fprintf('loading data...\n')
    load([directory,filename],'source')    

    % load mask
    load(sprintf('%s/derivatives/group/stat_pli-all_sourceMD.mat',dir_data),'stat_mask'); %#ok<*LOAD>

    % get sourcemodel    
    [atlas,sourcemodel] = get_sourcemodel;

    % define channels
    chans = cat(1,source.label(stat_mask.stat(stat_mask.inside)>2),{'MD'});
    
    % record trialinfo
    trialinfo = source.trialinfo(:,3);

    % get phase of ROI channels
    cfg         = [];
    cfg.output  = 'fourier';
    cfg.method  = 'wavelet';
    cfg.width   = 6;
    cfg.toi     = -0.8:0.05:0.8;
    cfg.foi     = 5:0.5:20;
    cfg.pad     = 'nextpow2';
    cfg.channel = chans;
    freq        = ft_freqanalysis(cfg,source); clear source
    freqs       = freq.freq;
    times       = freq.time;
    label       = freq.label;

    % extract key info
    phase = single(angle(freq.fourierspctrm));
    clear freq
        
    % compute connectivity for every trial and channel
    [connect_all,connect_beh,~,connect_perm] = getPhaseLockingValue(phase,trialinfo);
    
    % store group PLV
    group_all(pp,:,:) = mean(connect_all,1);
    group_beh(pp,:,:) = mean(connect_beh,1);
        
    % store permutation data 
    [~,peakfreq] = max(mean(mean(connect_all(:,:,times<0),3),1));
    vals(pp,:) = squeeze(mean(mean(connect_perm(:,peakfreq,times<0,:),3),1));
    obs(pp,:) = cat(1,zeros(size(vals,2)-1,1),1);
    pf(pp,:) = zeros(1,size(vals,2))+round(freqs(peakfreq));
    
    % get mean phase angle/vector
    mean_angle = abs(mean(exp(1i.*(circ_mean(phase(:,:,peakfreq,times<0),[],4)))));
    mean_angle(mean_angle==1) = NaN; 
    [~,peakchan] = max(mean_angle);
    res_phase = circ_mean(circ_mean(phase(:,peakchan,peakfreq,times<0),[],4),[],2);
    plv_angle(pp,1) = angle(mean(exp(1i.*res_phase)));
    plv_length(pp,1) = abs(mean(exp(1i.*res_phase)));
end

% prepare angular data
f2_c_tbl = table({'1','2','3','4','5','6','mean'}',[plv_angle; circ_mean(plv_angle)],[plv_length; abs(mean(exp(1i.*plv_angle)))],...
                 'VariableNames',{'subj','angle','length'});
writetable(f2_c_tbl,[dir_repo,'/source_data/fig_2.xlsx'],'Sheet','fig_2c');

% get mean values
avg_all = interp2(squeeze(mean(group_all)),4);
avg_beh = interp2(squeeze(mean(group_beh)),4);
freqs = linspace(freqs(1),freqs(end),size(avg_all,1));
times = linspace(times(1),times(end),size(avg_all,2));

% convert data
xls_data{1} = tfr_to_longform(avg_all,freqs,times);
xls_data{2} = tfr_to_longform(avg_beh,freqs,times);

% write tfr xls data for figure 2a
f2_a_tbl = table(xls_data{1}(:,1),xls_data{1}(:,2),xls_data{1}(:,3),...
    'VariableNames',{'all_plv','freq','time'});
writetable(f2_a_tbl,[dir_repo,'/source_data/fig_2.xlsx'],'Sheet','fig_2a')

% write tfr xls data for supplementary figure 7
s8_a_tbl = table(xls_data{1}(:,1),xls_data{1}(:,2),xls_data{1}(:,3),...
    'VariableNames',{'all_plv','freq','time'});
s8_b_tbl = table(xls_data{2}(:,1),xls_data{2}(:,2),xls_data{2}(:,3),...
    'VariableNames',{'beh_plv','freq','time'});
writetable(s8_a_tbl,[dir_repo,'/source_data/supp_fig_9.xlsx'],'Sheet','supp_fig_8a')
writetable(s8_b_tbl,[dir_repo,'/source_data/supp_fig_9.xlsx'],'Sheet','supp_fig_8b')

% write histogram xls data for figure 2d
f2_d_tbl = array2table(cat(1,vals,pf,obs(1,:))','VariableNames',...
             {'vals_pp1','vals_pp2','vals_pp3','vals_pp4','vals_pp5','vals_pp6',...
             'freq_pp1','freq_pp2','freq_pp3','freq_pp4','freq_pp5','freq_pp6',...
             'obs_boolean'});
writetable(f2_d_tbl,[dir_repo,'/source_data/fig_2.xlsx'],'Sheet','fig_2d')

%% Record Data for PSI Participant Plots
% load main effect
load(sprintf('%s/derivatives/group/stat_psi-MD.mat',dir_data),'stat');

% cycle through participants
for pp = 1 : 6
    
    % load data
    load(sprintf('%s/derivatives/sub-%02.0f/ephys/sub-%02.0f_ephys-psiMD_contrast.mat',dir_data,pp,pp),'data_all','perm_all')
    
    % restrict data
    chan_idx = stat{1}.posclusterslabelmat(stat{1}.inside)==1;
    [~,peak_idx] = max(nanmean(data_all.psi_z(chan_idx,freqs<=20,1)));
    
    % get key data    
    perm_avg = squeeze(nanmean(nanmean(perm_all(chan_idx,peak_idx,:),1),2));
    obs_avg = squeeze(nanmean(nanmean(data_all.psi(chan_idx,peak_idx,:),1),2));
    
    % store vals
    vals(:,pp) = cat(1,obs_avg,perm_avg);
end

% add boolean
vals(:,7) = cat(1,1,zeros(200,1));

% save table
tbl = array2table(vals,'variableNames',{'vals_pp1','vals_pp2','vals_pp3','vals_pp4','vals_pp5','vals_pp6','obs_boolean'});
writetable(tbl,[dir_repo,'/source_data/fig_2.xlsx'],'Sheet','fig_2f')
