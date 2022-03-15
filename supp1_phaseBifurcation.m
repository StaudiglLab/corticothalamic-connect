%% Prepare Workspace
% this script runs additional analyses required for the plots of figure 1

% clear workspace
clearvars
close all
clc

% define data directory
dir_data = 'F:/meg-ieeg_data_v2/';
dir_repo = 'C:/Users/ra34fod/github/corticothalamic-connect/';

% add toolboxes
addpath(genpath(sprintf('%ssubfunctions',dir_repo)))

% define participants
npp = 6;

%% Compute High-Res PBI TFR for Thalamus
% cycle through participants
group_md = [];
group_ant = [];
for pp = 1 : 6
                
    % define dataset
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,pp);
    filename = sprintf('sub-%02.0f_ephys-bipolar.mat',pp);
    
    % load data
    load([directory,filename],'ephys')
    
    % record trialinfo
    trialinfo = ephys.trialinfo(:,3);

    % define channels
    if pp ~= 6
        chans = {'MD','ANT'};
    else
        chans = 'MD';
    end
    
    % get phase of ROI channels
    cfg         = [];
    cfg.output  = 'fourier';
    cfg.method  = 'wavelet';
    cfg.width   = 6;
    cfg.toi     = -0.8:0.025:0.8;
    cfg.foi     = 5:0.5:20;
    cfg.pad     = 'nextpow2';
    cfg.channel = chans;
    freq        = ft_freqanalysis(cfg,ephys); clear ephys
    freqs       = freq.freq;
    times       = freq.time;
    label       = freq.label;

    % extract key info
    phase = single(angle(freq.fourierspctrm));
    
    % update user
    fprintf('\ncalculating phase bifurcation index...\n')

    % split phase data
    phase_hits = phase(trialinfo==1,:,:,:);
    phase_misses = phase(trialinfo==0,:,:,:);

    % get key variables
    ntimes = size(phase_hits,4);

    % get trials
    A = phase_hits;
    B = phase_misses;
    X = cat(1,A,B);

    % get PLV
    Ap = permute(abs(mean(exp(1i.*A))),[2 3 4 1]);
    Bp = permute(abs(mean(exp(1i.*B))),[2 3 4 1]);
    Xp = permute(abs(mean(exp(1i.*X))),[2 3 4 1]);

    % get phase bifurcation
    pbf = (Ap - Xp) .* (Bp - Xp);
    
    % store group
    group_md(:,:,pp) = pbf(ismember(freq.label,'MD'),:,:);
    if pp ~= 6
        group_ant(:,:,pp) = pbf(ismember(freq.label,'ANT'),:,:);
    else
    end
end

% --- prep. for figure 1 / supp. fig. 2 --- %
% store MD data
vals  = interp2(nanmean(group_md,3),4);
times = linspace(times(1),times(end),size(vals,2));
freqs = linspace(freqs(1),freqs(end),size(vals,1));
xls_data{1} = tfr_to_longform(vals,freqs,times);

% store ANT data
vals  = interp2(nanmean(group_ant,3),4);
xls_data{2} = tfr_to_longform(vals,freqs,times);

% store MD > ANT data
vals  = interp2(nanmean(group_md(:,:,1:5)-group_ant,3),4);
xls_data{3} = tfr_to_longform(vals,freqs,times);

% write xls data (for supplementary figure 2)
tbl = table(xls_data{1}(:,1),xls_data{2}(:,1),xls_data{3}(:,1),xls_data{1}(:,2),xls_data{1}(:,3),...
    'VariableNames',{'md_pbi','ant_pbi','md_vs_ant_pbi','freq','time'});
writetable(tbl,[dir_repo,'/source_data/supp_fig_2.xlsx'])

% hold onto table
fig_1c_tbl = tbl;


% --- supp. fig. 3 --- %
% store MD data
vals  = interp2(nanmean(group_md(:,:,[1 3 4 5 6]),3),4);
times = linspace(times(1),times(end),size(vals,2));
freqs = linspace(freqs(1),freqs(end),size(vals,1));
xls_data{1} = tfr_to_longform(vals,freqs,times);

% write xls data (for supplementary figure 2)
tbl = table(xls_data{1}(:,1),xls_data{1}(:,2),xls_data{1}(:,3),...
    'VariableNames',{'md_pbi','freq','time'});
writetable(tbl,[dir_repo,'/source_data/supp_fig_3.xlsx'])

%% Compute Alpha Phase Time Series
% cycle through participants
vals = []; obs = []; pf = [];
for pp = 1 : 6
                
    % define dataset
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,pp);
    filename = sprintf('sub-%02.0f_ephys-bipolar.mat',pp);
    
    % load data
    load([directory,filename],'ephys')
    
    % record trialinfo
    trialinfo = ephys.trialinfo(:,3);
    
    % filter signal
    ephys = ft_preprocessing(struct('bpfilter','yes','bpfreq',[6 9],'channel','MD'),ephys);
    
    % -------- %
    % split signal
    signal = reshape(cell2mat(ephys.trial),[numel(ephys.time{1}) numel(ephys.time)]);
    phase_hits = mean(signal(:,trialinfo==1),2);
    phase_misses = mean(signal(:,trialinfo==0),2);
    
    % get angle
    phase_hits = sin(angle(hilbert(phase_hits)));
    phase_misses = sin(angle(hilbert(phase_misses)));
    
%     % -------- %    
%     % cycle through trials and filter
%     signal = [];
%     for trl = 1 : numel(ephys.trial)
%         idx = ismember(ephys.label,'MD');
%         signal(trl,:) = angle(hilbert(ephys.trial{trl}(idx,:)));
%     end
%     
%     % split phase data
%     phase_hits = smooth(cos(circ_mean(signal(trialinfo==1,:,:,:))));
%     phase_misses = smooth(cos(circ_mean(signal(trialinfo==0,:,:,:))));
%     
%     % -------- %
    % store data    
    idx = ephys.time{1}>=-0.8&ephys.time{1}<=0.8;
    pp_pos = (1:2) + ((pp-1)*2);
    vals(:,pp_pos) = cat(2,phase_hits(idx),phase_misses(idx));
    times = ephys.time{1}(idx);
end

% save data
fig_1d_tbl = array2table(cat(2,vals,times'),'VariableNames',...
             {'vals_pp1_hits','vals_pp1_misses','vals_pp2_hits','vals_pp2_misses',...
             'vals_pp3_hits','vals_pp3_misses','vals_pp4_hits','vals_pp4_misses',...
             'vals_pp5_hits','vals_pp5_misses','vals_pp6_hits','vals_pp6_misses',...
             'times'});

%% Compute Observed vs. Chance Distribution Visualisation
% cycle through participants
vals = []; obs = []; pf = [];
for pp = 1 : 6
                
    % define dataset
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,pp);
    filename = sprintf('sub-%02.0f_ephys-bipolar.mat',pp);
    
    % load data
    load([directory,filename],'ephys')
    
    % record trialinfo
    trialinfo = ephys.trialinfo(:,3);
    
    % get phase of ROI channels
    cfg         = [];
    cfg.output  = 'fourier';
    cfg.method  = 'wavelet';
    cfg.width   = 6;
    cfg.toi     = -0.8:0.025:0;
    cfg.foi     = 6:0.5:10;
    cfg.pad     = 'nextpow2';
    cfg.channel = 'MD';
    freq        = ft_freqanalysis(cfg,ephys); clear ephys
    freqs       = freq.freq;
    times       = freq.time;
    label       = freq.label;

    % extract key info
    phase = single(angle(freq.fourierspctrm));
    
    % update user
    fprintf('\ncalculating phase bifurcation index...\n')

    % split phase data
    phase_hits = phase(trialinfo==1,:,:,:);
    phase_misses = phase(trialinfo==0,:,:,:);

    % get key variables
    ntimes = size(phase_hits,4);

    % get trials
    A = phase_hits;
    B = phase_misses;
    X = cat(1,A,B);

    % get PLV
    Ap = permute(abs(mean(exp(1i.*A))),[2 3 4 1]);
    Bp = permute(abs(mean(exp(1i.*B))),[2 3 4 1]);
    Xp = permute(abs(mean(exp(1i.*X))),[2 3 4 1]);

    % get phase bifurcation
    pbf = squeeze(mean((Ap - Xp) .* (Bp - Xp),3));
     
    % randomise conditions
    nperms = 1000;
    pbfc = zeros(nperms,size(pbf,2));
    for c = 1 : nperms

        % circshift data
        for trl = 1 : size(A,1); A(trl,:,:,:) = circshift(A(trl,:,:,:),randperm(ntimes,1),4); end
        for trl = 1 : size(B,1); B(trl,:,:,:) = circshift(B(trl,:,:,:),randperm(ntimes,1),4); end
        X = cat(1,A,B);

        % get PLV
        Ap = permute(abs(mean(exp(1i.*A))),[2 3 4 1]);
        Bp = permute(abs(mean(exp(1i.*B))),[2 3 4 1]);
        Xp = permute(abs(mean(exp(1i.*X))),[2 3 4 1]);

        % get phase bifurcation
        pbfc(c,:) = squeeze(mean((Ap - Xp) .* (Bp - Xp),3));

        % update user
        if mod(c,nperms/10)==0
            fprintf('%d%% complete...\n',round((c/nperms)*100))
        end
    end
    
    % find peak
    pbfz = (pbf-mean(pbfc))./std(pbfc);
    [~,idx] = max(pbfz);
    
    % plot
    figure;
    histogram(pbfc(:,idx),(-4:0.5:8).*0.001)
    xline(pbf(idx))
    xlim([-4 8]*0.001)
    
    % store data    
    vals(pp,:) = cat(1,pbf(idx),pbfc(:,idx));
    pf(pp,:) = zeros(nperms+1,1) + round(freqs(idx));
    obs(pp,:) = cat(1,1,zeros(nperms,1));
end

% save data
fig_1e_tbl = array2table(cat(1,vals,pf,obs(1,:))','VariableNames',...
             {'vals_pp1','vals_pp2','vals_pp3','vals_pp4','vals_pp5','vals_pp6',...
             'peakfreq_pp1','peakfreq_pp2','peakfreq_pp3','peakfreq_pp4','peakfreq_pp5','peakfreq_pp6',...
             'obs_boolean'});

%% Compute High-Res PBI TFR for mPFC
% get sourcemodel
[atlas,sourcemodel] = get_sourcemodel('frontal');

% cycle through participants
group_pbi = []; trl_count = [];
for pp = 1 : 6
                
    % define dataset
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,pp);
    filename = sprintf('sub-%02.0f_ephys-source.mat',pp);
    
    % load data
    load([directory,filename],'ephys')

    % use all channels
    label = ephys.label(atlas.mpfc_mask(sourcemodel.inside)>1);

    % select trials
    ephys = ft_selectdata(struct('channel',{label}),ephys);
    
    % get trialinfo
    trialinfo = ephys.trialinfo(:,3);
    
    % store trial counts
    trl_count(pp,:) = [sum(trialinfo==1) sum(trialinfo==0)];
    
    % get phase of ROI channels
    cfg         = [];
    cfg.output  = 'fourier';
    cfg.method  = 'wavelet';
    cfg.width   = 6;
    cfg.toi     = -0.8:0.025:0.8;
    cfg.foi     = 5:0.5:20;
    cfg.pad     = 'nextpow2';
    freq        = ft_freqanalysis(cfg,ephys); clear ephys
    freqs       = freq.freq;
    times       = freq.time;
    label       = freq.label;
    trialinfo   = freq.trialinfo(:,3);
    
    % update user
    fprintf('\ncalculating phase bifurcation index...\n')

    % extract key info
    phase = single(angle(freq.fourierspctrm));
    
    % split phase data
    phase_hits = phase(trialinfo==1,:,:,:);
    phase_misses = phase(trialinfo==0,:,:,:);

    % get key variables
    ntimes = size(phase_hits,4);

    % get trials
    A = phase_hits;
    B = phase_misses;
    X = cat(1,A,B);

    % get PLV
    Ap = permute(abs(mean(exp(1i.*A))),[2 3 4 1]);
    Bp = permute(abs(mean(exp(1i.*B))),[2 3 4 1]);
    Xp = permute(abs(mean(exp(1i.*X))),[2 3 4 1]);

    % get phase bifurcation
    pbf = squeeze(mean((Ap - Xp) .* (Bp - Xp))); % avg over chans
  
    % store group    
    group_pbi(:,:,pp) = pbf;
end

% cycle through participants
hlth_pbi = []; hlth_count = [];
for hlth = 1 : 12
                    
    % load data
    directory = sprintf('F:/meg-ieeg_magdeburg/derivatives/hlth-%02.0f/ephys/',hlth);%sprintf('%s/derivatives/hlth-%02.0f/ephys/',dir_data,hlth);
    filename = sprintf('%s/hlth-%02.0f_ephys-sourceMAG_5Hz.mat',directory,hlth);
    load(filename,'ephys')

    % use all channels
    label = ephys.label(atlas.frontal_mask(sourcemodel.inside)>1);

    % select trials
    ephys = ft_selectdata(struct('channel',{label}),ephys);
    
    % filter
    %ephys = ft_preprocessing(struct('hpfilter','yes','hpfreq',5),ephys);
    
    % get trialinfo
    trialinfo = ephys.trialinfo(:,4);
    hlth_count(hlth,:) = [sum(trialinfo==1) sum(trialinfo==0)];
    
    % get phase of ROI channels
    cfg         = [];
    cfg.output  = 'fourier';
    cfg.method  = 'wavelet';
    cfg.width   = 6;
    cfg.toi     = -0.8:0.025:0.8;    
    cfg.foi     = 5:0.5:20;
    cfg.pad     = 'nextpow2';
    freq        = ft_freqanalysis(cfg,ephys); clear ephys
    freqs       = freq.freq;
    times       = freq.time;
    label       = freq.label;
    trialinfo   = freq.trialinfo(:,4);
    
    % update user
    fprintf('\ncalculating phase bifurcation index...\n')

    % extract key info
    phase = single(angle(freq.fourierspctrm));
    
    % split phase data
    phase_hits = phase(trialinfo==1,:,:,:);
    phase_misses = phase(trialinfo==0,:,:,:);

    % get key variables
    ntimes = size(phase_hits,4);

    % get trials
    A = phase_hits;
    B = phase_misses;
    X = cat(1,A,B);

    % get PLV
    Ap = permute(abs(mean(exp(1i.*A))),[2 3 4 1]);
    Bp = permute(abs(mean(exp(1i.*B))),[2 3 4 1]);
    Xp = permute(abs(mean(exp(1i.*X))),[2 3 4 1]);

    % get phase bifurcation
    pbf = squeeze(mean((Ap - Xp) .* (Bp - Xp))); % avg over chans

    % store group
    hlth_pbi(:,:,hlth) = pbf;
end

% update user
fprintf('\ndone...\n\n')

% store source data
xls_data = {};
vals  = interp2(nanmean(group_pbi,3),4);
times = linspace(times(1),times(end),size(vals,2));
freqs = linspace(freqs(1),freqs(end),size(vals,1));
xls_data{1} = tfr_to_longform(vals,freqs,times);

% write xls data (for fig. 1f)
fig_1f_tbl = table(xls_data{1}(:,1),xls_data{1}(:,2),xls_data{1}(:,3),...
    'VariableNames',{'mpfc_pbi','freq','time'});

% store control data
vals  = interp2(nanmean(hlth_pbi,3),4);
times = linspace(times(1),times(end),size(vals,2));
freqs = linspace(freqs(1),freqs(end),size(vals,1));
xls_data{2} = tfr_to_longform(vals,freqs,times);

% write xls data (for fig. 1f)
fig_1g_tbl = table(xls_data{2}(:,1),xls_data{2}(:,2),xls_data{2}(:,3),...
    'VariableNames',{'control_pbi','freq','time'});

%% Compute High-Res PBI TFR w/o ERP
% cycle through participants
group_md = [];
for pp = 1 : 6
                
    % define dataset
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,pp);
    filename = sprintf('sub-%02.0f_ephys-bipolar.mat',pp);
    
    % load data
    load([directory,filename],'ephys')
    
    % record trialinfo
    trialinfo = ephys.trialinfo(:,3);
    
    % subtract erp
    tml = ft_timelockanalysis([],ephys);
    for trl = 1 : numel(ephys.trial)
        ephys.trial{trl} = ephys.trial{trl} - tml.avg;
    end
    
    % get phase of ROI channels
    cfg         = [];
    cfg.output  = 'fourier';
    cfg.method  = 'wavelet';
    cfg.width   = 6;
    cfg.toi     = -0.8:0.025:0.8;
    cfg.foi     = 5:0.5:20;
    cfg.pad     = 'nextpow2';
    cfg.channel = 'MD';
    freq        = ft_freqanalysis(cfg,ephys); clear ephys
    freqs       = freq.freq;
    times       = freq.time;
    label       = freq.label;

    % extract key info
    phase = single(angle(freq.fourierspctrm));
    
    % update user
    fprintf('\ncalculating phase bifurcation index...\n')

    % split phase data
    phase_hits = phase(trialinfo==1,:,:,:);
    phase_misses = phase(trialinfo==0,:,:,:);

    % get key variables
    ntimes = size(phase_hits,4);

    % get trials
    A = phase_hits;
    B = phase_misses;
    X = cat(1,A,B);

    % get PLV
    Ap = permute(abs(mean(exp(1i.*A))),[2 3 4 1]);
    Bp = permute(abs(mean(exp(1i.*B))),[2 3 4 1]);
    Xp = permute(abs(mean(exp(1i.*X))),[2 3 4 1]);

    % get phase bifurcation
    pbf = (Ap - Xp) .* (Bp - Xp);
    
    % store group
    group_md(:,:,pp) = pbf;
end

% get sourcemodel
[atlas,sourcemodel] = get_sourcemodel('mpfc');

% cycle through participants
group_pbi = []; trl_count = []; group_erp = [];
for pp = 1 : 6
                
    % define dataset
    directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',dir_data,pp);
    filename = sprintf('sub-%02.0f_ephys-source.mat',pp);
    
    % load data
    load([directory,filename],'ephys')

    % use all channels
    label = ephys.label(atlas.mpfc_mask(sourcemodel.inside)>1);

    % select trials
    ephys = ft_selectdata(struct('channel',{label}),ephys);
    
    % get trialinfo
    trialinfo = ephys.trialinfo(:,3);
      
    % subtract erp
    tml = ft_timelockanalysis([],ephys);
    group_erp(pp,:,1) = mean(tml.avg);
    tml = ft_timelockanalysis(struct('trials',find(trialinfo==1)),ephys);
    %tml = ft_timelockbaseline(struct('baseline',[-0.8 0]),tml);
    group_erp(pp,:,2) = mean(tml.avg);
    tml = ft_timelockanalysis(struct('trials',find(trialinfo==0)),ephys);
    %tml = ft_timelockbaseline(struct('baseline',[-0.8 0]),tml);
    group_erp(pp,:,3) = mean(tml.avg);
    for trl = 1 : numel(ephys.trial)
        if trialinfo(trl,1) == 1
            ephys.trial{trl} = ephys.trial{trl} - group_erp(pp,:,2);
        elseif trialinfo(trl,1) == 0
            ephys.trial{trl} = ephys.trial{trl} - group_erp(pp,:,3);
        end
    end
    
    % store trial counts
    trl_count(pp,:) = [sum(trialinfo==1) sum(trialinfo==0)];
    
    % get phase of ROI channels
    cfg         = [];
    cfg.output  = 'fourier';
    cfg.method  = 'wavelet';
    cfg.width   = 6;
    cfg.toi     = -0.8:0.025:0.8;
    cfg.foi     = 5:0.5:20;
    cfg.pad     = 'nextpow2';
    freq        = ft_freqanalysis(cfg,ephys); clear ephys
    freqs       = freq.freq;
    times       = freq.time;
    label       = freq.label;
    trialinfo   = freq.trialinfo(:,3);
    
    % update user
    fprintf('\ncalculating phase bifurcation index...\n')

    % extract key info
    phase = single(angle(freq.fourierspctrm));
    
    % split phase data
    phase_hits = phase(trialinfo==1,:,:,:);
    phase_misses = phase(trialinfo==0,:,:,:);

    % get key variables
    ntimes = size(phase_hits,4);

    % get trials
    A = phase_hits;
    B = phase_misses;
    X = cat(1,A,B);

    % get PLV
    Ap = permute(abs(mean(exp(1i.*A))),[2 3 4 1]);
    Bp = permute(abs(mean(exp(1i.*B))),[2 3 4 1]);
    Xp = permute(abs(mean(exp(1i.*X))),[2 3 4 1]);

    % get phase bifurcation
    pbf = squeeze(mean((Ap - Xp) .* (Bp - Xp))); % avg over chans
  
    % store group    
    group_pbi(:,:,pp) = pbf;
end

% cycle through participants
hlth_pbi = []; hlth_count = []; hlth_erp = [];
for hlth = 1 : 12
                    
    % load data
    directory = sprintf('%s/derivatives/hlth-%02.0f/ephys/',dir_data,hlth);
    filename = sprintf('%s/hlth-%02.0f_ephys-sourceMAG.mat',directory,hlth);
    load(filename,'ephys')

    % use all channels
    label = ephys.label(atlas.mpfc_mask(sourcemodel.inside)>1);

    % select trials
    ephys = ft_selectdata(struct('channel',{label}),ephys);
    
    % filter
    ephys_filt = ft_preprocessing(struct('lpfilter','yes','lpfreq',5),ephys);
    
    % get trialinfo
    trialinfo = ephys.trialinfo(:,4);
    hlth_count(hlth,:) = [sum(trialinfo==1) sum(trialinfo==0)];
     
    % subtract erp
    tml = ft_timelockanalysis([],ephys);
    hlth_erp(hlth,:,1) = mean(tml.avg);
    tml = ft_timelockanalysis(struct('trials',find(trialinfo==1)),ephys_filt);
    %tml = ft_timelockbaseline(struct('baseline',[-0.8 0]),tml);
    hlth_erp(hlth,:,2) = mean(tml.avg);
    tml = ft_timelockanalysis(struct('trials',find(trialinfo==0)),ephys_filt);
    %tml = ft_timelockbaseline(struct('baseline',[-0.8 0]),tml);
    hlth_erp(hlth,:,3) = mean(tml.avg);
    for trl = 1 : numel(ephys.trial)
        if trialinfo(trl,1) == 1
            ephys.trial{trl} = ephys.trial{trl} - group_erp(pp,:,2);
        elseif trialinfo(trl,1) == 0
            ephys.trial{trl} = ephys.trial{trl} - group_erp(pp,:,3);
        end
    end
    
    % get phase of ROI channels
    cfg         = [];
    cfg.output  = 'fourier';
    cfg.method  = 'wavelet';
    cfg.width   = 6;
    cfg.toi     = -0.8:0.025:0.8;    
    cfg.foi     = 5:0.5:20;
    cfg.pad     = 'nextpow2';
    freq        = ft_freqanalysis(cfg,ephys); clear ephys
    freqs       = freq.freq;
    times       = freq.time;
    label       = freq.label;
    trialinfo   = freq.trialinfo(:,4);
    
    % update user
    fprintf('\ncalculating phase bifurcation index...\n')

    % extract key info
    phase = single(angle(freq.fourierspctrm));
    
    % split phase data
    phase_hits = phase(trialinfo==1,:,:,:);
    phase_misses = phase(trialinfo==0,:,:,:);

    % get key variables
    ntimes = size(phase_hits,4);

    % get trials
    A = phase_hits;
    B = phase_misses;
    X = cat(1,A,B);

    % get PLV
    Ap = permute(abs(mean(exp(1i.*A))),[2 3 4 1]);
    Bp = permute(abs(mean(exp(1i.*B))),[2 3 4 1]);
    Xp = permute(abs(mean(exp(1i.*X))),[2 3 4 1]);

    % get phase bifurcation
    pbf = squeeze(mean((Ap - Xp) .* (Bp - Xp))); % avg over chans

    % store group
    hlth_pbi(:,:,hlth) = pbf;
end

% update user
fprintf('\ndone...\n\n')

% concatenate data w/ erp
mdX = table2array(fig_1c_tbl);
X = cat(1,mdX(:,[1 4 5]),table2array(fig_1f_tbl),table2array(fig_1g_tbl));

% store patient MD data
vals  = interp2(nanmean(group_md-nanmean(group_md(:)),3),4);
times = linspace(times(1),times(end),size(vals,2));
freqs = linspace(freqs(1),freqs(end),size(vals,1));
X(end+1:end+numel(vals),:) = tfr_to_longform(vals,freqs,times);

% store patient source data
vals  = interp2(nanmean(group_pbi-nanmean(group_pbi(:)),3),4);
X(end+1:end+numel(vals),:) = tfr_to_longform(vals,freqs,times);

% store control source data
vals  = interp2(nanmean(hlth_pbi,3),4);
X(end+1:end+numel(vals),:) = tfr_to_longform(vals,freqs,times);

% add additional columns
X(:,4) = [ones(size(fig_1f_tbl,1)*3,1);zeros(size(fig_1f_tbl,1)*3,1)]; % w/ vs. w/o ERP
X(:,5) = [ones(size(fig_1f_tbl,1),1);ones(size(fig_1f_tbl,1),1);zeros(size(fig_1f_tbl,1),1);...
          ones(size(fig_1f_tbl,1),1);ones(size(fig_1f_tbl,1),1);zeros(size(fig_1f_tbl,1),1)]; % patient vs. control ERP
X(:,6) = [ones(size(fig_1f_tbl,1),1);zeros(size(fig_1f_tbl,1),1);zeros(size(fig_1f_tbl,1),1);
          ones(size(fig_1f_tbl,1),1);zeros(size(fig_1f_tbl,1),1);zeros(size(fig_1f_tbl,1),1)]; % ieeg vs. source

% write xls data (for fig. 1f)
supp_6_tbl = array2table(X,'VariableNames',{'vals','freq','time','hasERP','isPatient','isiEEG'});

%% Save Data
writetable(fig_1c_tbl,[dir_repo,'/source_data/fig_1.xlsx'],'Sheet','figure_1c')
writetable(fig_1d_tbl,[dir_repo,'/source_data/fig_1.xlsx'],'Sheet','figure_1d')
writetable(fig_1e_tbl,[dir_repo,'/source_data/fig_1.xlsx'],'Sheet','figure_1e')
writetable(fig_1f_tbl,[dir_repo,'/source_data/fig_1.xlsx'],'Sheet','figure_1f')
writetable(fig_1g_tbl,[dir_repo,'/source_data/fig_1.xlsx'],'Sheet','figure_1g')

writetable(supp_6_tbl,[dir_repo,'/source_data/supp_fig_6.csv']) % save as csv because it's massive
