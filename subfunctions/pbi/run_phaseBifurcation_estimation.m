function run_phaseBifurcation_estimation(filename)
  
% load data
load(filename,'ephys')

% define channels
if numel(ephys.label)>6
    
    % get sourcemodel    
    [atlas,sourcemodel] = get_sourcemodel;
    
    % get atlas channels
    label = ephys.label(atlas.inside(sourcemodel.inside));
else
    
    % use all channels
    label = ephys.label;
end

% if healthy, add filter to match patients
if strncmpi(filename,'hlth',4)
    ephys = ft_preprocessing(struct('hpfilter','yes','hpfreq',5),ephys);
end

% get phase using variable width wavelets
[phase,trialinfo,label,freqs,times] = get_tfr(ephys,label);

% if healthy, change trialinfo
if strncmpi(filename,'hlth',4)
    trialinfo = ephys.trialinfo(:,4);
end

% if source channels exists 
if numel(label)>10
    
    % get phase bifurcation for mediodorsal thalamus
    [pbi_source,pbi_raw,pbi_dist] = get_phaseBifurcation(phase,trialinfo,1:numel(label));

    % create frequency structure
    freq = struct('label',{label},'freq',freqs,'time',times,...
                  'dimord','chan_freq_time','cfg',[],'pbi',pbi_source,'pbi_raw',pbi_raw,'pbi_dist',mean(pbi_dist,4));

    % save data
    save([filename(1:end-17),'_phaseBif-source.mat'],'freq') 
end

% if mediodorsal thalamus contact exists 
if any(ismember(ephys.label,'MD'))
    
    % get indices of MD
    idx = ismember(label,'MD');
    
    % get phase bifurcation for mediodorsal thalamus
    [pbi_MD,pbi_raw,pbi_dist] = get_phaseBifurcation(phase,trialinfo,idx);

    % create frequency structure
    freq = struct('label',{label(idx)},'freq',freqs,'time',times,...
                  'dimord','chan_freq_time','cfg',[],'pbi',pbi_MD,...
                  'pbi_raw',pbi_raw,'pbi_dist',mean(pbi_dist,4));

    % save data
    save([filename(1:end-17),'phaseBif-MD.mat'],'freq') 
end

% if anterior thalamus contact exists 
if any(ismember(ephys.label,'ANT'))
    
    % get phase bifurcation for anterior thalamus
    [pbi_ANT,pbi_raw,pbi_dist] = get_phaseBifurcation(phase,trialinfo,ismember(label,'ANT')); 
        
    % create freq structure
    freq = struct('label',{{'ANT'}},'freq',freqs,'time',times,...
                  'dimord','chan_freq_time','cfg',[],'pbi',pbi_ANT,...
                  'pbi_raw',pbi_raw,'pbi_dist',mean(pbi_dist,4));
    
    % save data
    save([filename(1:end-17),'phaseBif-ANT.mat'],'freq') 
end    
end
