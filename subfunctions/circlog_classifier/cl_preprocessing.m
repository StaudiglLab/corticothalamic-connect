function [phase,chan_labels] = cl_preprocessing(cfg,eeg)

% preprocessing Fieldtrip-formatted data and return sine + cosine of signal

% add missing fields
if ~isfield(cfg,'channel'); cfg.channel = 'all'; end
if ~isfield(cfg,'latency'); cfg.latency = 'all'; end
if ~isfield(cfg,'frequency'); cfg.frequency = 'all'; end
if ~isfield(cfg,'dimreduction'); cfg.dimreduction = 'none'; end
if ~isfield(cfg.dimreduction,'method') && ~strcmpi(cfg.dimreduction,'none'); cfg.dimreduction.method = 'pca'; end
if ~isfield(cfg.dimreduction,'nfeatures') && ~strcmpi(cfg.dimreduction,'none'); cfg.dimreduction.nfeatures = 20; end

% restrict to selected channels
eeg = ft_selectdata(struct('channel',cfg.channel),eeg);

% filter
if ~strcmpi(cfg.frequency,'all')
    eeg = ft_preprocessing(struct('bpfilter','yes','bpfreq',cfg.frequency),eeg);
end

% if dimenson reduction requested
if ~strcmpi(cfg.dimreduction,'none')

    % convert to components
    comps = ft_componentanalysis(struct('method',cfg.dimreduction.method,...
                                        'numcomponent',cfg.dimreduction.nfeatures),eeg);

    % extract data
    signal = reshape(cell2mat(comps.trial),[numel(comps.label),numel(comps.time{1}),numel(comps.time)]);  
    
else
    
    % extract data
    signal = reshape(cell2mat(ephys.trial),[numel(ephys.label),numel(ephys.time{1}),numel(ephys.time)]);  
end

% extract phase/power
fprintf('\ncalculating sine/cosine of signal...\n')
phase = zeros(size(signal,1)*2,size(signal,2),size(signal,3));
chan_labels = cell(size(signal,1)*2,1);
for trl = 1 : size(signal,3)
    for chan = 1 : size(signal,1)

        % calculate phase
        phase(chan,:,trl) = cos(angle(hilbert(signal(chan,:,trl))));
        phase(chan+size(signal,1),:,trl) = sin(angle(hilbert(signal(chan,:,trl))));
        
        % record channel label
        if trl == 1
            chan_labels{chan} = [comps.topolabel{chan},'_cosine'];
            chan_labels{chan+size(signal,1)} = [comps.topolabel{chan},'_sine'];
        end
    end
end

% restrict to time window of interest
if ~strcmpi(cfg.latency,'all')
    time_idx = comps.time{1}>=cfg.latency(1)&comps.time{1}<=cfg.latency(2);
    phase = phase(:,time_idx,:);
end

fprintf('cl_preprocessing complete...\n')
