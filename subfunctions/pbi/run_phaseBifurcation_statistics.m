function run_phaseBifurcation_statistics(directory,roi_str,participantType)
  
% define number of patients and controls
if strcmpi(roi_str,'ANT'); npp = 5; else; npp = 6; end
if nargin<3; participantType = 'patient'; end

% predefine cell for group data
group_data = {};
    
% cycle through participants 
switch participantType
    
    % get files for patients
    case 'patient'
        for pp = 1 : npp

            % load data
            pp_directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',directory,pp);
            filename = sprintf('sub-%02.0f_phaseBif-%s.mat',pp,roi_str);
            load([pp_directory,filename],'freq')

            % add data to group
            group_data{end+1} = freq;
            clear freq
        end

    % get files for patients
    case 'control'
        
        % cycle through participants 
        for hlth = 1 : 12

            % load mag data
            pp_directory = sprintf('%s/derivatives/hlth-%02.0f/ephys/',directory,hlth);
            filename = sprintf('hlth-%02.0f_phaseBif-%s.mat',hlth,roi_str);
            load([pp_directory,filename],'freq')

            % add data to group
            group_data{end+1} = freq;
            clear freq
        end   
end

% get grand average of participants
cfg                 = [];
cfg.keepindividual  = 'yes';
cfg.parameter       = {'pbi','pbi_raw','pbi_dist'};
freq                = ft_freqgrandaverage(cfg,group_data{:});

% if source data, convert data type to source
if strcmpi(roi_str,'source')
    
    % get sourcemodel
    [atlas,sourcemodel] = get_sourcemodel();
    
    % create source structure
    source = struct('dim',sourcemodel.dim,'inside',atlas.inside(:),'pos',sourcemodel.pos.*10,...
                    'label',{freq.label},'pbi',zeros(size(freq.pbi,1),numel(atlas.inside)),...
                    'unit','mm','dimord','rpt_pos','transform',atlas.transform);
                
    % add in measure
    if strcmpi(participantType,'patient')
        source.pbi(:,source.inside) = mean(mean(freq.pbi(:,:,freq.freq<=14,freq.time<0),4),3);
    elseif strcmpi(participantType,'control') % restrict to patient cluster
        source.pbi(:,source.inside) = mean(mean(freq.pbi(:,:,freq.freq>=11&freq.freq<=13,freq.time<0),4),3);
    end
    
    % rename source to fit code
    freq = source;
end

% set random seed
rng(1)

% define null hyp
null_hyp     = freq;
null_hyp.pbi = zeros(size(null_hyp.pbi));

% define stat design
design      = zeros(2,size(freq.pbi,1)*2);
design(1,:) = repmat(1:size(freq.pbi,1),[1 2]);
design(2,:) = [ones(1,size(freq.pbi,1)),ones(1,size(freq.pbi,1))+1];

% define general stat config structure
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 'all';
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.parameter           = 'pbi';
cfg.design              = design;
cfg.statistic           = 'ft_statfun_depsamplesT';  
cfg.tail                = 1;

% if is source
if strcmpi(roi_str,'source')
    
    % run stats
    stat                    = ft_sourcestatistics(cfg,freq,null_hyp);
else
    % get time-frequency domain specific stats
    cfg.latency             = [-0.8 0];
    cfg.frequency           = [6 14];
    stat                    = ft_freqstatistics(cfg,freq,null_hyp);
    
end

% get bayes at peak
addpath('C:\Users\ra34fod\github\bayesFactor\')
pidx = find(max(stat.stat(:))==stat.stat(:));
if ~strcmpi(roi_str,'source')
    [~,y,z] = ind2sub(size(stat.stat),pidx);
    tidx = stat.time(z); fidx = stat.freq(y);
    stat.bayes = bf.ttest(squeeze(freq.pbi(:,:,knnsearch(freq.freq',fidx),knnsearch(freq.time',tidx))));
else
    stat.bayes = bf.ttest(squeeze(freq.pbi(:,pidx)));
end

% get effect size at peak
stat.dz = abs(stat.stat ./ sqrt(size(freq.pbi,1)));

% report result
fprintf('\n--- Statistics: "%s" for "%s" ---\nt(%d): %3.3f\np: %3.3f\nBF: %3.3f\nd: %3.3f\n\n',roi_str,participantType,size(freq.pbi,1)-1,stat.stat(pidx),stat.posclusters(1).prob,stat.bayes,stat.dz(pidx))

% save outcome
if strcmpi(participantType,'control')
    save(sprintf('%s/derivatives/group/stat_pbi-%s_controls.mat',directory,roi_str),'stat','freq','-v7.3');
else
    save(sprintf('%s/derivatives/group/stat_pbi-%s.mat',directory,roi_str),'stat','freq','-v7.3');
end

end

