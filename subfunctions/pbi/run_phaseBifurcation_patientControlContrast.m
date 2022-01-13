function run_phaseBifurcation_patientControlContrast(directory)
  
% predefine cell for group data
group_patient = {};
    
% get files for patients
for pp = 1 : 6

    % load data
    pp_directory = sprintf('%s/derivatives/sub-%02.0f/ephys/',directory,pp);
    filename = sprintf('sub-%02.0f_phaseBif-source.mat',pp);
    load([pp_directory,filename],'freq')

    % add data to group
    group_patient{end+1} = freq;
    clear freq
end

% predefine cell for group data
group_control = {};
    
% cycle through participants 
for hlth = 1 : 12

    % load mag data
    pp_directory = sprintf('%s/derivatives/hlth-%02.0f/ephys/',directory,hlth);
    filename = sprintf('hlth-%02.0f_phaseBif-source.mat',hlth);
    load([pp_directory,filename],'freq')

    % add data to group
    group_control{end+1} = freq;
    clear freq
end   

% get grand average of participants
cfg                 = [];
cfg.keepindividual  = 'yes';
cfg.parameter       = {'pbi'};
freq1                = ft_freqgrandaverage(cfg,group_patient{:});
freq2                = ft_freqgrandaverage(cfg,group_control{:});

% get sourcemodel
[atlas,sourcemodel] = get_sourcemodel();

% create source structure
source1 = struct('dim',sourcemodel.dim,'inside',atlas.inside(:),'pos',sourcemodel.pos.*10,...
                'label',{freq1.label},'pbi',zeros(size(freq1.pbi,1),numel(atlas.inside)),...
                'unit','mm','dimord','rpt_pos','transform',atlas.transform);
source2 = struct('dim',sourcemodel.dim,'inside',atlas.inside(:),'pos',sourcemodel.pos.*10,...
                'label',{freq2.label},'pbi',zeros(size(freq2.pbi,1),numel(atlas.inside)),...
                'unit','mm','dimord','rpt_pos','transform',atlas.transform);

% add in measure
source1.pbi(:,source1.inside) = mean(mean(freq1.pbi(:,:,freq1.freq<=14,freq1.time<0),4),3);
source2.pbi(:,source2.inside) = mean(mean(freq2.pbi(:,:,freq1.freq<=14,freq1.time<0),4),3);

% rename source to fit code
freq1 = source1;
freq2 = source2;

% set random seed
rng(1)

% define stat design
design      = zeros(2,size(freq1.pbi,1)+size(freq2.pbi,1));
design(1,:) = [1:size(freq1.pbi,1) 1:size(freq2.pbi,1)];
design(2,:) = [ones(1,size(freq1.pbi,1)) ones(1,size(freq2.pbi,1))+1];

% define general stat config structure
cfg                     = [];
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 2000;
cfg.ivar                = 2;
%cfg.uvar                = 1;
cfg.parameter           = 'pbi';
cfg.design              = design;
cfg.tail                = 1;
cfg.statistic           = 'ft_statfun_indepsamplesT';  
stat                    = ft_sourcestatistics(cfg,freq1,freq2);

% get bayes at peak
addpath('C:\Users\ra34fod\github\bayesFactor\')
pidx = find(max(stat.stat(:))==stat.stat(:));
stat.bayes = bf.ttest(squeeze(freq1.pbi(:,pidx))-squeeze(freq2.pbi(:,pidx)));

% get effect size at peak
stat.dz = abs(stat.stat ./ sqrt(size(freq1.pbi,1)));

% report result
fprintf('\n--- Statistics: "%s" for "%s" ---\nt(%d): %3.3f\np: %3.3f\nBF: %3.3f\nd: %3.3f\n\n',roi_str,participantType,size(freq1.pbi,1)-1,stat.stat(pidx),stat.posclusters(1).prob,stat.bayes,stat.dz(pidx))

% save outcome
if strcmpi(participantType,'control')
    save(sprintf('%s/derivatives/group/stat_pbi-%s_controls.mat',directory,roi_str),'stat','freq','-v7.3');
else
    save(sprintf('%s/derivatives/group/stat_pbi-%s.mat',directory,roi_str),'stat','freq','-v7.3');
end

% save source plot
if strcmpi(roi_str,'source')
    
    % get sourcemodel information
    [atlas,sourcemodel] = get_sourcemodel;

    % get group average
    cfg                 = [];
    cfg.keepindividual  = 'yes';
    cfg.parameter       = {'pbi','pbi_raw','pbi_dist'};
    freq1                = ft_freqgrandaverage(cfg,group_data{:});
    
    % create source structure
    source = struct('dim',sourcemodel.dim,'inside',atlas.inside(:),'pos',sourcemodel.pos.*10,...
                    'label',{freq1.label},'pbi',zeros(numel(atlas.inside),1),...
                    'unit','mm','dimord','rpt_pos','transform',atlas.transform);
                
    % add in measure
    if strcmpi(participantType,'patient')
        source.pbi(source.inside) = mean(mean(mean(freq1.pbi(:,:,freq1.freq>=10&freq1.freq<=12,freq1.time>=-0.5&freq1.time<=-0.3),4),3));
    elseif strcmpi(participantType,'control') % restrict to patient cluster
        source.pbi(source.inside) = mean(mean(mean(freq1.pbi(:,:,freq1.freq>=10&freq1.freq<=12,freq1.time>=-0.4&freq1.time<=-0.2),4),3));
    end
    
%     % get measure of interest
%     pow = mean(source.pbi,1);
% 
%     % create source data structure
%     source              = [];
%     source.dim          = sourcemodel.dim;
%     source.pos          = sourcemodel.pos;
%     source.inside       = atlas.inside;
%     source.pow          = zeros(numel(atlas.inside),1);
%     source.pow(source.inside) = pow(source.inside);
%     source.dimord       = 'rpt_pos';
%     source.pos          = source.pos.*10;
%     source.unit         = 'mm';

    % reshape source
%     source.pow = reshape(source.pow,source.dim);

    % add transformation matrix
    source.transform        = [1,0,0,-91;
                               0,1,0,-127;
                               0,0,1,-73;
                               0,0,0,1];
                       
    % export stat
    cfg                 = [];
    cfg.parameter       = 'pbi';               % specify the functional parameter to write to the .nii file
    cfg.filename        = ['C:/Users/ra34fod/github/corticothalamic-connect/source_data/pbi_',participantType,'.nii'];  % enter the desired file name
    cfg.filetype        = 'nifti';
    cfg.vmpversion      = 2;
    cfg.datatype        = 'float';
    ft_volumewrite(cfg,source);      % be sure to use your interpolated source data   
    reslice_nii(cfg.filename,cfg.filename,[1 1 1])

end

end

