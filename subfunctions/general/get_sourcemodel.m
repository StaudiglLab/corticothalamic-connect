function [atlas,sourcemodel] = get_sourcemodel(roi)
    
% load brodmann atlas and source template
atlas = ft_read_mri('whole_brain_mask.nii');
load('C:\Users\ra34fod\github\fieldtrip\template\sourcemodel\standard_sourcemodel3d10mm.mat','sourcemodel');

% interpolate atlas with sourcemodel
cfg                 = [];
cfg.parameter       = 'anatomy';
atlas               = ft_sourceinterpolate(cfg,atlas,sourcemodel);

% find values for any BA region
atlas.anatomy = atlas.anatomy > 0;
atlas.inside  = atlas.anatomy == 1;

% get roi if requested
if nargin == 1;
    
    % load AAL
    aal = ft_read_mri('C:\Users\ra34fod\github\fieldtrip\template\atlas\aal\ROI_MNI_V4.nii');
    
    % interpolate aal with sourcemodel
    cfg                 = [];
    cfg.parameter       = 'anatomy';
    cfg.interpmethod    = 'nearest';
    aal                 = ft_sourceinterpolate(cfg,aal,sourcemodel);

    % switch based on requested roi
    switch roi
        case 'occipital'
            aal.anatomy(isnan(aal.anatomy)) = 0;
            aal.anatomy(aal.anatomy<5001 | aal.anatomy>5302) = 0;
        case 'occipital_ipsi'
            aal.anatomy(isnan(aal.anatomy)) = 0;
            aal.anatomy(aal.anatomy~=5001) = 0;
        case 'mpfc'
            aal.anatomy(aal.anatomy~=2111 & aal.anatomy~=2112 & aal.anatomy~=2702 & aal.anatomy~=2701) = 0;
        otherwise
            warning('mask "%s" does not exist...',roi)
    end
    
    % add to atlas
    atlas.([roi,'_mask']) = aal.anatomy;
end
