
% load layout
lay=ft_prepare_layout(struct('layout','4D248.lay'));
chan_count = zeros(numel(lay.label),1);

for i = [1 2 4 5 6 7]
    filename = sprintf('F:/meg-ieeg_magdeburg_OLD/derivatives/sub-%02.0f/ephys/sub-%02.0f_ephys-visrej3.mat',i,i);
    load(filename)
    
    % cycle through channels
    for j = 1 : numel(ephys.label)
        idx = strcmpi(lay.label,ephys.label{j});
        chan_count(idx) = chan_count(idx) + 1;
    end
    
    % get chan+trial removal percentage
    crp(i) = 248 - numel(ephys.label);
    trp(i) = 432 - numel(ephys.trial);
end

% create fake struct
tml = [];
tml.avg = chan_count(1:248);
tml.time = 1;
tml.dimord = 'chan_time';
tml.label = lay.label(1:248);

% inflate layout
lay.pos(:,2) = lay.pos(:,2)-0.02;
lay.pos = lay.pos*1.05;

% save data
save('C:/Users/ra34fod/github/corticothalamic-connect/source_data/layout.mat','lay')
tbl = table(tml.label,tml.avg,'VariableNames',{'MEG_Label','participant_count'});
writetable(tbl,'C:/Users/ra34fod/github/corticothalamic-connect/source_data/supp_fig_6.xlsx');

% plot result
cmap = brewermap(7,'OrRd');
cfg = [];
cfg.layout = layb;
cfg.zlim = [-0.1 5.9];
cfg.style = 'fill';
cfg.gridscale = 400;
cfg.contournum = 6;
cfg.colormap = cmap;
cfg.interpolation = 'nearest';
cfg.comment = 'no';
ft_topoplotER(cfg,tml);

imagesc();
colormap(cmap);