function run_connectivity_estimation(directory,subj,metric)
      
% if no metric supplied, use PLV
if nargin == 2; metric = 'PLV'; end

% define savefile name
fileout = sprintf('%ssub-%02.0f_pliconnect-sourceMD_%s.mat',directory,subj,metric);

% get phase data for source and MD data
data = internal_prepareData(directory,subj,'source','MD',metric);

% run source-MD connectivity
internal_connect(data,fileout,metric)

end 
    
function data_out = internal_prepareData(directory,subj,roiA,roiB,metric)
    
% update user
fprintf('\npreparing sub-%02.0f data for %s-%s connectivity...\n',subj,roiA,roiB)
    
% determine if any data is source
issource = strcmpi(roiA,'source') | strcmpi(roiB,'source');
    
% prepare filename   
if issource
    filename = sprintf('sub-%02.0f_ephys-%s%s.mat',subj,roiA,roiB);
    varname = 'source'; % define variable name of previously saved data
else
    filename = sprintf('sub-%02.0f_ephys-bipolar.mat',subj);
    varname = 'ephys'; % define variable name of previously saved data
end

% load data
fprintf('loading data...\n')
load([directory,filename],varname)    

% get phase data
if issource
    
    % get sourcemodel    
    [atlas,sourcemodel] = get_sourcemodel;

    % define channels
    chans = cat(1,source.label(atlas.inside(sourcemodel.inside)),{roiB});
    
    % get phase
    switch metric
        case 'imagCoh'; [phase,trialinfo,label,freqs,times,fspc] = get_tfr(source,chans);
        otherwise; [phase,trialinfo,label,freqs,times] = get_tfr(source,chans); 
    end
    
else
    chans = {'MD','ANT'};
    switch metric
        case 'imagCoh'; [phase,trialinfo,label,freqs,times,fspc] = get_tfr(source,chans);
        otherwise; [phase,trialinfo,label,freqs,times] = get_tfr(source,chans); 
    end 
    
end

% organise output
data_out.phase = phase;
data_out.trialinfo = trialinfo;
data_out.label = label;
data_out.freqs = freqs;
data_out.times = times;
if exist('fspc','var'); data_out.fourierspctrm = fspc; end

end

function internal_connect(data,save_name,metric)
 
% compute connectivity for every trial and channel
switch metric
    case 'PLV'; [connect_all,connect_beh,connect_lag] = getPhaseLockingValue(data.phase,data.trialinfo);
    case 'TimeResolvedPLV'; [connect_all,connect_hits,connect_misses] = getTimeResolvedPLV(data.phase,data.trialinfo,data.times);
    case 'imagCoh'; [connect_all,connect_beh] = getImagCoherence(data);
end

% store data
data = struct('freq',data.freqs,'time',data.times,'label',{data.label(1:end-1)},'cfg',[],...
              'dimord','chan_freq_time','all',connect_all);
          
% add additional data
if exist('connect_lag','var');      data.lag = connect_lag; end
if exist('connect_beh','var');      data.beh = connect_beh; end
if exist('connect_hits','var');     data.hits = connect_hits; end
if exist('connect_misses','var');   data.misses = connect_misses; end

% fix time for time-resolved PLV
if size(connect_all,3) == 2
    data.time = [-0.4 0.4];
end
          
% save data
fprintf('\nsaving data...\n')
save(save_name,'data')
fprintf('done...\n\n')

end
