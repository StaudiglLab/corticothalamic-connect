function [plv_all,plv_diff,plv_lag,plv_perm,plv_hits,plv_misses] = getPhaseLockingValue(phase,idx_beh)
    
% update user
fprintf('\ncalculating phase locking value...\n')

% --- get absolute plv --- %
% get source phase details
As = phase(:,1:end-1,:,:);

% get thalamic phase details
At = repmat(phase(:,end,:,:),[1 size(As,2) 1 1]);

% get observed plv (drop first dimension)
plv_obs = permute(abs(mean(exp(1i.*(As-At)))),[2 3 4 1]);

% predefine permuted plv
nperms = 100;
plv_perm = zeros([size(plv_obs),nperms]);

% cycle through permutations
for p = 1 : nperms

    % shuffle trials in thalamus data
    At = At(randperm(size(At,1)),:,:,:);
    
    % get permuted plv
    plv_perm(:,:,:,p) = permute(abs(mean(exp(1i.*(As-At)))),[2 3 4 1]);

    % update user
    if mod(p,nperms/10)==0
        fprintf('absolute estimate: %d%% complete...\n',round((p/nperms)*100))
    end
end

% get z-score of observed PLV
plv_all = (plv_obs - mean(plv_perm,4)) ./ std(plv_perm,[],4);
plv_perm(:,:,:,end+1) = plv_obs;
clear plv_obs A*

% --- get behavioural difference  in PLV --- %
if nargout >= 2
    
    % splvt phase data
    phase_hits = phase(idx_beh==1,:,:,:);
    phase_misses = phase(idx_beh==0,:,:,:);

    % get condition with smaller n
    n_hits = size(phase_hits,1);
    n_misses = size(phase_misses,1);
    min_trl = min([n_hits,n_misses]);
    
    % get phase locking for hits
    As = phase_hits(:,1:end-1,:,:);
    At = repmat(phase_hits(:,end,:,:),[1 size(As,2) 1 1]);
    
    % get phase locking for misses
    Bs = phase_misses(:,1:end-1,:,:);
    Bt = repmat(phase_misses(:,end,:,:),[1 size(Bs,2) 1 1]);
        
    % cycle through permutations
    for p = 1 : nperms
        
        % get subset indices
        idx_misses = randperm(size(phase_misses,1),min_trl);
        idx_hits = randperm(size(phase_hits,1),min_trl);
                
        % get phase result
        plv_misses(:,:,:,p) = permute(abs(mean(exp(1i.*(Bs(idx_misses,:,:,:)-Bt(idx_misses,:,:,:))))),[2 3 4 1]);
        plv_hits(:,:,:,p) = permute(abs(mean(exp(1i.*(As(idx_hits,:,:,:)-At(idx_hits,:,:,:))))),[2 3 4 1]);
        
        % get phase lag
        misses_angle = permute(angle(mean(exp(1i.*(Bs(idx_misses,:,:,:)-Bt(idx_misses,:,:,:))))),[2 3 4 1]);
        hits_angle = permute(angle(mean(exp(1i.*(As(idx_hits,:,:,:)-At(idx_hits,:,:,:))))),[2 3 4 1]);
        plv_lag(:,:,:,p) = circ_dist(hits_angle,misses_angle);
        
        % update user
        if mod(p,nperms/10)==0
            fprintf('contrast estimate: %d%% complete...\n',round((p/nperms)*100))
        end
    end
    
    % determine difference
    plv_diff = mean(plv_hits,4) - mean(plv_misses,4);
    plv_hits = mean(plv_hits,4);
    plv_misses = mean(plv_misses,4);
end



