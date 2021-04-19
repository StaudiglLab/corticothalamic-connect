function [pbfz,pbf,pbfc] = get_phaseBifurcation(phase,idx_beh,idx_chan)
    
% update user
fprintf('\ncalculating phase bifurcation index...\n')
    
% split phase data
phase_hits = phase(idx_beh==1,idx_chan,:,:);
phase_misses = phase(idx_beh==0,idx_chan,:,:);

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
    
% randomise conditions
nperms = 1000;
pbfc = zeros([size(pbf),nperms]);
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
    pbfc(:,:,:,c) = (Ap - Xp) .* (Bp - Xp);
        
    % update user
    if mod(c,nperms/10)==0
        fprintf('%d%% complete...\n',round((c/nperms)*100))
    end
end

% z-transform result
pbfz = (pbf - nanmean(pbfc,4)) ./ nanstd(pbfc,[],4);
