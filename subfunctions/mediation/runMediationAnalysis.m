function stat = runMediationAnalysis(X,y)

% check outcome variable
if numel(unique(y)) == 2; useLogit = true;
else; useLogit = false;
end

% cycle through betas
for beta = 1 : 2
    
    % run IV -> MV GLM (always use linear here)
    a = int_runModel(X(:,beta),X(:,~ismember([1 2],beta)),false);
    
    % run full model
    tmp = int_runModel(X(:,[beta find(~ismember([1 2],beta))]),y,useLogit);
    b = tmp(2);
    
    % package result
    stat{beta}.a = a;
    stat{beta}.b = b;
    stat{beta}.cdash = tmp(1);
    stat{beta}.ab = (stat{beta}.a.*stat{beta}.b) ./ sqrt(stat{beta}.a.^2 + stat{beta}.b.^2 + 1);
    stat{beta}.d = tmp(2) - tmp(1); % difference: indirect path > direct path
end
 
end

function tval = int_runModel(X,y,useLogit)

% run model
if useLogit; [~,~,stat] = mnrfit(X,y);
else; [~,~,stat] = glmfit(X,y);
end

% get betas
tval = stat.t(2:end);

% fix sign bug in mnrfit
if useLogit
    for i = 1 : size(X,2)
        [~,~,~,tstat] = ttest2(X(y==2,i),X(y==1,i));
        if sign(tval(i)) ~= sign(tstat.tstat)
            tval(i) = tval(i)*-1;
        end
    end
end
%             
% % y-standardise coefficents
% if numel(beta) > 1
%     a = (beta(1).^2).*var(X(:,1));
%     b = (beta(2).^2).*var(X(:,2));
%     cov_ab = cov(X(:,1),X(:,2));
%     ab = (2.*beta(1)*beta(2)).*cov_ab(2,1);
%     beta = beta ./ sqrt(a + b + ab + (pi.^2/3));
% else
%     beta = beta ./ sqrt((beta.^2).*var(X) + (pi.^2/3));
% end

end

