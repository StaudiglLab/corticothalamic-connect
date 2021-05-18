function stat = runMediationAnalysis(X,y)

% check outcome variable
if numel(unique(y)) == 2; useLogit = true;
else; useLogit = false;
end

% run IV -> MV GLM (always use linear here)
a = int_runModel(X(:,1),X(:,2),false);

% run IV -> DV logisitic
c = int_runModel(X(:,1),y,true);

% run full model
tmp = int_runModel(X,y,true);
b = tmp(2);
cdash = tmp(1);

% package result
stat.a = a;
stat.b = b;
stat.c = c;
stat.cdash = cdash;
stat.ab = (stat.a.*stat.b) ./ sqrt(stat.a.^2 + stat.b.^2 + 1);
stat.d = stat.c - stat.cdash;
 
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

