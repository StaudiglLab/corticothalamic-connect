function s = sem(x,dim)

if nargin == 1; dim = 1; end
s = std(x,[],dim) ./ sqrt(size(x,dim));