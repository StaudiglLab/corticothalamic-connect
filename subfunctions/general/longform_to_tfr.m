function [vals,freqs,times] = longform_to_tfr(longformData,freqs,times)

% get number of times/freqs
nfreqs = numel(unique(freqs));
ntimes = numel(unique(times));

% predefine [vals] output
vals = zeros(nfreqs,ntimes);

% get unique times/freqs;
ut = unique(times);
uf = unique(freqs);

% cycle through each data point
for n = 1 : numel(longformData)
    
    % find 2d idx
    i = find(freqs(n)==uf);
    j = find(times(n)==ut);
    
    % add value to tfr matrix
    vals(i,j) = longformData(n); %#ok<FNDSB>
end

% update freqs/times to unique values
freqs = uf;
times = ut;