function data = tfr_to_longform(tfr,freqs,times)

% duplicate dimensions of time and freq
freqs = repmat(freqs',[1 size(times,2)]);
times = repmat(times,[size(freqs,1) 1]);

% reshape vals matrix
data = cat(2,tfr(:),freqs(:),times(:));