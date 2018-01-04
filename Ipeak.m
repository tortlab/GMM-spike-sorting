function [peak_idx]= Ipeak(x)

peak_idx  = sum(findpeaks(x/max(x)));

end