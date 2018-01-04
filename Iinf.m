function [inflection_idx]= Iinf(x)

[pk,idx] = findpeaks(abs(diff(x)));
inflection_idx = sum((x(idx)/max(x)));
end