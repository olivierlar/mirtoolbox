function c = center(x)
if isempty(x)
    c = [];
else
    if find(isnan(x))
        warning('NaNs detected. Centering does not handle NaNs and will return NaN.');
    end
    c = x - repmat(mean(x),[size(x,1),1,1]);
end