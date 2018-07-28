function objDims = dims(obj)
% returns the dimensions for the BAscalar
objDims = size(obj.Coef) > 1;
end