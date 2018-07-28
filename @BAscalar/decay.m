function coefDecay = decay(obj)
% returns norm of last Taylor coefficient collapsed onto first variable.


% Written S. Kepley 06/2017	
% Updated for BAscalar coefficient type 07/2017 

% TODO: Add support for measuring decay in other dimensions.

% ---------------------- INPUT ----------------------
% obj is a BAscalar
% varargin is dimension in which to measure the decay

% ---------------------- OUTPUT ----------------------
% decay of the coefficients in the required dimension

if numel(obj) > 1 % obj is a vector of BAscalars
    rowDecay = arrayfun(@(j)obj(j).decay,1:numel(obj));
    coefDecay = reshape(rowDecay,size(obj));
elseif strcmp(obj.CoefType,'BAscalar')
    coefDecay = obj.Coef(end).norm();
else
    switch obj.SurfaceDimension
        case 0
            coefDecay = abs(obj.Coef(1));
        case 1
            if strcmp(obj.Weight,'ones')
                coefDecay = abs(obj.Coef(end));
            else
                coefDecay = obj.Weight*abs(obj.Coef(end));
            end
        case 2
            if strcmp(obj.Weight,'ones')
                finalCoef = BAscalar(obj.Coef(end,:),obj.Modes(2:end));
                coefDecay = finalCoef.norm();
            else
                finalCoef = BAscalar(obj.Coef(end,:),obj.Modes(2:end),obj.Weight(2:end));
                coefDecay = finalCoef.norm();
            end
        case 3
            if strcmp(obj.Weight,'ones')
                finalCoef = BAscalar(obj.Coef(end,:,:),obj.Modes(2:end));
                coefDecay = finalCoef.norm();
            else
                finalCoef = BAscalar(obj.Coef(end,:,:),obj.Modes(2:end),obj.Weight(2:end));
                coefDecay = finalCoef.norm();
            end
    end
    if isa(obj.Coef,'intval')
        coefDecay = mid(coefDecay);
    end
end
end