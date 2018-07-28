function append(obj,nextCoef,varargin)
% Appends a new coefficient to given dimension (default is 1st dimension). 

% Written by S. Kepley 05/9/17 
% Updated for multiple surface dimensions 06/24/17
% Updated for BAscalar coefficients 07/02/17

% ---------------------- INPUT ----------------------
% obj is a BAscalar
% nextCoef is a coefficient array of size [1,obj.SurfaceDimension-1]

% ---------------------- OUTPUT ----------------------
% obj (BAscalar): The rescaled object

if strcmp(obj.CoefType,'BAscalar') % append BAscalar coefficient
    obj.Coef(end+1) = nextCoef; 
    obj.Modes(1) = obj.Modes(1) + 1;
    
else % append double/intval array
    if isa(nextCoef,'BAscalar')
        nextCoef = nextCoef.Coef;
    end
    
    try
        switch obj.SurfaceDimension
            case 1
                obj.Coef(end+1) = nextCoef;
            case 2
                obj.Coef(end+1,:) = nextCoef;
            case 3
                obj.Coef(end+1,:,:) = nextCoef;
            otherwise
                error('Not yet implemented for higher dimension')
        end
        obj.Modes(1) = obj.Modes(1) + 1;
    catch
        error('Appended coefficient must have same dimension as the surface being appended to')
    end
end
end