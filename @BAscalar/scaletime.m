function scaleTime(obj,L)
% rescales time

% Written by S. Kepley 05/2017
% Updated for BAscalar coefs 07/2017

% ---------------------- INPUT ----------------------
% obj (BAscalar): A Taylor series with 1st coordinate (time) scaled to t = 1
% L (double): The new time rescaling 

% ---------------------- OUTPUT ----------------------
% obj (BAscalar): The rescaled object

if numel(obj) > 1 % vectorized norm
    %                 scale_matr = BAscalar(repmat(bsxfun(@power,abs(L),[0:obj(1).Modes(1)-1]'),[1,obj(1).Modes(2)]));
    for j = 1:length(obj)
        obj(j).scaleTime(L);
    end
elseif obj.Modes(1) ==1
    error('obj.Modes should not equal 1') % rule out constant with respect to time
elseif isa(obj.Coef,'BAscalar')
    for k = 1:obj.Modes(1)
        obj.Coef(k).scaleTime(L^(k-1));
    end
else
    switch obj.SurfaceDimension
        case 1
            obj.Coef = abs(L)*obj.Coef;
        case 2
            obj.Coef = repmat(bsxfun(@power,abs(L),[0:obj.Modes(1)-1]'),[1,obj.Modes(2)]).*obj.Coef;
        case 3
            scaleCoef = @(j)L^(j-1).*obj.Coef(j,:,:);
            for j = 1:obj.Modes(1)
                obj.Coef(j,:,:) = scaleCoef(j);
            end
        otherwise
            error('Not implemented for higher dimension')
    end
end
end