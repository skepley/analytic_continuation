function invObj = inv(obj,varargin)
% Compute b = 1/A for specified A by solving the differential equation u' = -A'u^2 

% TODO: Add validation and make the argument recursive.

% Written by S. Kepley 05/2017
% Updated surfaceDimension 06/2017
% Support for BAscalar coef 07/2017

if strcmp(obj.CoefType,'BAscalar')
    const = obj.Coef(1).Coef(1);
else
    const = obj.Coef(1);
end

if abs(const) < 1e-13 % check if constant term is near 0. 
    disp('BAscalar is not invertible. Inverse should not be trusted')
end

switch obj.SurfaceDimension
    case 0 % inverse of a real scalar
        a0 = obj.Coef(1);
        invObj = BAscalar(1/a0);
        
    case 1 % inverse of 1d power series
        ddtObj = obj.ds;
        a0 = obj.Coef(1);
        invObj = BAscalar(1/a0,obj.Modes);
        for m = 1:obj.Modes - 1
            uu = invObj*invObj;
            Aprime = ddtObj.Coef(m:-1:1);
            newCoef = -(1/m)*dot(uu.Coef(1:m),Aprime);
            invObj.Coef(m+1) = newCoef;
        end

    case 2 % inverse of 2d power series
        deg = size(obj.Coef);
        ddtObj = obj.dt;
        ddtObj.Coef = [ddtObj.Coef(2:end,:);zeros(1,deg(2))]; % This padding shouldn't be necessary.
        a0 = BAscalar(obj.Coef(1,:));
        invObj = BAscalar(inv(a0),[1,a0.Modes]);
        for m = 1:obj.Modes(1)-1
            uu = invObj*invObj;
            Aprime = BAscalar(ddtObj.Coef(1:m,:),uu.Modes);
            newCoef = -(1/m)*mtimes(uu,Aprime,'Recursion');
            invObj.append(newCoef);
        end
    otherwise
        error('not yet implemented')
end
end

