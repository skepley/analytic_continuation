function varargout = sqrt(obj,varargin)
% compute square roots for BAscalars

% Written by S. Kepley 05/2017
% Updated to compute by automatic differentiation and added lower dimensions 06/2017


% ---------------------- INPUT ----------------------
% obj = A is a BAscalar
% varargin = 1 or -1 to specify which branch of sqrt. (default is 1)

% ---------------------- OUTPUT ----------------------
% nargout = 1: Output is x satisfying x^2 - A = 0 in the Banach algebra.
% nargout = 2: Output is [x,u] where u satisfies u*x = 1 in the Banach algebra.

if strcmp(obj.CoefType,'intval')
%     disp('sqrt is untested with interval valued BAscalars')
end

if abs(obj.Coef(1,1)) < 1e-9
    disp('constant term < 1e-9. square root should not be trusted')
end

if nargin ==2
    branch = varargin{1};
else
    branch = 1;
end

switch obj.SurfaceDimension
    case 0 % sqrt of a real scalar
        a0 = abs(obj.Coef(1));
        sqrtObj = BAscalar(branch*sqrt(a0));
        invObj = BAscalar(1/sqrtObj.Coef(1));
        
    case 1 % sqrt of 1d power series
        dObj = obj.ds;
        a0 = abs(obj.Coef(1));
        sqrtObj = BAscalar(branch*sqrt(a0),1,1);
        invObj = BAscalar(1/sqrtObj.Coef(1),1,1);
        for m = 1:obj.Modes - 1
            dA = BAscalar(dObj.Coef(1:m),invObj.Modes);
            dsqrt = dA*invObj; % obj'* 
            uu = invObj*invObj;
            sqrtNext = .5*(1/m)*dsqrt.Coef(m);
            invNext = .5*(-1/m)*mtimes(uu,dsqrt,'Recursion');
            sqrtObj.append(sqrtNext);
            invObj.append(invNext);
        end
    case 2 % sqrt of 2d power series
        dObj = obj.dt;
        dObj.Coef = dObj.Coef(2:end,:);
        a0 = BAscalar(obj.Coef(1,:));
        [sqrtInit,invInit] = sqrt(a0);
        sqrtObj = BAscalar(sqrtInit,[1,a0.Modes]);
        invObj = BAscalar(invInit,[1,a0.Modes]);
        for m = 1:obj.Modes(1)-1
            dA = BAscalar(dObj.Coef(1:m,:),invObj.Modes);
            dsqrt = dA*invObj;
            uu = invObj*invObj;
            sqrtNext = .5*(1/m)*dsqrt.Coef(m,:);
            invNext = .5*(-1/m)*mtimes(uu,dsqrt,'Recursion');
            sqrtObj.append(sqrtNext);
            invObj.append(invNext);
        end
    otherwise
        error('sqrt not implemented for this dimension')
end

if nargout == 1
    varargout{1} = sqrtObj;
elseif nargout ==2
    varargout{1} = sqrtObj;
    varargout{2} = invObj;
end
end





% Old version only for 2D and didn't use automatic differentiation.

% sqrt_ = zeros(obj.Modes);
% sqrt_(1,1) = sqrt(obj.Coef(1,1)); %initial condition
% D = 2*sqrt_(1,1); %scalar denominator
% for k = 1:obj.Modes(2)-1
% 	sqrt_(1,k+1) = (-1/D)*(recursive_conv(sqrt_(1,1:k+1),sqrt_(1,1:k+1),'1d') - obj.Coef(1,k+1));
% end
%
% D = 2*sqrt_(1,:); %row denominator
% for k = 1:obj.Modes(1)-1;
% 	Xnum = -(recursive_conv(sqrt_(1:k+1,:),sqrt_(1:k+1,:),'row') - obj.Coef(k+1,:));
% 	sqrt_(k+1,:) = taylor_divide(Xnum,D);
% end
% sqrtObj = BAscalar(sqrt_);
% end