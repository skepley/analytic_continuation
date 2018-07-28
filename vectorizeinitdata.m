function ddt = vectorizeinitdata(t,y,varargin)
% vectorize initial conditions for the Lorenz system. 

% Written by S.K. 11/2015
% Implemented support for other parameters 03/2018

switch nargin
	case 2
		params = [28,10,8/3];
	case 3
		params = varargin{1};
end

% Pass ddt to ode45
idx = 1:3:length(y);
idy = 2:3:length(y);
idz = 3:3:length(y);

X = y(idx);
Y = y(idy);
Z = y(idz);
ddt = nan(size(y));

xdot = params(2)*(Y-X);
ydot = X.*(params(1) - Z) - Y;
zdot = X.*Y - params(3)*Z;

ddt(idx) = xdot;
ddt(idy) = ydot;
ddt(idz) = zdot;
end

