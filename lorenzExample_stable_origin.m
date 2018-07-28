%% ----------------------------------------- LORENZ EXAMPLE -----------------------------------------
% Demonstration of rigorous analytic continuation of the (slow) stable manifold at the origin for the Lorenz system. Initial local manifold is computed with analytic error bound ~1e-14 and stored in 
% the file origin_local_data.mat. The local manifold boundary is parameterized by lifting 8 piecewise linear segments from the model space through the local parameterization, and then the global
% manifold is parameterized by rigorously integrating this boundary. The resulting atlas is plotted with the corresponding error. 

% EXAMPLE:
% compact box: [-100,100]x[-100,100]x[-40,120]
% tau: 1 time unit
% truncation: (M,N) = [40,40]
% initial subdivision: 8 arcs
% maximum time-step: .1
% subdivision threshold: 1 decimal point per time-step
% number of subarcs: 4

% RESULTS: 
% Maximum error: 1.2229e-6
% Number of charts: 2009
% Computation time: ~4 hours

load origin_local_data  % local manifold data 
localError = r_minus; 
parameter = [28,10,8/3];
 

%% ----------------------------------------- INTEGRATION PARAMETERS -----------------------------------------
% Define compact box to confine global manifold. Stop advection of boundary arcs which leave this box (only in the x,y directions). Choosing this takes some care as the stable manifold in this example
% will leave small boxes only to return a short time later. The regions which are removed may appear as gaps in the final manifold when surrounding portions have returned (and are thus 
% retained). 
boxCenter = [0,0,40];
boxRadius = [100,100,80];

stopTime = -ones(1,8); % Maximum integration time for boundary regions which do not exit the box. 

N = 40; % spatial truncation
M = 40; % temporal truncation
modes = [M,N]; 
% The choice of M,N are important for a good run and determining good choices is a difficult problem. If your error is growing rapidly or you can't integrate long enough, these 
% parameters are the most likely suspects. As a general (but not exact) rule:
% Increasing M --- longer timesteps, fewer charts, reduced error due to floating point roundoff (fewer validations), and increased error due to truncation. 
% Increasing N --- fewer subdivisions, tighter root bounds for radii polynomials, reduced truncation error. 
% Increasing either has a dramatic effect on runtime and efficiency. The "game" is to take M as large as possible, increasing N to match in order to control errors until the tradeoff 
% in computational efficiency is no longer worth it. This requires a fair bit of guesswork and heuristics on the user end. 

% Single timestep integration and subdivision
errTol = 1; % precision loss threshold (decimal places per time-step). Arcs with greater error propagation are subdivided.
maxTau = .1; % maximal tau for a single timestep
numSubDivs = 4; % number of manifold divisions to perform when needed

% shrinkwrap local parameterization to spatial truncation
[Ps1,r1Local] = shrinkwrap(P(:,:,1),[N,N],localError);
[Ps2,r2Local] = shrinkwrap(P(:,:,2),[N,N],localError);
[Ps3,r3Local] = shrinkwrap(P(:,:,3),[N,N],localError);
if max([r1Local,r2Local,r3Local]) > 1e-13
    disp('Local parameterization error is large - Increase N or rescale eigenvectors')
end
P = zeros(N,N,3);
P(:,:,1) = Ps1;
P(:,:,2) = Ps2; 
P(:,:,3) = Ps3;
rLocal = [r1Local,r2Local,r3Local];

%% ----------------------------------------- INTEGRATE BOUNDARIES -----------------------------------------
Arc = cell(1,8); % Initialize chart atlas with initial boundary division into 8 piecewise linear subarcs.

% SouthEast face of manifold
if stopTime(1) < 0
    timespan = [0,stopTime(1)];
    disp('Integrating SouthEast Face')
    
    %  lift parameterization
    segmentFrom = [-1;-1];
    segmentTo = [-1;0];
    [gamma,initError] = liftsegment(segmentFrom,segmentTo,P,localError);
    
    % integrate
    Arc{1} = lorenzbox([gamma(1,:);gamma(2,:);gamma(3,:)],initError,parameter,modes,timespan,errTol,numSubDivs,boxCenter,boxRadius);
end

% SouthWest face of manifold
if stopTime(2) < 0
    timespan = [0,stopTime(2)];
    disp('Integrating SouthWest Face')
    
    % lift parameterization
    segmentFrom = [-1;0];
    segmentTo = [-1;1];
    [gamma,initError] = liftsegment(segmentFrom,segmentTo,P,localError);
    
    % integrate
    Arc{2} = lorenzbox([gamma(1,:);gamma(2,:);gamma(3,:)],initError,parameter,modes,timespan,errTol,numSubDivs,boxCenter,boxRadius);
end

% WestSouth face of manifold
if stopTime(3) < 0
    timespan = [0,stopTime(3)];
    disp('Integrating WestSouth Face')
    
    % lift parameterization
    segmentFrom = [-1;1];
    segmentTo = [0;1];
    [gamma,initError] = liftsegment(segmentFrom,segmentTo,P,localError);
    
    % integrate
    Arc{3} = lorenzbox([gamma(1,:);gamma(2,:);gamma(3,:)],initError,parameter,modes,timespan,errTol,numSubDivs,boxCenter,boxRadius);
end

% WestNorth face of manifold
if stopTime(4) < 0
    timespan = [0,stopTime(4)];
    disp('Integrating WestNorth Face')
    
    % lift parameterization
    segmentFrom = [0;1];
    segmentTo = [1;1];
    [gamma,initError] = liftsegment(segmentFrom,segmentTo,P,localError);
    
    % integrate
    Arc{4} = lorenzbox([gamma(1,:);gamma(2,:);gamma(3,:)],initError,parameter,modes,timespan,errTol,numSubDivs,boxCenter,boxRadius);
end

% NorthWest face of manifold
if stopTime(5) < 0
    timespan = [0,stopTime(5)];
    disp('Integrating NorthWest Face')
    
    % lift parameterization
    segmentFrom = [1;1];
    segmentTo = [1;0]; % [-1;1]
    [gamma,initError] = liftsegment(segmentFrom,segmentTo,P,localError);
    
    % integrate
    Arc{5} = lorenzbox([gamma(1,:);gamma(2,:);gamma(3,:)],initError,parameter,modes,timespan,errTol,numSubDivs,boxCenter,boxRadius);
end

% NorthEast face of manifold
if stopTime(6) < 0
    timespan = [0,stopTime(6)];
    disp('Integrating NorthEast Face')
    
    % lift parameterization
    segmentFrom = [1;0];
    segmentTo = [1;-1];
    [gamma,initError] = liftsegment(segmentFrom,segmentTo,P,localError);
    
    % integrate
    Arc{6} = lorenzbox([gamma(1,:);gamma(2,:);gamma(3,:)],initError,parameter,modes,timespan,errTol,numSubDivs,boxCenter,boxRadius);
end

% EastNorth face of manifold
if stopTime(7) < 0
    timespan = [0,stopTime(7)];
    disp('Integrating EastNorth Face')
    
    % lift parameterization
    segmentFrom = [1;-1];
    segmentTo = [0;-1]; % [-1;1]
    [gamma,initError] = liftsegment(segmentFrom,segmentTo,P,localError);
    
    % integrate
    Arc{7} = lorenzbox([gamma(1,:);gamma(2,:);gamma(3,:)],initError,parameter,modes,timespan,errTol,numSubDivs,boxCenter,boxRadius);
end

% EastSouth face of manifold
if stopTime(8) < 0 
    timespan = [0,stopTime(8)];
    disp('Integrating EastSouth Face')
    
    % lift parameterization
    segmentFrom = [0;-1];
    segmentTo = [-1;-1];
    [gamma,initError] = liftsegment(segmentFrom,segmentTo,P,localError);
    
    % integrate
    Arc{8} = lorenzbox([gamma(1,:);gamma(2,:);gamma(3,:)],initError,parameter,modes,timespan,errTol,numSubDivs,boxCenter,boxRadius);
end

%% ----------------------------------------- GET GLOBAL ATLAS DATA -----------------------------------------
atlas = cell(size(Arc));
for j = 1:length(Arc)
   atlas{j} = [Arc{j}{:}]; % vectorize charts 
end
totalChart = numel([atlas{:}]);

% find maximum error
thisMaxError = zeros(1,length(Arc));
for j = 1:length(Arc)
    thisMaxError(j) = max([atlas{j}.ErrorBound]);
end
maxError = max(thisMaxError);

%% Plot the local and global manifold
whichFace = true(1,8); % [SE SW WS WN NW NE EN ES]

% restrict view to our compact box
viewBoxCenter = boxCenter;
viewBoxRadius = boxRadius;
ax(1:2:6) = viewBoxCenter - viewBoxRadius;
ax(2:2:6) = viewBoxCenter + viewBoxRadius;

% setup figure
getRes = get(0,'screensize');
screenRes = getRes(3:4);
figAnchor = .1*screenRes;
figSize = .8*screenRes;

close all
myfig = figure('OuterPosition',[figAnchor,figSize]);
hold on
colorbar

% plot manifold patches
tdivs = 400; % density of evaluation nodes in time 
sdivs = [25,1]; % spatial density for evaluating:  [minimum nodes, minimum manifold distance]
% increase tdivs,sdivs(1) to smooth out plots and remove gaps due to sparse evaluation. This significantly increases computation time though.

for j = 1:8
    if whichFace(j)
        lorenzpatch(Arc{j},sdivs,tdivs);
    end
end
plotlocal(P,[0,1,0],[0,0,1]) % plot local manifold
refOrbit % plot a reference orbit on the attractor

% restrict to the box
view(-34,22);
axis tight
zz = axis;
axis([ax(1:4),zz(5:6)]);
title({sprintf('Max C^0 error: %d',max(thisMaxError)),sprintf('Number of charts: %d',totalChart)})




