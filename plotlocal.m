function plotlocal(P,boundaryColor,interiorColor,varargin)
% add local stable/unstable manifold plots

% Written by S. Kepley 03/2018

% parse input and varargin
p = inputParser;
addRequired(p,'P')
addRequired(p,'boundaryColor')
addRequired(p,'interiorColor')
addParameter(p,'facealpha',1);
addParameter(p,'boundarywidth',.5);

parse(p,P,boundaryColor,interiorColor,varargin{:})
faceAlpha = p.Results.facealpha;
boundaryWidth = p.Results.boundarywidth;

% unpack node args
numBoundaryNode = 50;
numInteriorNode = 50;

% boundary evaluation nodes
N = size(P,1) - 1;
M = size(P,2) - 1;
parmBoundary = linspace(-1,1,numBoundaryNode);
xB = [parmBoundary,ones(1,numBoundaryNode),fliplr(parmBoundary),-1*ones(1,numBoundaryNode)];
yB = [-1*ones(1,numBoundaryNode),parmBoundary,ones(1,numBoundaryNode),fliplr(parmBoundary)];
XB = bsxfun(@power,xB,(0:N)');
YB = bsxfun(@power,yB',(0:M));


% interior evaluation nodes
parmInterior = linspace(-1,1,numInteriorNode);
[xI,yI] = meshgrid(parmInterior,parmInterior);
data = [reshape(xI,[],1),reshape(yI,[],1)];
XI = data(:,1);
YI = data(:,2);
face = delaunay(XI,YI); % index triples for delaunay triangulation in space-time.
interiorNode = data;

%% plot local manifold interior
interior = lorenztimestep();
interior.Coord = BAscalar(squeeze(P(:,:,1)));
interior.Coord(2) = BAscalar(squeeze(P(:,:,2)));
interior.Coord(3) = BAscalar(squeeze(P(:,:,3)));
[x,y,z] = interior.eval(interiorNode); % cooresponding coordinates for each vertex
vertex = [x,y,z];
patch('Faces', face, 'Vertices', vertex,'FaceColor',interiorColor,'EdgeColor','none','FaceAlpha',faceAlpha);

% plot stable manifold boundary
boundary = zeros(3,4*numBoundaryNode);
for j = 1:3
    AA = squeeze(P(:,:,j))*XB;
    BB = YB'.*AA;
    boundary(j,:) = real(sum(BB,1));
end
X = boundary(1,:);
Y = boundary(2,:);
Z = boundary(3,:);
plot3(X,Y,Z,'Color',boundaryColor,'LineWidth',boundaryWidth); % initial conditions
end

