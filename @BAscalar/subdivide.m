function varargout = subdivide(obj,parmRange,varargin)
% rigorously subdivide a surface given by a Taylor series into subsurfaces

% Written by S. Kepley 04/2017

% ---------------------- INPUT ----------------------
% obj (BAscalar): vector of BAscalars parameterizered on [-1,1]^d (space only) or [-1,1]^{d+1} (space/time)
% parmRange (double): m-by-(2d) array = [row1;row2;...;rowm] where each row is of the form [s11,s12,s21,s22,...,sm1,sm2] corresponds to the closed rectangle [s11,s12]x[s21,s22]x...x[sm1,sm2] in R^d
% varargin{1} = t0 (double): Use if obj is space/time parameterization to subdivide obj(s,t0) which is parameterizered on [-1,1]^d

% ---------------------- OUTPUT ----------------------
% subDomain (double): m-by-(2d) array of material coordinates with respect to [-1,1]^d (allows tracking of multiple subdivisions in the larger spatial domain) 
% subSurface (BAscalar): m-by-length(obj) array where subSurface(i,j) = obj(j) defined on Ri and rescaled to [-1,1]^d

% parse input and varargin
p = inputParser;
addRequired(p,'obj')
addRequired(p,'parmRange')
addOptional(p,'t0',[])
addParameter(p,'domain',[]);

parse(p,obj,parmRange,varargin{:})
t0 = p.Results.t0;
objDomain = p.Results.domain;

%% Evaluate obj(s,t0) if specified
if isempty(t0)
    spaceParm = obj;
else
    t0 = varargin{1};
    spaceParm = obj.fixtime(t0);
end

spaceDimension = spaceParm.SurfaceDimension; 
if size(parmRange,2) ~= 2*spaceDimension
    disp('we should stop here')
    error('parmRange should have 2*dimension-many columns')
end

%% Get spatial domain for parent surface
if nargout == 2
    if isempty(objDomain) % domain is [-1,1]^d
        spaceDomain = repmat([-1,1],1,spaceDimension);
    else
        spaceDomain = objDomain;
    end
    
    switch spaceDimension
        case 1
            % obtain new material coordinates
            rescale = @(s)(mean(spaceDomain) + .5*s*diff(spaceDomain));
            newMTC = rescale(parmRange);
        otherwise
            error('not yet implemented for higher dimensions')
    end
end

%% Reparameterize surface into subsurfaces
switch spaceDimension
    case 1 % spaceParm is a coefficient vector for a 1-d spatial parameterization. 
        if length(spaceParm) ==1 % spaceParm is a single BAscalar
            numSubSurface = size(parmRange,1); % number of subSurfaces
            if numSubSurface == 1 % parmRange specifies a single subSurface to be parameterized centered at (s1+s2)/2 with radius (s2-s1)/2.
                newCenter = intval(mean(parmRange)); % midpoint
                reScaleBy = intval(newCenter - parmRange(1)); % radius
                N = length(spaceParm.Coef);
                V = pascal(N);
                shiftOperator = intval(zeros(N));
                centerPowers = newCenter.^(0:N-1);
                for j = 1:N
                    shiftOperator(j,j:end) = V(j,1:N-j+1).*centerPowers(1:N-j+1)*reScaleBy^(j-1);
                end
                newParm = (shiftOperator*spaceParm.Coef')';
                
                
            else % each row of parmRange specifies a subsurface
                % if isequal(obj.CoefType,'intval')
                newParm = midrad(zeros(numSubSurface,spaceParm.Modes),0);
                for j = 1:numSubSurface
                    newParm(j,:) = spaceParm.subdivide(parmRange(j,:));
                end
            end
            
        else % spaceParm is a vector of BAscalars
            objLength = length(spaceParm);
            newParm{objLength} = spaceParm(objLength).subdivide(parmRange);
            for j = 1:objLength-1
                newParm{j} = spaceParm(j).subdivide(parmRange);
            end
        end
        
    otherwise % spaceParm is a higher dimension surface
        error('not yet implemented for higher dimensions')
end

switch nargout
    case 1
        varargout{1} = newParm; % no coordinate transform
    case 2
        varargout{1} = newMTC; % new material coordinates
        varargout{2} = newParm;
end
end


