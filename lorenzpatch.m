function varargout = lorenzpatch(timestepArray,sdivs,tdivs,varargin)
% Adds patch plot of given atlas of charts

% Written by S.K. 01/2017
if nargin ==4
    alphaType = varargin{1};
elseif nargin ==5
    alphaType = varargin{1};
    tdir = varargin{2};
elseif nargin ==6
    alphaType = varargin{1};
    tdir = varargin{2};
    patchColor = varargin{3};
else
    alphaType = 'distance';
    tdir = -1;
end

gcf;
warning('off','MATLAB:delaunay:DupPtsDelaunayWarnId') % suppress duplicate vertices in triangulation warning
% plotting parameters
spaceNodeMin = sdivs(1);
spaceNodeDist = sdivs(2);

if length(tdivs) == 1 %tdivs is the number of time steps
    % get max Time
    maxTime = 0;
    for j = 1:length(timestepArray)
        for k = 1:length(timestepArray{j})
            maxTime = max([maxTime,abs(timestepArray{j}(k).TimeSpan(1) + tdir*timestepArray{j}(k).Tau)]);
        end
    end
    timespan = [0,max([tdir,tdir*maxTime])];
    timeNodes = linspace(timespan(1),timespan(2),tdivs);
else % tdivs is a vector of time steps
    timeNodes = tdivs;
end

if iscell(timestepArray)
    face = [];
    vertex = [];
    colorData = [];
    alphaData = [];
    for j = 1:length(timestepArray)
        for k = 1:length(timestepArray{j})
            [x,y,z] = timestepArray{j}(k).eval([linspace(-1,1,spaceNodeMin)',zeros(spaceNodeMin,1)]);
            dist = sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2);
            numNodes = max([spaceNodeMin,round(sum(dist)./spaceNodeDist)]);
            obj = timestepArray{j}(k);
            spaceNodes = obj.rk45subdiv([obj.TimeSpan(1),obj.TimeSpan(1) + tdir*obj.Tau],numNodes,[4*numNodes,50]);
            [v,f,cData] = timestepArray{j}(k).tri(spaceNodes,timeNodes);
            vertexCount = size(vertex,1);
            face = [face;f + vertexCount];
            vertex = [vertex;v];
            colorData = [colorData;cData];
        end
    end
    if strcmp(alphaType,'distance')
        alphaMin = .01;
        alphaMax = .5;
        distMax = 185;
        alphaByDistance = max(abs(vertex'))'/distMax;
        vAlpha = alphaByDistance;
        alphaData = 64*alphaMax*(1 - (1-alphaMin)*vAlpha);
        % p = patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData,'FaceColor','interp','EdgeColor', 'none','FaceVertexAlphaData',alphaData,'FaceAlpha','interp');
        if exist('patchColor')
            p = patch('Faces', face, 'Vertices', vertex,'FaceColor',patchColor,'EdgeColor','none','FaceVertexAlphaData',alphaData,'FaceAlpha','interp','AlphaDataMapping','direct');                
        else
            p = patch('Faces', face, 'Vertices', vertex,'FaceVertexCData', colorData,'FaceColor','interp','EdgeColor','none','FaceVertexAlphaData',alphaData,'FaceAlpha','interp','AlphaDataMapping','direct');
        end
    elseif strcmp(alphaType,'time')
        alphaMin = .2;
        alphaMax = .6;
        timeEnd = timespan(2);
        alphaByTime = colorData/timeEnd;
        vAlpha = alphaByTime;
        alphaData = 64*(alphaMax*(1-vAlpha) + alphaMin*vAlpha);
        if exist('patchColor')
            p = patch('Faces', face, 'Vertices', vertex,'FaceColor',patchColor,'EdgeColor','none','FaceVertexAlphaData',alphaData,'FaceAlpha','interp','AlphaDataMapping','direct');
        else
            p = patch('Faces', face, 'Vertices', vertex,'FaceVertexCData', colorData,'FaceColor','interp','EdgeColor','none','FaceVertexAlphaData',alphaData,'FaceAlpha','interp','AlphaDataMapping','direct');
        end
    elseif strcmp(alphaType,'none')
        if exist('patchColor')
            p = patch('Faces', face, 'Vertices', vertex,'FaceColor',patchColor,'EdgeColor','none');
        else
            p = patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData, 'FaceColor','interp','EdgeColor', 'none');
        end
    elseif strcmp(alphaType,'fixed')
        fixedAlpha = .6;
        if exist('patchColor')
            p = patch('Faces', face, 'Vertices', vertex,'FaceColor',patchColor, 'EdgeColor', 'none','FaceAlpha',fixedAlpha);
        else
            p = patch('Faces', face, 'Vertices', vertex, 'FaceVertexCData', colorData, 'FaceColor','interp','EdgeColor', 'none','FaceAlpha',fixedAlpha);
        end
    end
    fprintf('Number of Faces: %d \n',size(face,1));
    if nargout == 1
        varargout{1} = p;
    elseif nargout ==3;
        varargout{1} = vertex;
        varargout{2} = face;
        varargout{3} = colorData;
    end
end
warning('on','MATLAB:delaunay:DupPtsDelaunayWarnId') % turn warnings back on

end
