function arcObj = lorenzbox(initData,initError,parameter,modes,timespan,errTol,numSubDivs,boxCenter,boxRadius)
% integrates an initial local parameterization boundary segment until a given time or until it leaves a specified
% box. 

% Written by S. Kepley 06/2017
% Added soft crash for arcs which fail their validation 03/2018

% INPUT: 
% InitError is the local Paramterization error bound. 
% numSubDivs is % the number of subdivisions to perform when error tolerance is exceeded. 

t0 = timespan(1); % initial time
tf = timespan(2); % final time
tdir = sign(tf-t0); % direction to integrate

% evaluation nodes to determine arc segments which have exited the box.
numSpaceNodes = 100; % 
spaceNodes = linspace(-1,1,numSpaceNodes);

if isa(initData,'double') % initial integration 
    Xinit = initData(1,:);
    Yinit = initData(2,:);
    Zinit = initData(3,:);
    currentTime = {t0,.1};
    numT = 1;
    initObj = lorenztimestep();
    initObj.ErrorBound = initError;
    arcObj = {onesteplorenz(Xinit,Yinit,Zinit,parameter,modes,currentTime,errTol,numSubDivs,initObj,tdir,initError)};
else % continue a previous integration by extending the integration time
   arcObj = initData;
   numT = length(initData);
end
%%
stopflag = false;
while ~isempty(arcObj{numT}) && ~stopflag
    clear newStep
    newStep = lorenztimestep();
    for j = 1:length(arcObj{numT})
        stopflag = true;
        thisArc = arcObj{numT}(j);
        
        % check if integration time is exceeded
        Tcurrent = thisArc.TimeSpan(1) + tdir*thisArc.Tau;
        if tdir*Tcurrent <  tdir*tf && ~isequal(thisArc.ErrorBound,Inf) % integration time not yet exceeded and last validation didn't fail
            
            if thisArc.InitialError > 1e-2
                disp('here')
            end
            stopflag = false; % not stopping this round
            currentTime = {Tcurrent,.1};
            
            %  CUT OUT PIECES WHICH HAVE LEFT THE BOX
            [x,y,z] = thisArc.eval([spaceNodes',tdir*ones(numSpaceNodes,1)]);
            img = [x,y,z];
            if length(boxRadius) ==3
                sidx = find(abs(img(:,1) - boxCenter(1)) < boxRadius(1) & abs(img(:,2) - boxCenter(2)) < boxRadius(2)); % & abs(img(:,3) - boxCenter(3)) < boxRadius(3));
            else
                sidx = 1:numSpaceNodes;
            end
            if any(sidx)
                parentId = thisArc;
                parmRange = spaceNodes([max([sidx(1)-1,1]),min([sidx(end)+1,numSpaceNodes])]);
                thisBoundary = thisArc.taumap(tdir);
                
                if isequal(parmRange,[-1,1]) % entire arc remains inside box
                    Xinit = thisBoundary(1).Coef;
                    Yinit = thisBoundary(2).Coef;
                    Zinit = thisBoundary(3).Coef;
                    newError = thisArc.ErrorBound;
                else % cut out the portion which left
                    Xsub = thisBoundary(1).subdivide(parmRange);
                    Ysub = thisBoundary(2).subdivide(parmRange);
                    Zsub = thisBoundary(3).subdivide(parmRange);
                    
                    % shinkWrap to floats and get new initial Error
                    [Xinit,XsubError] = shrinkwrap(Xsub,length(Xsub),thisArc.ErrorBound);
                    [Yinit,YsubError] = shrinkwrap(Ysub,length(Ysub),thisArc.ErrorBound);
                    [Zinit,ZsubError] = shrinkwrap(Zsub,length(Zsub),thisArc.ErrorBound);
                    newError = max([XsubError,YsubError,ZsubError]);
                end
                
                % Then try to integrate the remaining portion
                newStep = [newStep,onesteplorenz(Xinit,Yinit,Zinit,parameter,modes,currentTime,errTol,numSubDivs,parentId,tdir,newError)];
            else
                disp('This arc has exited')
            end
        else
            newStep = [newStep,thisArc]; % done integrating this arc (integration time exceeded or error = Inf). Carry these arcs down
        end
    end
    arcObj{numT+1} = newStep(2:end);
    getTime = [newStep.TimeSpan];
    endTime = getTime(1:2:end) - getTime(2:2:end); % ending integration times for each segment.
    fprintf('Current Time: %.4g \n',(max(endTime)))
    numT = length(arcObj);
end
arcObj = arcObj(1:end-1);
end

