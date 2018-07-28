function oneStep = onesteplorenz(Xinit,Yinit,Zinit,parameter,modes,currentTime,errTol,numSubDivs,parentId,tdir,initError,varargin)
% Returns a single timestep for the advected image of the given initial conditions in the Lorenz example (subdividing if necessary). 

% Written by S.K. 06/2017

% INPUT: 
% Xinit, Yinit, Zinit: Row vectors of length N specifying the Taylor coefficients for the initial curve: gamma(t) = (X(t),Y(t),Z(t)).
% parameter = [rho,sigma,beta]: vector specifying the Lorenz parameters.
% modes = [M,N]: non-negative integers specifying the truncation size.
% currentTime = {t_0,tauMax}: t_0 is the initial time for this timestep, tauMax specifies the largest timeStep to allow.  
% errTol = epsilon: single step precision loss threshold (in decimal digits). Arc will be subdivided if propagation error exceeeds this. 
% numSubDivs = K: Number of subarcs to use in case of subdivision. 
% parentId: handle object which refers to the immediate predecessor of this arc in the atlas. 
% tdir: integration time direction specified as either 1 (unstable manifold) or -1 (stable manifold).
% initError: Maximum ell^1 validation error from last timestep. 

% OUPUT: 
% oneStep: A vector of objects of type = lorenztimestep. 

obj = lorenztimestep(Xinit,Yinit,Zinit,parameter,modes,'timeSpan',currentTime); % attempts to integrate the entire arc for 1 timestep
obj.InitialError = initError;
switch nargin
    case 11
        lastTry = 0;
    case 12
        lastTry = varargin{1};
    otherwise
        disp('Something broke')
end

obj.validate(initError); % validate this timestep
switch tdir
    case -1 % stable manifold goes backwards
        maxErr = errTol;
        if obj.ErrorProp < maxErr || abs(obj.ErrorProp - lastTry) < 1 % Error propagation does not exceed threshold
            oneStep = obj;
            oneStep.Parent = parentId;
        elseif isequal(obj.ErrorBound,Inf)
            sprintf('Validation failed after %.4d time units. Increase order or decrease subdivision threshold to integrate longer',-obj.TimeSpan(2))
            oneStep = obj;
            oneStep.Parent = parentId;
        else
            try % error propagation exceeds threshold so subdivide.
                disp('subdividing');
                divNodes = obj.rk45subdiv([obj.TimeSpan(1),obj.TimeSpan(1) + tdir*3*obj.Tau],numSubDivs);

                xSubs = obj.InitCoef(1).subdivide([divNodes(1:end-1)',divNodes(2:end)'],0);
                ySubs = obj.InitCoef(2).subdivide([divNodes(1:end-1)',divNodes(2:end)'],0);
                zSubs = obj.InitCoef(3).subdivide([divNodes(1:end-1)',divNodes(2:end)'],0);
                
                % shrinkwrap
                newError = 0;
                for j = 1:numSubDivs
                    [xSubArcs(j,:),xError] = shrinkwrap(xSubs(j,:),modes(2),initError);
                    [ySubArcs(j,:),yError] = shrinkwrap(ySubs(j,:),modes(2),initError);
                    [zSubArcs(j,:),zError] = shrinkwrap(zSubs(j,:),modes(2),initError);
                    newError = max([newError,xError,yError,zError]);
                end
                
                oneStep = onesteplorenz(xSubArcs(1,:),ySubArcs(1,:),zSubArcs(1,:),parameter,modes,currentTime,errTol,numSubDivs,parentId,tdir,newError,obj.ErrorProp);
                for j = 2:numSubDivs
                    oneStep = [oneStep,onesteplorenz(xSubArcs(j,:),ySubArcs(j,:),zSubArcs(j,:),parameter,modes,currentTime,errTol,numSubDivs,parentId,tdir,newError,obj.ErrorProp)];
                end
            catch ME
                error('something broke')
            end
        end
    case 1
        errVStime = obj.ErrorProp/obj.Tau
        if errVStime < errTol ||  (errVStime > .9*lastTry && lastTry > 0)
            oneStep = obj;
            oneStep.Parent = parentId;
        else
            try
                disp('subdividing');
                divNodes = obj.rk45subdiv([obj.TimeSpan(1),obj.TimeSpan(1) + tdir*3*obj.Tau],numSubDivs);
                xSubArcs = obj.InitCoef(1).subdivide([divNodes(1:end-1)',divNodes(2:end)'],0);
                ySubArcs = obj.InitCoef(2).subdivide([divNodes(1:end-1)',divNodes(2:end)'],0);
                zSubArcs = obj.InitCoef(3).subdivide([divNodes(1:end-1)',divNodes(2:end)'],0);
                oneStep = onesteplorenz(xSubArcs(1,:),ySubArcs(1,:),zSubArcs(1,:),parameter,modes,currentTime,errTol,numSubDivs,parentId,tdir,initError,errVStime);
                for j = 2:numSubDivs
                    oneStep = [oneStep,onesteplorenz(xSubArcs(j,:),ySubArcs(j,:),zSubArcs(j,:),parameter,modes,currentTime,errTol,numSubDivs,parentId,tdir,initError,errVStime)];
                end
            catch ME
                disp('something broke')
            end
        end
end
end

