classdef lorenztimestep < handle
% An analytic arc segment and its evolution in time under the Lorenz flow parameterized as a 2-variable Taylor series. 

% Written by S.K. 05/2016

% FIXED: 
% Fast operator inversion added 02/2017
% Deprecated the lorenz_arc class 12/2016

% TO DO: 
% Delegate functions inside class constructor. constructor should do no work
% Implement as a subclass of TaylorTimestep superclass
% Add support for vectorized method calls
% Implement surface area computation

    properties
        Coord; % vector of BAscalars for each coordinate of this chart 
        Modes; % truncation size in [space,time] directions
        Tau=1; % vector field (time) rescaling
        TimeSpan; % interval in time for this chart
        MTCrange = [-1,1]; % material (spatial) coordinates this chart
        ErrorBound; % analytic error bound 
    end
	
	properties(Hidden = true)
        Parameter = [28,10,8/3]; %[rho,sigma,beta]
		Weight = 'ones';
		MaxTau; % maximal tau for this timestep 
        InitCoef; % (BAscalar) form of initial data
        RadiiPoly; % Radii polynomial data for validation of this step 
        ErrorProp; % decimal digits of precision loss for this step
        SubArcDepth; % keep track of sub-division depth and break out of infinite sub-division loops
        InitialError = 0; % error at the start of this timestep
		Parent; % handle for the preceding timestep
	end
    
    methods
        %% =================================== CLASS CONSTRUCTOR METHODS ===================================
        function obj = lorenztimestep(Xinit,Yinit,Zinit,parameter,modes,varargin)
		% class constructor
            if(nargin > 0)
                p = inputParser;
                addRequired(p,'Xinit')
                addRequired(p,'Yinit')
                addRequired(p,'Zinit')
                addRequired(p,'parameter')
                addRequired(p,'modes')
                addParameter(p,'timeSpan',0)
                
                %parse varargs
                parse(p,Xinit,Yinit,Zinit,parameter,modes,varargin{:})
                obj.TimeSpan = p.Results.timeSpan;
                
                % set properties
                obj.Parameter = parameter;
                obj.Modes = modes;
                obj.InitCoef = [BAscalar(Xinit,obj.Modes);BAscalar(Yinit,obj.Modes);BAscalar(Zinit,obj.Modes)];
                switch length(modes)
                    case 1 % initial date is a point
                        obj.Coord = [BAscalar(Xinit);BAscalar(Yinit);BAscalar(Zinit)];
                        generatecoef(obj); % No scaling on IVP solver.
                    case 2 % initial data is an arc
                        obj.Coord = [BAscalar(Xinit,[1,modes(2)]);BAscalar(Yinit,[1,modes(2)]);BAscalar(Zinit,[1,modes(2)])];
                        generatecoef(obj);
                        if length(obj.TimeSpan) == 1 || isa(obj.TimeSpan,'cell') 
                            row_decay = obj.decay;
                            L = (1e-16/row_decay)^(1/(obj.Modes(1)-1));
                            if length(obj.TimeSpan) == 1 % rescale coefficients automatically based on coefficient (in time) decay
                                obj.scaletime(L);
                            else
                                obj.scaletime(min([L,obj.TimeSpan{2}])); % take a timestep based on coefficient decay but not to exceed Timespan{2} units
                            end
                        elseif isa(obj.TimeSpan,'double') % take a fixed timestep
                            obj.scaletime(diff(obj.TimeSpan));
                        end
                end
            end
        end
        
        function generatecoef(obj)
		% generate the coefficients for this timestep by recursively solving the differential equation
            switch length(obj.Modes)
                case 1
                    for m = 1:obj.Modes - 1
                        XZ = recursive_conv(obj.Coord(1).Coef,obj.Coord(3).Coef,'1d');
                        XY = recursive_conv(obj.Coord(1).Coef,obj.Coord(2).Coef,'1d');
                        obj.Coord(1).Coef(m+1) = (obj.Parameter(2)/m)*(obj.Coord(2).Coef(m) - obj.Coord(1).Coef(m));
                        obj.Coord(2).Coef(m+1) = (1/m)*(obj.Parameter(1)*obj.Coord(1).Coef(m) - XZ - obj.Coord(2).Coef(m));
                        obj.Coord(3).Coef(m+1) = (1/m)*(XY - obj.Parameter(3)*obj.Coord(3).Coef(m));
                    end
                case 2
                    for m = 1:(obj.Modes(1)-1)
                        % cauchy products for non-linear terms
                        XZ = conv2([zeros(size(obj.Coord(1).Coef)-[0,1]),obj.Coord(1).Coef],obj.Coord(3).Coef,'valid');
                        XY = conv2([zeros(size(obj.Coord(1).Coef)-[0,1]),obj.Coord(1).Coef],obj.Coord(2).Coef,'valid');
                        obj.Coord(1).Coef(m+1,:) = (obj.Parameter(2)/m)*(obj.Coord(2).Coef(m,:) - obj.Coord(1).Coef(m,:));
                        obj.Coord(2).Coef(m+1,:) = (1/m)*(obj.Parameter(1)*obj.Coord(1).Coef(m,:) - XZ - obj.Coord(2).Coef(m,:));
                        obj.Coord(3).Coef(m+1,:) = (1/m)*(XY - obj.Parameter(3)*obj.Coord(3).Coef(m,:));
                    end
            end
            for jj = 1:3
                obj.Coord(jj).Modes = obj.Modes;
            end
        end
        
        function ddt = diff(obj)
			% return time derivative for this chart
            ddt = [obj.Coord(1).dt;obj.Coord(2).dt;obj.Coord(3).dt];
        end

        function coefs = taumap(obj,Tau)
            % evaluation of time-tau map
            coefs = [obj.Coord(1).fixtime(Tau);obj.Coord(2).fixtime(Tau);obj.Coord(3).fixtime(Tau)];
        end
                   
        function scaletime(obj,L)
		% rescale time for this chart by L
		
            obj.Coord(1).scaletime(L/obj.Tau); % new coefficients correspond to image under the time-rescaled flow. 
            obj.Coord(2).scaletime(L/obj.Tau);
            obj.Coord(3).scaletime(L/obj.Tau);
            obj.Tau = L;
            if isa(obj.TimeSpan,'cell')
                obj.TimeSpan = [obj.TimeSpan{1},L];
            else
                obj.TimeSpan(2) = L;
            end
        end
                
        function objNorm = norm(obj)
		% returns the max norm for this chart
            objNorm = max(norm(obj.Coord));
        end
        
        function decay = decay(obj)
		% returns the coefficient decay for this chart 
            decay = max([obj.Coord(1).decay,obj.Coord(2).decay,obj.Coord(3).decay]);
        end
        
        
        %% =================================== EVALUATION METHODS ===================================
        
        function varargout = eval(obj,data)
		% evaluation method for obtaining image of this chart in phase space. 
		
		
            % INPUT: 
			% data: m-by-(k+1) array of the form [S1,S2,...,Sk,T];
            x = real(obj.Coord(1).eval(data));
            y = real(obj.Coord(2).eval(data));
            z = real(obj.Coord(3).eval(data));
            switch nargout
                case 1
                    varargout{1} = [x;y;z];
                case 3
                    varargout{1} = x;
                    varargout{2} = y;
                    varargout{3} = z;
            end
        end
        
        function [vertex,face,cData] = tri(obj,spaceNodes,timeNodes)
		% evaluate this chart on a Delauney triangulation of space-time grid 
		
            % INPUT:
			% sNodes, tNodes are vectors or REAL time and RELATIVE space coordinates 
            
            % pick out timeNodes which lie in this timestep
            tdir = sign(timeNodes(end)-timeNodes(1)); % direction of the flow
            initTime = obj.TimeSpan(1); % initial time for this timestep
            finalTime = initTime + tdir*obj.Tau; % final time for this timestep
            timeSlice = tdir*initTime < tdir*timeNodes & tdir*timeNodes < tdir*finalTime; % index for tNodes which lie in this timestep
            
            if ~any(timeSlice)
                face = [];
                vertex = [];
                cData = [];
                return
            end
           
            % real time coordinates
            if tdir < 0
                realTime = [initTime,timeNodes(timeSlice),max([finalTime,timeNodes(end)])]; % 
            else
                realTime = [initTime,timeNodes(timeSlice),min([finalTime,timeNodes(end)])]; 
            end
            materialTime = obj.mtc(realTime); % relative time coordinates
            if numel(materialTime) ==1
                numVertex = length(spaceNodes);
                data = [reshape(spaceNodes,[],1),materialTime*ones(numVertex,1)];
                [X,Y,Z] = obj.eval(data);
                vertex = [X,Y,Z];
                face = [(1:numVertex-2)',(2:numVertex-1)',(3:numVertex)'];
                cData = realTime*ones(numVertex,1);
            else
                % triangulate flow values on valid space-time coordinates
                [ss,tt] = meshgrid(spaceNodes,materialTime);
                data = [reshape(ss,[],1),reshape(tt,[],1)];
                s = data(:,1);
                t = data(:,2);
                try
                    face = delaunay(s,t); % index triples for delaunay triangulation in space-time.
                    [X,Y,Z] = obj.eval(data); % coordinates for each vertex
                    vertex = [X,Y,Z];
                    cData = repmat(realTime',length(spaceNodes),1);
                catch ME
                    face = [];
                    vertex = [];
                    cData = [];
                    disp('Empty triangles returned. Something is wrong')
                end
            end
        end
        
        
        function MTC = mtc(obj,T)
		% scale a vector of real time coordinates, T, to material time coordinates in [-1,1] relative to this time step.
            MTC = (T - obj.TimeSpan(1))/obj.Tau;
        end
		
		function nodes = rk45subdiv(obj,timeSpan,numArcs,varargin)
		% returns spatial sub-division nodes by using RK45 to estimate the (spatial) Lyapunov exponents locally (in time). 
            switch nargin
                case 3
                    gridDensity = [2500,200]; % use a dense space-time grid for numerics 
                otherwise
                    gridDensity = varargin{1};
            end
			% approximate subdivision using RK45 heuristics
			s = linspace(-1,1,gridDensity(1));
			T = linspace(timeSpan(1),timeSpan(2),gridDensity(2));
			evl = @(arc,s)polyval(fliplr(arc),s); % evaluation of polynomial arc
			
			
			%RK45 test integration
			lorenz = @(t,y)vectorizeinitdata(t,y);
			ICs = nan(3*length(s),1);
			ICs(1:3:end) = evl(obj.InitCoef(1).Coef(1,:),s)';
			ICs(2:3:end) = evl(obj.InitCoef(2).Coef(1,:),s)';
			ICs(3:3:end) = evl(obj.InitCoef(3).Coef(1,:),s)';
			[~,sol] = ode45(lorenz,T,ICs);
			PX = sol(:,1:3:end)';
			PY = sol(:,2:3:end)';
			PZ = sol(:,3:3:end)';
			
			% approximate lyapunov exponents
			dX = diff(diff(PX),1,2);
			dY = diff(diff(PY),1,2);
			dZ = diff(diff(PZ),1,2);
			d_tot = abs(dX) + abs(dY) + abs(dZ);
			dS_tot = sum(d_tot,2);
			intS = [0,cumsum(dS_tot)'];
			davg = intS(end)/numArcs;
			nodes = [-1,nan(1,numArcs)];
			for j = 1:numArcs
				next_ind = find(intS > davg,1);
				nodes(j+1) = mean([s(next_ind-1),s(next_ind)]);
				intS = intS - davg;
			end
			nodes(end) = 1;
		end	
          
        %% =================================== VALIDATION METHODS ===================================
        
        function Feval = zeroMap(obj,varargin)
		% compute F for polynomial approximation
		
			if nargin > 1
				numType = varargin{1}; 
			else
				numType = 'intval';
			end
            
			switch numType
                case 'double' % non-rigorous. Use only for testing/debugging 
                    % add/subtract gamma has no effect on the finite part so it is omitted
                    F1 = obj.Coord(1).dt - obj.Tau*obj.Parameter(2)*shift(obj.Coord(2)-obj.Coord(1));
                    F2 = obj.Coord(2).dt - obj.Tau*shift(obj.Parameter(1)*obj.Coord(1) - mtimes(obj.Coord(1),obj.Coord(3),'Full') - obj.Coord(2));
                    F3 = obj.Coord(3).dt - obj.Tau*shift(mtimes(obj.Coord(1),obj.Coord(2),'Full') - obj.Parameter(3)*obj.Coord(3));
                    Feval = [F1;F2;F3];
                case 'intval' % rigorous enclosure for big F map of obj
					intvalCoord = obj.Coord.intval;
                    zeroOrder = [BAscalar(intvalCoord(1).Coef(1,:),obj.Modes);BAscalar(intvalCoord(2).Coef(1,:),obj.Modes);BAscalar(intvalCoord(3).Coef(1,:),obj.Modes)];
                    F1 = intvalCoord(1).dt + obj.InitCoef(1).intval - zeroOrder(1) - obj.Tau*obj.Parameter(2)*shift(intvalCoord(2) - intvalCoord(1));
                    F2 = intvalCoord(2).dt + obj.InitCoef(2).intval - zeroOrder(2) - obj.Tau*shift(obj.Parameter(1)*intvalCoord(1) - mtimes(intvalCoord(1),intvalCoord(3),'Full') - intvalCoord(2));
                    F3 = intvalCoord(3).dt + obj.InitCoef(3).intval - zeroOrder(3) - obj.Tau*shift(mtimes(intvalCoord(1),intvalCoord(2),'Full') - obj.Parameter(3)*intvalCoord(3));
                    Feval = [F1;F2;F3];
            end
        end
        
        function rPoly = validate(obj,varargin)
		% Computes validated error bounds for analytic error (ell^1) of the polynomial chart 
            
            % ----------------------- F,DF, A, and A_dagger -----------------------
            if(nargin>1)
                initialError = varargin{1};
            else
                initialError = 0;
            end
            MN = prod(obj.Modes);

			% F = obj.zeroMap('double');
            F = obj.zeroMap();
            
            % prime operator given by its action on vectors by Ah = h'
            Idprime = BAoperator(diag(repmat([1,1:obj.Modes(1)-1],1,obj.Modes(2))),obj.Modes);
            
            % shift operator given by its action on vectors by Ah = eta(h)
            shift_subdiag = ones(1,MN);
            shift_subdiag(obj.Modes(2):obj.Modes(2):end) = zeros(1,obj.Modes(1));
            shift_matrix = diag(shift_subdiag);
            Idshift = BAoperator([shift_matrix(end,:);shift_matrix(1:end-1,:)],obj.Modes);
            
            % finite part of DF as 9 BAoperators (components of the Jacobian)
            DF(1,1) = Idprime + obj.Tau*obj.Parameter(2)*Idshift; %DxF1
            DF(1,2) = -obj.Tau*obj.Parameter(2)*Idshift; %DyF1
            DF(1,3) = BAoperator(zeros(MN),obj.Modes); %DzF1
            DF(2,1) = -obj.Tau*obj.Parameter(1)*Idshift + obj.Tau*leftTimesOperator(obj.Coord(3).shift); %DxF2
            DF(2,2) = Idprime + obj.Tau*Idshift; %DyF2
            DF(2,3) = obj.Tau*leftTimesOperator(obj.Coord(1).shift); %DzF2
            DF(3,1) = -obj.Tau*leftTimesOperator(obj.Coord(2).shift); %DxF3
            DF(3,2) = -obj.Tau*leftTimesOperator(obj.Coord(1).shift); %DyF3
            DF(3,3) = Idprime + obj.Tau*obj.Parameter(3)*Idshift; %DzF3
            
            DF_matrix = DF.block;
            A_matrix = inv(DF); % inverse operator 
            
            % A to BAoperator
            for j = 1:3
                for k = 1:3
                    A(j,k) = BAoperator(A_matrix(1+(j-1)*MN:j*MN,1+(k-1)*MN:k*MN),obj.Modes);
                end
            end
            
            %% ----------------------- Radii Polynomial Bounds -----------------------
            r_star = max([100*initialError,1e-13]); % initial guess for Lipschitz bound
            Anorms = [sum(arrayfun(@(j)A(1,j).norm,1:3));sum(arrayfun(@(j)A(2,j).norm,1:3));sum(arrayfun(@(j)A(3,j).norm,1:3))];
            normA = max(Anorms); % operator norm ||A||

            % Y0
            F_MN = [BAscalar(F(1).Coef,obj.Modes),BAscalar(F(2).Coef,obj.Modes),BAscalar(F(3).Coef,obj.Modes)];
            y1 = A(1,1)*F_MN(1) + A(1,2)*F_MN(2) + A(1,3)*F_MN(3); %[A*F]_1
            Y_1 = y1.norm;
            
%             y2 = mtimes(obj.Coord(1),obj.Coord(3),'Full'); %full X*Z including spillover terms
            y2 = BAscalar(F(2).Coef);
            AF_2 = (A(2,1)*F_MN(1) + A(2,2)*F_MN(2) + A(2,3)*F_MN(3)); %convolution without spillover terms
            y2.Coef(1:obj.Modes(1),1:obj.Modes(2)) = AF_2.Coef;
            Y_2 = y2.norm;
            
%             y3 = mtimes(obj.Coord(1),obj.Coord(2),'Full'); %full X*Y including spillover terms
            y3 = BAscalar(F(3).Coef);
            AF_3 = (A(3,1)*F_MN(1) + A(3,2)*F_MN(2) + A(3,3)*F_MN(3)); %convolution without spillover terms
            y3.Coef(1:obj.Modes(1),1:obj.Modes(2)) = AF_3.Coef;
            Y_3 = y3.norm;
            
            Y0 = max(sup([Y_1,Y_2,Y_3])) + initialError;
            
            % Z0
            ADF_matrix = A_matrix*DF_matrix;
            Id_MN = eye(size(ADF_matrix));
            Id_ADF = Id_MN - ADF_matrix;
            for j = 1:3
                for k = 1:3
                    I_ADF(j,k) = BAoperator(Id_ADF(1+(j-1)*MN:j*MN,1+(k-1)*MN:k*MN),obj.Modes);
                end
            end
            IADFnorms = [sum(arrayfun(@(j)I_ADF(1,j).norm,1:3));sum(arrayfun(@(j)I_ADF(2,j).norm,1:3));sum(arrayfun(@(j)I_ADF(3,j).norm,1:3))];
            Z0 = max(IADFnorms);

            % Z1
            Z1 = (obj.Tau./obj.Modes(1))*max([2*obj.Parameter(1),obj.Parameter(2) + obj.Coord(3).norm + 1 + obj.Coord(1).norm, obj.Coord(2).norm + obj.Coord(1).norm + obj.Parameter(3)]);
            
            % Z2
            Z2 = (obj.Tau*r_star)*max(1/obj.Modes(1),normA);
            
            % radii polynomial
            rPoly.bounds = [Z2,Z1,Z0,Y0];
            rPoly.coefs = [Z2,Z0+Z1 - 1,Y0];
            try
                rPoly.radius = solveradiipoly(rPoly.coefs); % interval quadratic equation solver
            catch ME
                error('Radii polynomial solver encountered a problem')
            end
%             rPoly.radius = min(roots(rPoly.coefs));
            obj.RadiiPoly = rPoly.bounds;
            
            if imag(rPoly.radius) ~= 0 
                obj.ErrorBound = Inf;
				warning('Radii polynomial has no real roots')
            elseif rPoly.radius < 0
                obj.ErrorBound = Inf;
                warning('Radii polynomial has a NEGATIVE root. Something is wrong')
            else
                obj.ErrorBound = rPoly.radius;
%                 fprintf('Error: %.4d \n', [rPoly.radius])
            end
			
            if initialError > 0
                obj.ErrorProp = log10(obj.ErrorBound/initialError); % decimal places of precision lost on this step
            else
                obj.ErrorProp = log10(obj.ErrorBound/eps(1));
            end
            if obj.ErrorProp < 0
                error('Error can never decrease in a timestep. Something is wrong')
            end
        end
		
		%% ===================================DEPRECATED METHODS ===================================
		% These should never be called anymore and will return an error if they are. Remove these if no calls occur for some time. 
		
		
		function plot(obj,s,t,varargin) % deprecated 03/2018
			error('This method is deprecated') 
			
            % if nargin > 3
                % plot_color = varargin{1};
            % else
                % plot_color = 'g';
            % end
            % [x,y,z] = obj.mesh_eval(s,t);
            % x = x'; y = y'; z = z';
            % plot3(x,y,z,'LineWidth',1.2,'Color',plot_color);
        end
		
		
		function newStep = advect(obj,timeDirection,errTol,numSubDivs,maxSubDiv) % deprecated 03/2018
		% this should never be called 
			error('This method is deprecated')

			% % errTol: maximum error per time unit
			% % timeDirection: 1 or -1 for direction to flow in time 
			% % numSubDivs: number of subdivisions to take if error tolerance is exceeded
			% % maxSubDiv: maximum number of subdivision attempts before integration stops
			% obj.validate(obj.InitialError);
			% maxErr = obj.Tau*errTol;
			% if obj.ErrorProp < maxErr
				% newStep = obj;
			% elseif obj.SubArcDepth == maxSubDiv
				% set(obj,'Tau',0)
				% newStep = obj;
			% else
				% disp('subdividing');
				% divNodes = obj.rk45subdiv([obj.TimeSpan(1),obj.TimeSpan(1) + timeDirection*5*obj.Tau],numSubDivs);
				% xSubArcs = obj.InitCoef(1).subdivide([divNodes(1:end-1)',divNodes(2:end)'],0);
				% ySubArcs = obj.InitCoef(2).subdivide([divNodes(1:end-1)',divNodes(2:end)'],0);
				% zSubArcs = obj.InitCoef(3).subdivide([divNodes(1:end-1)',divNodes(2:end)'],0);
				% if obj.MaxTau < Inf;
					% timeSpan = {obj.TimeSpan(1),obj.MaxTau};
				% else
					% timeSpan = obj.TimeSpan(1);
				% end
				% for j = 1:numSubDivs
					% newStep(j) = lorenztimestep(xSubArcs(j,:),ySubArcs(j,:),zSubArcs(j,:),obj.Parameter,obj.Modes,'timeSpan',timeSpan);
					% newStep(j).SubArcDepth = obj.SubArcDepth + 1;
					% newStep(j).MTCrange = mean(obj.MTCrange) + divNodes(j,j+1)*diff(obj.MTCrange);
				% end
			% end
		end % end advect
		
		
		function patch(obj,s,t,varargin) % deprecated 03/2018
			error('This method is deprecated. Use the lorenzpatch function')
			
            % % s,t are vectors of evaluation points. 
            % switch length(obj)
                % case 1
                    % if nargin > 3
                        % patch_color = varargin{1};
                    % else
                        % patch_color = 'g';
                    % end
                    % S = length(s); T = length(t);
                    % [Xs,Ys,Zs] = obj.mesh_eval(s,t(1));
                    % [Xs(S+1:S+T),Ys(S+1:S+T),Zs(S+1:S+T)] = obj.mesh_eval(s(end),t);
                    % [Xs(S+T+1:2*S+T),Ys(S+T+1:2*S+T),Zs(S+T+1:2*S+T)] = obj.mesh_eval(fliplr(s),t(end));
                    % [Xs(2*S+T+1:2*(S+T)),Ys(2*S+T+1:2*(S+T)),Zs(2*S+T+1:2*(S+T))] = obj.mesh_eval(s(1),fliplr(t));
                    % gcf;
                    % patch(Xs,Ys,Zs,patch_color,'EdgeColor',patch_color);
                % otherwise
                    % for j = 1:length(obj)
                        % obj(j).patch(s,t,varargin{:})
                    % end
            % end
        end  
		
		            
        function [x,y,z] = mesh_eval(obj,s,t) % deprecated 03/2018
			error('This method is deprecated')
			
            % x = real(obj.Coord(1).gridEval(s,t));
            % y = real(obj.Coord(2).gridEval(s,t));
            % z = real(obj.Coord(3).gridEval(s,t));
        end %end mesh_eval

    end %end methods
	
end %end class


