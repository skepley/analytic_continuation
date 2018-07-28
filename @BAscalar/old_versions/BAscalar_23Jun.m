classdef BAscalar < handle
    %% ----------------------------------- TO DO -----------------------------------
    
    % Change from handle class to value class.
    % Vectorize all methods properly.
    % Add sin,cos,exp with error bounds using automatic differentation.
    % verify sqrt works correctly and add error bounds.
    % Test methods with intval coefficients.
    % Change evaluation methods to ndgrid format.
    % Does fixtime/fixSpace need to be evaluated on intervals?
    % implement radii polynomials using polynom.rootbound
    
    % ----------------------------------- DONE (NEEDS ADDITIONAL TESTING) -----------------------------------
    % Add support for intval coefficients.
    % Add multiplication for intval BAscalars.
    
    %% ----------------------------------- Properties -----------------------------------
    properties
        Coef; % Taylor coefficients
        CoefType; % double or interval
        SurfaceDimension;
        Modes; %[M,N] where M is a scalar and N is a length-k vector
    end
    
    properties(Hidden = 1)
        Weight = 'ones'; % default is [1,1,.....1]
        % 		Norm;  commented 12/16. to be removed
        %         Decay; commented 12/16. to be removed
    end
    
    %% ----------------------------------- MATLAB FUNCTIONALITY AND CLASS METHODS -----------------------------------
    methods
        function obj = BAscalar(Coef,varargin)
            %class constructor
            if(nargin > 0)
                if isa(Coef,'cell') % vector of BAscalars
                    objLength = length(Coef);
                    obj(objLength) = BAscalar(Coef{objLength},varargin{:});
                    for j = 1:objLength-1
                        obj(j) = BAscalar(Coef{j},varargin{:});
                    end
                    obj = reshape(obj,size(Coef));
                    
                else % single BAscalar
                    switch nargin
                        case 1 % input is Coef of correct size
                            if isa(Coef,'BAscalar')
                                obj = Coef;
                            else
                                obj.Coef = Coef;
                                obj.CoefType = class(obj.Coef);
                                obj.Modes = size(Coef);
                                if size(Coef,1) == 1;
                                    obj.SurfaceDimension = length(obj.Modes) - 1; % number of spatial dimensions
                                else
                                    obj.SurfaceDimension = length(obj.Modes) - 1;
                                end
                            end
                            
                        otherwise % input is Coef and desired truncation size
                            obj.Modes = varargin{1};
                            obj.SurfaceDimension = length(obj.Modes) - 1;
                            if nargin > 2;
                                obj.Weight = varargin{2};
                            end
                            if isa(Coef,'double') || isa(Coef,'intval') % coefficients given as double or intval array
                                switch obj.SurfaceDimension
                                    case 0
                                        obj.Coef = Coef(1:min(end,obj.Modes));
                                    case 1
                                        obj.Coef = Coef(1:min(end,obj.Modes(1)),1:min(end,obj.Modes(2)));
                                    case 2
                                        if obj.Modes(1) == 1
                                            obj.Coef = Coef(1:min(end,obj.Modes(2)),1:min(end,obj.Modes(3)));
                                        else
                                            obj.Coef = Coef(1:min(end,obj.Modes(1)),1:min(end,obj.Modes(2)),1:min(end,obj.Modes(3)));
                                        end
                                    otherwise
                                        error('Not yet implemented')
                                end
                                obj.CoefType = class(obj.Coef);
                            elseif isa(Coef,'polynom') % coefficients given as Intlab polynom object.
                                if obj.SurfaceDimension ~= 1
                                    error('Not yet implemented')
                                else
                                    obj.Coef = reshape(flip(flip(Coef.c,1),2),obj.Modes(1),[]);
                                end
                                obj.CoefType = class(obj.Coef);
                            elseif isa(Coef,'BAscalar')
                                if length(Coef) == 1
                                    obj = BAscalar(Coef.Coef,obj.Modes);
                                else
                                    obj = Coef;
                                end
                            end
                            padCoef(obj);
                    end
                end
            end
        end
        
        function padCoef(obj)
            % pads a BAscalar with zeros to achieve truncation size consistent with obj.Modes
            if length(obj) > 1
                return
            elseif isequal(obj.Modes,size(obj.Coef))
                return
            end
            switch obj.SurfaceDimension
                case 0
                    if strcmp(obj.CoefType,'double')
                        coef = zeros(1,obj.Modes);
                    elseif strcmp(obj.CoefType,'intval')
                        coef = midrad(zeros(1,obj.Modes),0);
                    end
                    coef(1,1) = obj.Coef;
                case 1
                    if strcmp(obj.CoefType,'double')
                        coef = zeros(obj.Modes(1),obj.Modes(2));
                    elseif strcmp(obj.CoefType,'intval')
                        coef = midrad(zeros(obj.Modes(1),obj.Modes(2)),0);
                    end
                    coef(1:size(obj.Coef,1),1:size(obj.Coef,2)) = obj.Coef;
                case 2
                    if strcmp(obj.CoefType,'double')
                        coef = zeros(obj.Modes);
                    elseif strcmp(obj.CoefType,'intval')
                        coef = midrad(zeros(obj.Modes),0);
                    end
                    if obj.Modes(1) == 1
                        coef(1,1:size(obj.Coef,1),1:size(obj.Coef,2)) = obj.Coef;
                    else
                        coef(1:size(obj.Coef,1),1:size(obj.Coef,2),1:size(obj.Coef,3)) = obj.Coef;
                    end
                otherwise
                    error('padCoef not implemented for this surface dimension')
            end
            obj.Coef = coef;
        end
        
        % Have to spend more time to figure out how to vectorize the disp method properly.
        %         function disp(obj)
        %             if isequal(size(obj),[1,1]);
        %                 %                 fprintf('Coefficient Type: %s \nTruncation: %d \nSurface Dimension: %d',obj.CoefType,obj.Modes,obj.SurfaceDimension)
        %                 %                 sprintf('Truncation: %d',obj.Modes);
        %                 %                 sprintf('Surface Dimension: %d',obj.SurfaceDimension);
        %                 obj.Coef
        %             elseif size(obj,1) == 1 || size(obj,2) == 1;
        %                 %                 disp(obj(1,1).CoefType)
        %                 %                 disp(obj(1,1).Modes)
        %                 %                 disp(obj(1,1).SurfaceDimension)
        %                 for j = 1:length(obj)
        %                     obj(j).Coef
        %                 end
        %             end
        %         end
        
        function intvalObj = intval(obj)
            % returns an interval enclosure of a BAscalar with CoefType 'double'
            if length(obj) > 1
                intvalObj = repmat(BAscalar(0,obj(1).Modes),size(obj));
                for j = 1:length(obj)
                    intvalObj(j) = obj(j).intval;
                end
            else
                intvalObj = BAscalar(midrad(obj.Coef,0));
            end
        end
        
        function polyObj = intlabPoly(obj)
            % returns intlab polynomial object with intval coefficients corresponding to BAscalar.
            if obj.SurfaceDimension == 1
                [S,T] = meshgrid(0:obj.Modes(2)-1,0:obj.Modes(1)-1);
                exponObj = [reshape(S,[],1),reshape(T,[],1)];
                if strcmp(obj.CoefType,'double')
                    coefObj = reshape(midrad(obj.Coef,0),[],1);
                elseif strcmp(obj.CoefType,'intval')
                    coefObj = reshape(obj.Coef,[],1);
                end
                polyObj = polynom(coefObj,exponObj,{'s','t'});
            else
                error('intlabPoly not yet implemented for SurfaceDimension other than 1')
            end
        end
        
        function columnObj = col(obj)
            % returns Coefficients of BAsclar as a column vector (double) under the canonical isomorphism
            columnObj = reshape(obj.Coef,[],1);
        end
        
        
        %% ----------------------------------- BANACH ALGEBRA ARITHMETIC METHODS -----------------------------------
        
        
        
        
        function sumObj = plus(leftObj,rightObj)
            % Defines sums of BAscalars with BAscalars and Fscalars (double or intval).
            
            if isa(leftObj,'double') % Sum of BAscalar and (double) Fscalar (left or right)
                sumObj = BAscalar(leftObj + rightObj.Coef);
            elseif isa(rightObj,'double')
                sumObj = BAscalar(leftObj.Coef + rightObj);
                
            elseif isa(rightObj,'intval') % Sum of BAscalar and (intval) Fscalar (right only)
                sumObj = BAscalar(leftObj.Coef + rightObj);
                
            else % Sum of 2 BAscalars
                if leftObj.Modes == rightObj.Modes
                    sumObj = BAscalar(leftObj.Coef + rightObj.Coef);
                else
                    padUp = max([leftObj.Modes;rightObj.Modes]);
                    sumObj = BAscalar(leftObj.Coef,padUp) + BAscalar(rightObj.Coef,padUp);
                end
            end
        end
        
        function minusObj = uminus(obj)
            % returns negative of BAscalar
            minusObj = -1*obj;
        end
        
        function minusObj = minus(obj,rightObj)
            % Difference of BAscalar and Fscalar (intval or double) or another BAscalar.
            minusObj = obj + -rightObj;
        end
        
        function sumObj = sum(obj)
            % sums a vector/array of BAscalars
            sumObj = BAscalar(zeros(obj(1).Modes));
            for j = 1:numel(obj)
                sumObj = sumObj + obj(j);
            end
        end
        
%         function prodObj = times(L,R)
%             % THIS IS BROKEN
%             if length(L) ==1
%                 dim = sum(L.Modes>1);
%                 switch dim
%                     case 1
%                         prodObj = conv(cat(2,zeros(1,length(L.Coef)-1),L.Coef),R.Coef,'valid');
%                 end
%             else
%                 prod = arrayfun(@(j)L(j).*R(j),1:length(L),'UniformOutput',false);
%                 prodObj = BAscalar(prod);
%             end
%         end
        
        function dotProduct = dot(obj,rightFactor)
            % computes the Euclidean dot product where multiplication of individual coordinates is
            % Cauchy products
            if length(obj) ~= length(rightFactor)
                error('dot product requires vectors must be the same length')
            else
                %                 ptProductCell = arrayfun(@(j)obj(j)*rightFactor(j),1:length(obj),'UniformOutput',false);
                %                 ptProduct = cell2mat(ptProductCell);
                ptProduct(1) = obj(1)*rightFactor(1);
                for j = 2:length(obj)
                    ptProduct(j) = obj(j)*rightFactor(j);
                end
                dotProduct = ptProduct.sum;
            end
            
            
        end
        
        function operator = leftTimesOperator(obj)
            % returns a linear operator, A, whose action on any BAscalar, h, is given by convolution with obj. (i.e. for every h, A(h) = obj*h)
            if strcmp(obj.CoefType,'intval')
                warning('Lmult is untested with interval valued BAscalars')
            end
            padblock = [zeros(obj.Modes(1)-1,obj.Modes(2)-1),zeros(obj.Modes(1)-1,obj.Modes(2));zeros(obj.Modes(1),obj.Modes(2)-1),obj.Coef];
            %             Lmult(obj.Modes(1),obj.Modes(2)) = BAscalar(1:obj.Modes(1),1:obj.Modes(2)); % NO LONGER WORKS?
            Lmult(obj.Modes(1),obj.Modes(2)) = BAscalar(0,obj.Modes);
            for j = 1:obj.Modes(1)
                for k = 1:obj.Modes(2)
                    Lmult(obj.Modes(1)+1-j,obj.Modes(2)+1-k) = BAscalar(padblock(j:j+obj.Modes(1)-1,k:k+obj.Modes(2)-1));
                end
            end
            operator = BAoperator(Lmult,obj.Modes);
        end
        
        
        
        function sqrtObj = sqrt(obj)
            if strcmp(obj.CoefType,'intval')
                disp('sqrt is untested with interval valued BAscalars')
            end
            sqrt_ = zeros(obj.Modes);
            sqrt_(1,1) = sqrt(obj.Coef(1,1)); %initial condition
            if abs(obj.Coef(1,1)) < 1e-6
                disp('constant term < 1e-6. square root should not be trusted')
            end
            
            D = 2*sqrt_(1,1); %scalar denominator
            for k = 1:obj.Modes(2)-1
                sqrt_(1,k+1) = (-1/D)*(recursive_conv(sqrt_(1,1:k+1),sqrt_(1,1:k+1),'1d') - obj.Coef(1,k+1));
            end
            
            D = 2*sqrt_(1,:); %row denominator
            for k = 1:obj.Modes(1)-1;
                Xnum = -(recursive_conv(sqrt_(1:k+1,:),sqrt_(1:k+1,:),'row') - obj.Coef(k+1,:));
                sqrt_(k+1,:) = taylor_divide(Xnum,D);
            end
            sqrtObj = BAscalar(sqrt_);
        end
        
        function reObj = real(obj)
            % returns real part of coefficients
            if strcmp(obj.CoefType,'intval')
                warning('real part is untested with interval valued BAscalars')
            end
            reObj = BAscalar(real(obj.Coef));
        end
        
        function imObj = imag(obj)
            % returns imaginary part of coefficients
            if strcmp(obj.CoefType,'intval')
                warning('real part is untested with interval valued BAscalars')
            end
            imObj = BAscalar(imag(obj.Coef));
        end
        
        %% ----------------------------------- MISCELLANEOUS METHODS -----------------------------------
  
        function timeDecay = decay(obj)
            % returns norm of last Taylor coefficient (with respect to time).
            if numel(obj) > 1
                rowDecay = arrayfun(@(j)obj(j).decay,1:numel(obj));
                timeDecay = reshape(rowDecay,size(obj));
            else
                switch obj.SurfaceDimension
                    case 0
                        timeDecay = abs(obj.Coef(end));
                    case 1
                        if strcmp(obj.Weight,'ones')
                            finalTimeCoef = obj.Coef(end,:);
                            timeDecay = sum(abs(finalTimeCoef(:)));
                        else
                            finalTimeCoef = BAscalar(obj.Coef(end,:),obj.Modes(2:end),obj.Weight(2:end));
                            timeDecay = finalTimeCoef.norm();
                        end
                    case 2
                        if strcmp(obj.Weight,'ones')
                            finalTimeCoef = obj.Coef(end,:,:);
                            timeDecay = sum(abs(finalTimeCoef(:)));
                        else
                            finalTimeCoef = BAscalar(obj.Coef(end,:),obj.Modes(2:end),obj.Weight(2:end));
                            timeDecay = finalTimeCoef.norm();
                        end
                end
            end
        end
        
        function dObj_dt = dt(obj)
            % compute time derivative
            if obj.SurfaceDimension == 1 % space is 1-dim
                C = repmat((0:obj.Modes(1)-1),obj.Modes(2),1)';
                dObj_dt = BAscalar(C.*obj.Coef);
            elseif obj.SurfaceDimension ==2 % space is 2-dim
                C = zeros(obj.Modes);
                for j = 2:obj.Modes(1)
                    C(j,:,:) = (j-1)*ones(obj.Modes(2:3));
                end
                dObj_dt = BAscalar(C.*obj.Coef);
            end
        end
        
        function dObj_ds = ds(obj)
            % compute spatial derivative
            if obj.SurfaceDimension == 1
                C = repmat((1:obj.Modes(2)-1),obj.Modes(1),1);
                dObj_ds = BAscalar(C.*obj.Coef(:,2:obj.Modes(2)));
            else
                error('ds not implemented for this surface dimension')
            end
        end
        
        function int_ds = intds(obj,varargin)
            % evaluate definite or indefinite integral with respect to spatial variable
            if obj.SurfaceDimension == 1
                C = repmat(1./(1:obj.Modes(2)),obj.Modes(1),1);
                int_ds = BAscalar([zeros(obj.Modes(1),1),C.*obj.Coef]); % indefinite integral
                if (nargin > 1) % compute integral obj ds on [a,b]
                    bounds = varargin{1};
                    int_ds = int_ds.fix_space(bounds(2)) - int_ds.fix_space(bounds(1)); % evaluation definite integral
                end
            else
                error('ds not implemented for this surface dimension')
            end
        end
        
        function objNorm = norm(obj,varargin)
            % computes weighted ell-one norm for BAscalar
            if numel(obj) > 1 % vectorized norm
                objNorm = arrayfun(@(j)obj(j).norm(varargin{:}),1:numel(obj));
                objNorm = reshape(objNorm,size(obj));
            else
                if strcmp(obj.Weight,'ones')
                    if nargin == 2
                        normDim = varargin{1};
                        objNorm = sum(abs(obj.Coef),normDim);
                    else
                        objNorm = sum(abs(obj.Coef(:)));
                    end
                elseif obj.SurfaceDimension ~= 1
                    error('norm not implemented for this weight and surface dimension')
                else
                    weight_matrix = bsxfun(@(x,y)obj.Weight(1).^(x).*obj.Weight(2).^(y),0:obj.Modes(2)-1,(0:obj.Modes(1)-1)');
                    objNorm = sum(dot(weight_matrix,abs(obj.Coef)));
                end
            end
        end
        
        
        function shiftObj = shift(obj)
            % returns the value of shift(obj) where shift is the BAoperator which multiplies by t.
            shiftObj = BAscalar([zeros(1,obj.Modes(2));obj.Coef(1:obj.Modes(1)-1,:)]);
        end
        
        function coefArray = coef(obj)
            % returns a double or intval array of the coefficients of obj
            if length(obj) == 1
                coefArray = squeeze(obj.Coef);
            else
                cellCoef = arrayfun(@(j)obj(j).Coef,1:length(obj),'UniformOutput',false);
                % coefArray = cell2mat(cellCoef');
                coefArray = cellCoef;
            end
        end
        
        
        %% ----------------------------------- EVALUATION METHODS -----------------------------------
        
        function evalObj = gridEval(obj,varargin)
            switch nargin
                case 3
                    s = varargin{1};
                    t = varargin{2};
                    % evaluation is in meshgrid format. This should be changed to ndgrid and incorporated into the eval method.
                    
                    flipCoefs = fliplr(flip(obj.Coef)); %switch Coef to descending powers
                    evalSpatial = nan(length(s),obj.Modes(1));
                    for j = 1:obj.Modes(1)
                        evalSpatial(:,j) = polyval(flipCoefs(j,:),s);
                    end
                    evalObj = nan(length(t),length(s));
                    for k = 1:length(s)
                        evalObj(:,k) = polyval(evalSpatial(k,:),t);
                    end
                    
                case 4 % 2d evaluation is in ndgrid format
                    % TEST EVALUATION
                    % modes = [3,3,3]
                    % A = zeros(modes);
                    % A(1,2,1) = 1;
                    % A(1,1,2) = 1;
                    % A(2,1,1) = 1;
                    % A corresponds to P(x,y,z) = x + y + z
                    %
                    % a = BAscalar(A);
                    % s1 = [1,-2,3]
                    % s2 = [0,1]
                    % t = [2,3,4,5]
                    % checkme = a.gridEval(s1,s2,t)
                    
                    s1 = varargin{1};
                    s2 = varargin{2};
                    t = varargin{3};
                    
                    evalDims = [length(s1),length(s2),length(t)];
                    evalObj = nan(evalDims);
                    coef = zeros(length(s2),length(s1),obj.Modes(1));
                    for j = 1:obj.Modes(1)
                        pj = BAscalar(squeeze(obj.Coef(j,:,:)));
                        coef(:,:,j) = pj.gridEval(s1,s2);
                    end
                    
                    for k = 1:evalDims(1)
                        for l = 1:evalDims(2)
                            evalObj(k,l,:) = polyval(flip(squeeze(coef(l,k,:))),t);
                        end
                    end
            end
        end
        
        
        function evalObj = eval(obj,data)
            % data should be an m-by-(SurfaceDimension + 1) array of space-time coordinates (of the form data = [S,T]
            switch obj.SurfaceDimension
                case 0 % data should be an m length column vector
                    flipCoefs = flip(obj.Coef);
                    evalObj = polyval(flipCoefs,data);
                case 1 % data should be an m-by-2 array
                    flipCoefs = fliplr(flip(obj.Coef)); %switch Coeffs to descending powers
                    if size(data,2) > 2
                        data = data';
                    end
                    eval_s = nan(size(data,1),obj.Modes(1));
                    for j = 1:obj.Modes(1)
                        eval_s(:,j) = polyval(flipCoefs(j,:),data(:,1));
                    end
                    evalObj = nan(size(data,1),1);
                    for k = 1:size(data,1)
                        evalObj(k) = polyval(eval_s(k,:),data(k,2));
                    end
                otherwise
                    error('eval not implemented for this surface dimension')
            end
        end
        
        % THIS NEEDS TO BE IMPLEMENTED USING INTLAB POLYEVAL AND MERGED INTO THE EVAL METHOD.
        function evalObj = intvalEval(obj,s,t)
            % s is a vector of intervals, t is a single scalar.
            X = obj.fixtime(t);
            F = polynom(flip(X.Coef),'s');
            evalObj = polyval(F,s);
        end
        
        % THIS NEEDS TO BE MERGED INTO THE EVAL METHOD.
        function evalObj = fixSpace(obj,s)
            % collapse Taylor series onto fixed space evaluation, s
            space_vector = bsxfun(@power,s,(0:size(obj.Coef,2)-1)');
            evalObj = obj.Coef*space_vector;
        end
        
        
        %% ----------------------------------- METHOD GRAVEYARD -----------------------------------
        % These methods are 'probably' not used anymore and should be deleted if no errors arise within several weeks of removing them.
        
        % MOVED TO GRAVEYARD ON 12/21 - Does this function still get called? It doesn't look correct.
        
        
        % MOVED TO GRAVEYARD ON 12/21 - I don't think this is called anymore
        function fout = Tnorm(obj)
            fout = sum(abs(obj.Coef),2);
            % error('Uncomment Tnorm in the Method graveyard')
        end
        
        % MOVED TO GRAVEYARD ON 12/21 - I don't think this is called anymore
        function fout = Snorm(obj)
            fout = sum(abs(obj.Coef),1)';
            % error('Uncomment Snorm in the Method graveyard')
        end
        
        % MOVED TO GRAVEYARD ON 12/21 - Needs to be fixed for arbitrary dimension. Then it can replace the separate methods dt and ds? Or just leave it alone
        function dBAscalar = diff(obj,indvar)
            error('Change this call to .dt or remove diff from the Method graveyard')
            %             if indvar == 't'
            %                 dBAscalar = dt(obj);
            %             else
            %                 dBAscalar = ds(obj)
            %             end
        end
        
        % DEPRECATED ON 12/21. BAscalars with intval coefficients implemented.
        function objNorm = intval_norm(obj,varargin)
            error('intval_norm is deprecated. Convert BAscalar to intval coefficients and use the usual norm')
            % computes interval enclosure of weighted ell-one norm for BAscalar
            % if numel(obj) > 1 % vectorized norm
            % objNorm = arrayfun(@(j)obj(j).intval_norm,1:numel(obj));
            % objNorm = reshape(objNorm,size(obj));
            % else
            % if strcmp(obj.Weight,'ones')
            % intval_Coef = midrad(obj.Coef,0);
            % objNorm = sum(abs(intval_Coef(:)));
            % else
            % error('not yet implemented for other weights')
            % end
            % end
        end
        
        % DEPRECATED ON 12/21. Replace by scaleTime
        function scale_time(obj,L)
            error('scale_time has been replaced by scaleTime')
        end
        
        
    end % end methods
    
    %% STATIC METHODS
    methods(Static)
        function zarray = zeros(varargin)
            if nargin ==0
                zarray = BAscalar(0,[5,5]);
            else
                zarray = repmat(BAscalar(0,[5,5]),varargin{:});
            end
        end
    end % end static methods
end % end classdef





%         function image = evaluate(obj,s,t,varargin)
%             t_obj = obj.fixtime(t);
%             if(nargin > 3)
%                 out_type = varargin{1};
%             else
%                 out_type = 'numeric';
%             end
%
%             switch out_type
%                 case 'numeric'
%                     image = [polyval(fliplr(t_obj),s)];
%
%                 case 'interval'
%                     % return enclosure of arc at time t via interval arithmetic as m-by-4 matrix of corners
%                     sints = infsup(s(1:end-1),s(2:end));
%                     imBoxes = polyval(polynom(fliplr(t_obj{1})),sints);
%                     image = nan(numints,2);
%                     image(:,1) = inf(imBoxes)';
%                     image(:,2) = sup(imBoxes)';
%             end
%         end



