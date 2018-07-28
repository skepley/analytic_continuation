classdef BAoperator
    
	
	% Written by S.K. 05/2016
    
    % FIXED: 
	% Removed automatic construction of the Tensor object. This is slow and I don't think its used anymore. 09/2016
	% Fast operator inversion added 02/2017
	
	

    properties
        Tensor;
        Matrix;
        Modes;
    end
    
    properties(Hidden = 1)
        Weight='ones';
    end
    
    
    methods
        function obj = BAoperator(action,Modes)
            % class constructor
            % defines a linear operator in BA_MN (The bi-infinite Banach algebra truncated to [M,N] Modes) in one of two ways.
            % Action can be a canonical Matrix description of the linear operator
            % after vectors in BA are converted to columns. In this case
            % action is specified by a (M*N)-by-(M*N) Matrix.
            %
            % Alternatively, action is an array of dimension (M,N,M,N)
            % where (i,j,:,:) is the image of the ij^th basis vector untder
            % the linear transformation.
            if(nargin > 0)
                obj.Modes = Modes;
                if(isa(action,'double')) %INPUT IS A MATRIX
                    obj.Matrix = action;
                    obj.Tensor = [];
                    
                        % DON'T NEED THIS PART UNLESS .TENSOR IS CALLED
                        % SOMEWHERE WHICH SO FAR IT ISN'T
%                     blocks = reshape(action',[Modes,Modes]);
%                     obj.Tensor = cell(Modes(1),Modes(2));
%                     for m = 1:Modes(1)
%                         for n = 1:Modes(2)
%                             obj.Tensor{m,n} = squeeze(blocks(m,n,:,:));
%                         end
%                     end

                else
                    obj.Tensor = action; %INPUT IS AN ARRAY OF BASCALARS
                    M = zeros([Modes,Modes]);
                    for m = 1:Modes(1)
                        for n = 1:Modes(2)
                            M(m,n,:,:) = action(m,n).Coef;
                        end
                    end
                    obj.Matrix = reshape(M,prod(Modes)*ones(1,2))';
                end
            end
            
        end %class constructor
        
                
        function opNorm = norm(obj)
%             [idRow,idCol] = ndgrid(1:obj.Modes(1),1:obj.Modes(2));
%             norm_matrix = arrayfun(@(j,k)obj.Tensor(j,k).norm,idRow,idCol);
%             norm_matrix = arrayfun(@(j,k)sum(abs(obj.Tensor{j,k}(:))),idRow,idCol);
% NEED TO MAKE SURE THAT MATRIX NORM BOUNDS THE TENSOR NORM

			switch size(obj,1)
				case 1
					if isequal(obj.Weight,'ones')
		%                 opNorm = max(norm_matrix(:));
						opNorm = norm(obj.Matrix,Inf);
					else
						weight_matrix = bsxfun(@(x,y)obj.Weight(1).^(x).*obj.Weight(2).^(y),0:obj.Modes(2)-1,(0:obj.Modes(1)-1)');
						opNorm = max(max(norm_matrix./weight_matrix));
					end
				otherwise
					coordNorms = arrayfun(@(i,j)obj(i,j).norm,size(obj,1),size(obj,2));
					rowNorms = sum(coordNorms,2);
					opNorm = max(rowNorms);
			end 
        end
        
        function opNorm = intval_norm(obj)
            intval_coef = midrad(obj.Matrix',0);
            if isequal(obj.Weight,'ones')
                opNorm = norm(intval_coef,Inf);
            else
                error('not yet implemented for other weights')
            end
        end
            
        
        function op_v = mtimes(obj,v)
            if isa(obj,'double')
                op_v = BAoperator(obj.*v.Matrix,v.Modes);
            elseif isa(v,'double')
                op_v = BAoperator(v.*obj.Matrix,obj.Modes);
            elseif isa(v,'BAscalar')
                if isequal(size(v),[1,1])
                    Lv = obj.Matrix*v.col;
                    op_v = BAscalar(reshape(Lv,v.Modes));
                elseif size(obj,1) ~= length(v)
                    error('Operator matrix and vector of BAscalars must have the same size')
                else
                    for j = 1:size(obj,1)
                        rowJ(1) = obj(j,1)*v(1);
                        for k = 2:size(obj,2)
                            rowJ(k) = obj(j,k)*v(k);
                        end
                        op_v(j) = sum(rowJ);
                    end
                    op_v = reshape(op_v,size(v));
                end
            elseif isa(v,'BAoperator')
                if isequal(obj.Modes,v.Modes)
                    op_v = BAoperator(obj.Matrix*v.Matrix,obj.Modes);
                else
                    error('Tensors must have the same covariant and contravariant rank')
                end
            end
        end % end mtimes

                
                                

        
        function op_sum = plus(obj,v)
            if isa(v,'double') && isequal(size(obj.Matrix),size(v))
                op_sum = BAoperator(obj.Matrix + v);
            elseif isa(v,'BAoperator') && isequal(obj.Modes,v.Modes)
                op_sum = BAoperator(obj.Matrix + v.Matrix,obj.Modes);
            else
                error('dimension mismatch')
            end
        end
        
        function op_minus = minus(obj,v)
            if isa(v,'double') && isequal(size(obj.Matrix),size(v))
                op_minus = BAoperator(obj.Matrix - v);
            elseif isa(v,'BAoperator') && isequal(obj.Modes,v.Modes)
                op_minus = BAoperator(obj.Matrix - v.Matrix,obj.Modes);
            else
                error('dimension mismatch')
            end
        end
        
        function opBlock = block(obj)
             % returns a block matrix for operator specified as a square array of BAscalars. s 
            opDims = size(obj);
            opCell = cell(opDims);
            for j = 1:prod(opDims)
                opCell{j} = obj(j).Matrix;
            end
            opBlock = cell2mat(opCell);
        end

        function invBlock = inv(obj)
            % inverts a square linear operator on product space of BAscalars (assumed to be in Taylor basis)
            opDims = size(obj);
            if opDims(1) ~= opDims(2)
                error('BAoperator array must be square')
            end
            opCell = cell(opDims);
            for j = 1:prod(opDims)
                opCell{j} = obj(j).Matrix;
            end
%             opBlock = cell2mat(opCell);
            opTril = block2tril(opCell); % convert to lower triangular (triL) basis
            invTril = opTril\eye(opDims(1)*prod(obj(1).Modes)); % fast inversion using triL structure
            invBlock = tril2block(invTril,opDims); % convert inverse back to block basis
        end
            
        
   
    end % end methods
end % end classdef

