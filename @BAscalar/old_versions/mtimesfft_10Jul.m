function productObj = mtimes(leftObj,rightObj,varargin)
%Define multiplication of BAscalars with BAscalars, Fscalars (double or intval) and BAoperators (multiplication in usual matrix/vector sense). Current implementation only supports dimension-1 surfaces.


% ---------------------- INPUT ----------------------
% leftObj BAscalar or BAoperator or double or intval: left factor
% rightObj BAscalar or BAoperator or double or intval: right factor
% varargin{1} = truncation type (string): Use to allow fast multiplication via FFT % when full products aren't needed

% ---------------------- OUTPUT ----------------------
% productObj: BAscalar corresponding to leftObj*rightObj


% Added support for 1D algebra 06/24/17
% Changed intval products to use FFT 07/06/17

if isa(leftObj,'BAoperator') % BAoperator acts on BAscalar (left only)
    switch size(obj,1)
        case 1
            Lv = leftObj.Matrix*rightObj.col;
            productObj = reshape(Lv,rightObj.Modes);
        otherwise
            if size(leftObj) ~= size(obj,1)*ones(1,2)
                error('Operator matrix and vector of BAscalars must have the same size')
            else
                % productCell = arrayfun(@(j)dot(leftObj(j,:),obj),1:size(obj,1),'UniformOutput',false);
                % productObj = cell2mat(productCell);
            end
    end
elseif isa(leftObj,'BAscalar') && isa(rightObj,'BAscalar') % Multiplication of 2 BAscalars (Cauchy product of analytic functions)
    
    % this is currently a catch all that needs to be fixed
    if isequal(leftObj.Modes,[]) || isequal(rightObj.Modes,[])
        disp('BAscalar should never have empty value for Modes')
        L = BAscalar(leftObj.Coef);
        R = BAscalar(rightObj.Coef);
        productObj = BAscalar(mtimes(L,R,varargin{:}));
        productObj.Modes = [];
        
    % one factor is a constant function
    elseif (leftObj.SurfaceDimension == 0 || rightObj.SurfaceDimension == 0)
        productObj = BAscalar(leftObj.Coef*rightObj.Coef);
        
    % products require factors of the same dimension
    elseif ~isequal(leftObj.SurfaceDimension,rightObj.SurfaceDimension)
        error('Cauchy products for factors of different dimensions is not defined yet')
        
    % left/right surfaceDimension match. Set product surfaceDimension equal
    else
        surfaceDimension = leftObj.SurfaceDimension;
    end
    
    if (nargin > 2); % If both factors have double CoefType default behavior is to truncate the product to the same size as obj. Otherwise call with varargin.
        truncSize = varargin{1}; % truncation types: {'Fixed','Recursion','Full',ArraySize}
        % 'Fixed' requires two BAscalars with double CoefType of identical size.
        % 'Recursion' gives fast convolution for computing Taylor coefficient by recursion. Returns only the (M+1)-st coefficient in 1st dimension. Note: Object returned in this case is a double, not a BAscalar.
        % 'Full' produces the full Cauchy product including higher order.
        % 'FFT' same as Fixed but uses FFT for the product. 
        % ArraySize = [M,N1,N2,...] produces products truncated to the size specified.
        
    elseif ~isa(leftObj.Coef,'intval') && ~isa(rightObj.Coef,'intval') && ~isequal(leftObj.Modes,rightObj.Modes)
        disp('BAscalars arent the same size. Using Full option.')
        truncSize = 'Full';
    else
        truncSize = 'Fixed';
    end
       
    switch surfaceDimension   
        case 1 % ---------------------------------------- 1D SURFACES ----------------------------------
            if strcmp(leftObj.CoefType,'intval') || strcmp(rightObj.CoefType,'intval')
                if strcmp(truncSize,'Recursion')
                    minMode = min([leftObj.Modes,rightObj.Modes]);
                    productObj = dot(leftObj.Coef(1:minMode),rightObj.Coef(end:-1:end - minMode + 1)); % Returns only the M^th coefficient of L*R
                else
                    % pad factors up to next power of 2
                    maxMode = leftObj.Modes + rightObj.Modes - 1; % Find degree of product.
                    nextPwrTwo = ceil(log2(maxMode));
                    padUp = 2^nextPwrTwo;
                    padLeft = BAscalar(leftObj,padUp);
                    padRight = BAscalar(rightObj,padUp);
                    
                    % compute product by intvalFFT
                    leftDFT = verifyfft(padLeft.Coef);
                    rightDFT = verifyfft(padRight.Coef);
                    productDFT = leftDFT.*rightDFT;
                    productCoef = real(verifyfft(productDFT,-1));
                    switch truncSize
                        case 'Fixed'
                            minMode = min([leftObj.Modes,rightObj.Modes]);
                            productObj = BAscalar(productCoef,minMode);
                        case 'Full'
                            productObj = BAscalar(productCoef,maxMode);
                        otherwise
                            productObj = BAscalar(productCoef,truncSize);
                    end
                end
            else
                switch truncSize
                    case 'Fixed' 
                        productSize = leftObj.Modes;
                        productCoef = conv([zeros(1,productSize - 1),leftObj.Coef],rightObj.Coef,'valid');
                        productObj = BAscalar(productCoef,productSize,surfaceDimension);
                    case 'Recursion'
                        productObj = dot(leftObj.Coef,flip(rightObj.Coef)); % Returns only the M^th coefficient of L*R
                    case 'Full'
                        productObj = BAscalar(conv(leftObj.Coef,rightObj.Coef));
                    otherwise % specify an integer truncation 
                        productCoef = conv(leftObj.Coef,rightObj.Coef); 
                        productObj = BAscalar(productCoef,truncSize);
                end
            end
        case 2 % ---------------------------------------- 2D SURFACES ----------------------------------
            if strcmp(leftObj.CoefType,'intval') || strcmp(rightObj.CoefType,'intval') % products for intval coefficients
                % pad to next power of 2
                maxMode = leftObj.Modes + rightObj.Modes - [1,1];
                nextPwrTwo = ceil(log2(maxMode));
                padUp = 2.^nextPwrTwo;
                padLeft = BAscalar(leftObj,padUp);
                padRight = BAscalar(rightObj,padUp);
                
                % compute product using intvall fft
                leftDFT = altfftn(padLeft.Coef);
                rightDFT = altfftn(padRight.Coef);
                productDFT = leftDFT.*rightDFT;
                productCoef = real(altifftn(productDFT));
                switch truncSize
                    case 'Fixed'
                        minMode = min([leftObj.Modes;rightObj.Modes]);
                        productObj = BAscalar(productCoef,minMode);
                    case 'Recursion'
                        minMode = min([leftObj.Modes;rightObj.Modes]);
                        productObj = productCoef(minMode(1),1:minMode(2)); % Returns an intval vector, not a BAscalar.
                    case 'Full'                        
                        productObj = BAscalar(productCoef,maxMode);
                    otherwise
                        productObj = BAscalar(productCoef,truncSize);
                end
                
            elseif isa(leftObj.Coef,'BAscalar') && isa(rightObj.Coef,'BAscalar')
                % compute coefficients of product
                if strcmp(truncSize,'Full')
                    minIdx = 1;
                    maxIdx = leftObj.Modes(1) + rightObj.Modes(1) - 1; % full cauchy product (in 1st variable) has deg m + n - 1.
                    subTruncSize = 'Full';
                elseif strcmp(truncSize,'Fixed') % Truncate Cauchy product to same size as factors in 1st variable and subModes in remainder.
                    minIdx = 1;
                    maxIdx = leftObj.Modes(1);
                    subTruncSize = min([leftObj.Modes(2:end),rightObj.Modes(2:end)]);
                elseif strcmp(truncSize,'Recursion')
                    minIdx = min([leftObj.Modes(1),rightObj.Modes(1)]);
                    maxIdx = minIdx;
                    subTruncSize = min([leftObj.Modes(2:end),rightObj.Modes(2:end)]);
                end
                
                % manually compute Cauchy product.
                for m = minIdx:maxIdx      
                    lowerIdx = (1 + max([0,m - rightObj.Modes(1)])):min([m,leftObj.Modes(1)]);
                    upperIdx = m + 1 - lowerIdx;
                    leftContribution = leftObj.Coef(lowerIdx);
                    rightContribution = rightObj.Coef(upperIdx);
                    productCoef(m) = dot(leftContribution,rightContribution,subTruncSize);
                end
                
                if strcmp(truncSize,'Fixed') || strcmp(truncSize,'Full')
                    productObj = BAscalar();
                    productObj.SurfaceDimension = leftObj.SurfaceDimension;
                    productObj.Coef = productCoef;
                    productModes = [productObj.Coef.Modes];
                    productObj.Modes = [maxIdx,max(productModes)];
                else
                    productObj = productCoef(m);
                end
                
            else % CoefType is double
                switch truncSize
                    case 'Fixed'
                        productCoef = conv2([zeros(leftObj.Modes-1),zeros(leftObj.Modes-[1,0]);zeros(leftObj.Modes-[0,1]),leftObj.Coef],rightObj.Coef,'valid');
                        productObj = BAscalar(productCoef,size(productCoef),surfaceDimension);
                    case 'FFT'
                        n = leftObj.Modes + rightObj.Modes - [1,1];
                        cauchyProduct = ifftn*(fftn(leftObj.Coef,n).*fftn(rightObj.Coef,n));
                        productObj = cauchyProduct(1:obj.Modes(1),1:obj.Modes(2));
                    case 'Recursion' % Object returned in this case is a double, not a BAscalar.
                        productObj = conv2([zeros(size(leftObj.Coef)-[0,1]),leftObj.Coef],rightObj.Coef,'valid');
                    case 'Full'
                        productObj = BAscalar(conv2(leftObj.Coef,rightObj.Coef));
                    otherwise % truncSize is a specfied size to truncate the product
                        if length(truncSize) == leftObj.SurfaceDimension;
                            productFull = mtimes(leftObj,rightObj,'Full');
                            productObj = BAscalar(productFull.Coef,truncSize);
                        else
                            error('Multiplication of these BAscalars is undefined')
                        end
                end
            end
        case 3 % ---------------------------------------- 3D SURFACES ----------------------------------
            switch truncSize
                case 'Fixed'
                    modes = leftObj.Modes;
                    padTo = leftObj.Modes + rightObj.Modes - [1,1,1];
                    fullCoef = ifftn(fftn(leftObj.Coef,padTo).*fftn(rightObj.Coef,padTo));
                    productObj = BAscalar(fullCoef(1:modes(1),1:modes(2),1:modes(3)));
                case 'Recursion'
                    modes = leftObj.Modes;
                    padTo = leftObj.Modes + rightObj.Modes - [1,1,1];
                    fullCoef = ifftn(fftn(leftObj.Coef,padTo).*fftn(rightObj.Coef,padTo));
                    productObj = BAscalar(squeeze(fullCoef(modes(1),1:modes(2),1:modes(3))));
            end
    end
elseif isa(leftObj,'BAscalar') % BAscalar multiplication with Fscalar (right)
    productObj = BAscalar(leftObj.Coef*rightObj,leftObj.Modes);
elseif isa(rightObj,'BAscalar')% BAscalar multiplication with Fscalar (left)
    try
        productObj = BAscalar(leftObj*rightObj.Coef,rightObj.Modes);
    catch ME
        if strcmp(ME.message, 'Undefined function ''times'' for input arguments of type ''BAscalar''.')
            disp('Left multiplication by intval is undefined. Use right multiplication only.')
        end
    end
end
end