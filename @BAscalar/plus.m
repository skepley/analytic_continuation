function sumObj = plus(leftObj,rightObj)
% Defines sums of BAscalars with BAscalars and Fscalars (double or intval).
	
% Written S. Kepley 05/2017
% Updated for BAscalar coefficients 07/2017

% ---------------------- INPUT ----------------------
% leftObj (BAscalar)
% rightObj (BAscalar)

% ---------------------- OUTPUT ----------------------
% sumObj (BAscalar) - sum of the input objects in the Banach algebra. 


	
if isa(leftObj,'double') % Sum of BAscalar and (double) Fscalar (left or right)
    sumObj = BAscalar(rightObj);
    if isa(sumObj.Coef,'BAscalar')
        sumObj.Coef(1).Coef(1) = sumObj.Coef(1).Coef(1) + leftObj;
    else
        sumObj.Coef(1) = sumObj.Coef(1) + leftObj;
    end
elseif isa(rightObj,'double') || isa(rightObj,'intval');  % Sum of BAscalar and (double or intval) Fscalar (right only)
    sumObj = BAscalar(leftObj);
    if isa(sumObj.Coef,'BAscalar')
        sumObj.Coef(1).Coef(1) = sumObj.Coef(1).Coef(1) + rightObj;
    else
        sumObj.Coef(1) = sumObj.Coef(1) + rightObj;
    end

else % Sum of 2 BAscalars
    if ~isequal(leftObj.SurfaceDimension,rightObj.SurfaceDimension)
       error('Addition for different surface Dimensions is not defined') 
    end
    
    if ~isa(leftObj.Coef,'BAscalar') && ~isa(rightObj.Coef,'BAscalar') % both summands have double or intval Coefs
        if leftObj.Modes == rightObj.Modes
            sumObj = BAscalar(leftObj.Coef + rightObj.Coef,leftObj.Modes);
        else % summands have non-equal size.
            padUp = max([leftObj.Modes;rightObj.Modes]);
            sumObj = BAscalar(leftObj.Coef,padUp) + BAscalar(rightObj.Coef,padUp);
        end
        
    elseif isa(leftObj.Coef,'BAscalar') && isa(rightObj.Coef,'BAscalar') % both summands have BAscalar Coefs
        if isequal(leftObj.Modes,rightObj.Modes) % summands have same size.
            sumObj = BAscalar(leftObj);
            for j = 1:leftObj.Modes{1}
                sumObj.Coef(j) = sumObj.Coef(j) + rightObj.Coef(j);
            end
        else % summands have non-equal size.
            commonModes = min([leftObj.Modes{1};rightObj.Modes{1}]);
            if leftObj.Modes{1} > commonModes % left summand is the large one.
                largerSummand = leftObj;
            else
                largerSummand = rightObj;
            end
            
            sumObj = BAscalar(largerSummand);
            for j = 1:commonModes
                sumObj.Coef(j) = sumObj.Coef(j) + rightObj.Coef(j); % sum over non-zero modes.
            end
        end
    else
        error('Addition of BAscalar CoefType and non-BAscalar CoefType not supported')
    end
end
end