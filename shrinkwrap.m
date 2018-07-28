function [doubleCoef,newError] = shrinkwrap(intvalCoef,truncationSize,initError)
% Analytic shrink-wrapping for truncated analytic functions. 

% Written by S.K. 11/2016

% INPUT: 
% intvalCoef: a vector or 2-dimensional array of intvals which define a 1 or 2 variable polynomial approximation of an analytic function.
% truncationSize: positive integer describing the truncation size for the output polynomial. 
% initError: double which is a bound for the analytic (ell^1) tail of the function approximated by intvalCoef. 

% OUTPUT: 
% doubleCoef: double representation for the polynomial specified by intvalCoef (after truncation)
% newError: Updated estimate for the analytic tail which includes truncation, initial error, and the "slack" from the non-truncated intval coefficients. 

switch length(truncationSize)
    case 1 % 1-dimensional polynomial
        intvalError = sum(rad(intvalCoef(1:truncationSize))); % shrink-wrap first 1:truncationSize coefficients and replace by floats.
        absTail = abs(intvalCoef(truncationSize+1:end));
        truncError = sum(sup(absTail)); % error from terms which are truncated. 
        doubleCoef = mid(intvalCoef(1:truncationSize));
        newError = intvalError + truncError + initError;
        if isa(newError,'intval')
            newError = sup(newError);
        end
    case 2
        intvalError = sum(sum(rad(intvalCoef(1:truncationSize(1),1:truncationSize(2))))); % shrink-wrap first 1:truncationSize coefficients and replace by floats.
        absTail = abs(intvalCoef(truncationSize(1)+1:end,truncationSize(2)+1:end));
        truncError = sum(sum(sup(absTail))); % error from terms which are truncated. 
        doubleCoef = mid(intvalCoef(1:truncationSize(1),1:truncationSize(2)));
        newError = intvalError + truncError + initError;
        if isa(newError,'intval')
            newError = sup(newError);
        end
end
end

