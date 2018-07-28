function objTailRatio = tailRatio(obj,tailBeginsAt)
% computes the norm of the higher order terms of a polynomial relative to the polynomial's norm

% Written by S. Kepley 04/2017
% Added support for intvals 08/2017

% INPUT:
% tailBeginsAt (double): vector of values specifying where the tail terms begin

objDimsidx = obj(1).dims();
Sz = size(obj(1).Coef);
Sz = Sz(objDimsidx);
if any(tailBeginsAt < 1) % Non-tail terms specified as fraction of total terms.
    fractionTerms = tailBeginsAt < 1;
    tailBeginsAt(fractionTerms) = round(tailBeginsAt(fractionTerms).*Sz(fractionTerms));
    tailBeginsAt(tailBeginsAt < 1) = ones(size(tailBeginsAt < 1));
end

if length(obj) > 1 % vectors of BAscalars
%     getRatio = @(j)obj(j).tailRatio(tailBeginsAt);
    objTailRatio(length(obj)) = obj(end).tailRatio(tailBeginsAt);
    for j = 1:length(obj)-1
       objTailRatio(j) = obj(j).tailRatio(tailBeginsAt);
    end
else
    switch length(tailBeginsAt)
        case 1 % obj is 1d
            objTail = BAscalar(obj.Coef(tailBeginsAt:end));
            
        case 2 % obj is 2d
            tailCoef = obj.Coef;
            tailCoef(1:tailBeginsAt(1)-1, 1:tailBeginsAt(2)-1) = zeros(1:tailBeginsAt(1)-1, 1:tailBeginsAt(2)-1);
            objTail = BAscalar(tailCoef);
            
        case 3 % obj is 3d
            objTail = BAscalar(obj.Coef(tailBeginsAt(1):end,tailBeginsAt(2):end,tailBeginsAt(3):end));
            
        otherwise
            error('objTailRatio not implemented for BAscalar with dimension greater than 3')
    end
    tailNorm = objTail.norm();
    objTailRatio = tailNorm./obj.norm;
end
if isa(objTailRatio,'intval')
    objTailRatio = sup(objTailRatio);
end
end
