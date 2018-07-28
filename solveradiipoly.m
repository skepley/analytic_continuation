function varargout = quadRadiiPoly(coef)
% Returns a rigorous bound on smallest root of a quadratic polynomial specified (usually) as [Z2, -(1-Z1), Y0].

polyCoef = intval(coef);
radiiPoly = polynom(polyCoef,'r'); % radii polynomial
rGuess = roots(mid(polyCoef)); % numerical estimate for true roots.
if nargout ==1
    rTrue = verifypoly(radiiPoly,min(rGuess));
    errorBound = sup(rTrue);
    varargout{1} = errorBound;
else
    rTrue(1) = verifypoly(radiiPoly,min(rGuess));
    rTrue(2) = verifypoly(radiiPoly,max(rGuess));
    varargout{1} = rTrue(1);
    varargout{2} = rTrue(2);
end
end

