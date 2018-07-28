function [gamma,maxError] = liftsegment(segmentFrom,segmentTo,P,localError)
% Lift the parameterized linear segment between two points through P with rigorous error.

% Written by S.K. 03/2018

% INPUT: 
% segmentFrom, segmentTo: Ordered pairs representing points in C which are the endpoints of the segment. 
% P: M-by-N-by-d array of coefficients for a 2-dimensional local manifold embedded in R^d. 
% localError: A vector of the analytic error bounds for each coordinate of P. 

% OUTPUT: 
% gamma: Taylor coefficients for the image of the semgent lifted through P. Each row is one, 1-dimensional Taylor series. 
% maxError: The error for the lifted arc in the max-norm on R^d including all shrinkwrapping and truncation errors. 


[M,N,numVar] = size(P);
deg = max([M,N]) -1;

% parameterized line segment
thisParm(:,1) = .5*(segmentFrom + segmentTo);
thisParm(:,2) = .5*(segmentTo - segmentFrom);
intvalParm = intval(thisParm);

% lift boundary through P
intvalLift = intval(zeros(numVar,deg +1));
for j = 1:numVar
    intvalLift(j,:) = liftparm(P(:,:,j),intvalParm(1,:),intvalParm(2,:));
end

% shrink-wrap lifted boundary
gamma = zeros(numVar,deg +1); % coefficients for lifted boundary
liftError = zeros(1,numVar);
for j = 1:numVar
    [gamma(j,:),liftError(j)] = shrinkwrap(intvalLift(j,:),N,localError);
end
maxError = max(liftError);
end

