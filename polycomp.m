function powerArray = polycomp(p,m)
% Given coefficients defining a single variable polynomial, returns coefficients for powers of p. 

% Written by S.K. 03/2015
% Added support for intLab polynomials 06/2017 

% INPUT: 
% p: vector of coefficients defining p in INCREASING order. 
% m: non-negative integer specifying largest power of p to compute. 

% OUTPUT: 
% powerArray: a matrix whose (k+1)^th row is the coefficients for p(s)^k in INCREASING order

deg = length(p) - 1;
if isa(p,'intval')
    P = polynom(fliplr(p)); % intLab has coefficients in decreasing order
    powerArray = cell(m,1);
    powerArray{1} = polynom(1);
    for j = 2:m+1
        powerArray{j} = P*powerArray{j-1};
    end
else
    M = zeros(m,deg*m + 1);
    M(1,1:deg+1) = p;
    for j = 2:m
        row = conv(p,M(j-1,:));
        M(j,:) = row(1:size(M,2));
    end
    powerArray = [zeros(1,size(M,2));M];
    powerArray(1,1) = 1;
end
end