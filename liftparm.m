function boundaryTaylorCoef = liftparm(P,theta1,theta2)
% Lift a single coordinate for a variable in parameter space through the local parameterized manifold into phase space. 

% Written by S.K. 05/2016

% INPUT:
% P: 2-dimensional array of Taylor coefficients defining an analytic function on the unit polydisc in C^2. So P: C^2 ---> R
% theta1,theta2: Taylor coefficients for coordinates which parameterize a manifold in the model space, Theta: R ---> C^2. 

% OUTPUT: 
% boundaryTaylorCoef: Taylor coefficients for the map, gamma: R ---> R with gamma = P o theta. 

[M,N] = size(P); % P has total degree (M-1)*(N-1)
T1 = polycomp(theta1,M-1); % matrix of powers for theta1
T2 = polycomp(theta2,N-1); % matrix of powers for theta2
deg = max([(M-1)*(length(theta1) -1), (N-1)*(length(theta2)-1)]); % Highest order term for gamma  is either (M-1)*deg(theta1) or (N-1)*deg(theta2)

if isa(T1,'cell') % interval coefficients
    boundaryTaylorCoef = intval(zeros(1,deg+1));
    for j=1:M
        for k = 1:N
            addPoly = P(j,k)*T1{k}*T2{j};
            addto = flip(addPoly.c);
            boundaryTaylorCoef(1:length(addto)) = boundaryTaylorCoef(1:length(addto)) + addto;
        end
    end
else % double coefficients  
    boundaryTaylorCoef = zeros(1,deg+1);   
    if isa(P,'intval')
        boundaryTaylorCoef = intval(boundaryTaylorCoef);
    end    
    for j=1:M
        for k = 1:N
            addto = P(j,k)*conv(T1(k,:),T2(j,:));
            addto = addto(1:min([deg+1,length(addto)])); % throw away extra zeros from conv.
            boundaryTaylorCoef(1:length(addto)) = boundaryTaylorCoef(1:length(addto)) + addto;
        end
    end
end
end