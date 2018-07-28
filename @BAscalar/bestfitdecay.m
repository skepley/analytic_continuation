function r  = bestfitdecay(obj)
% Returns the best approximate geometric sequence fitted to the given coefficents by linear least squares regression.

% Written by S.K. 01/2018

% 1D - If the coefs are of the form a = (A, Ar, Ar^2, ... Ar^{n-1}), then log a lies on the line y(x) = log A + log(r)*x. The best fit (log) linear line is the minimizer of the sum of squared residuals
% S = sum(r_i^2) = sum(log(a_i) - beta_0 - i*beta_1)^2. In 2D, the returned ratios minimize S = sum_j sum_i (log(a_ij) - beta_0 - j*beta_1 - i*beta_2)^2 where j is the column index, and i is the row index. 

if length(obj) > 1
    for j = 1:length(obj)
       r(j,:) = bestfitdecay(obj(j)); 
    end
else
    if obj.SurfaceDimension ==1 % coef is a 1D coefficient vector
        logCoef = log(abs(obj.Coef)); % best fit geometric sequence is linear best fit for log(sequence)
        N = length(logCoef);
        X = 0:N-1;
        M = [N,sum(X);sum(X),sum(X.^2)];
        rhs = [sum(logCoef);sum(X.*logCoef)];
        beta = M\rhs;
        % A0 = exp(beta(1)); % constant multiple for geometric sequence
        r = exp(beta(2)); % best fit ratio for geometric sequence
    elseif obj.SurfaceDimension == 2
        logCoef = log(abs(obj.Coef));
        [M,N] = size(logCoef);
        S = 0:N-1;
        T = (0:M-1);
        
        % 2D least squares regression
        X = zeros(M*N,3);
        X(:,1) = ones(M*N,1);
        X(:,2) = repmat(S',M,1);
        X(:,3) = reshape(repmat(T,N,1),[],1);
        XTX = X'*X; % get symmetric part of X
        
        Y = reshape(logCoef',[],1);
        XTY = X'*Y;
        beta = XTX\XTY;
        % A0 = exp(beta(1));
        r = exp(beta(2:3));
    else
        error('not implemented for higher dimensions')
    end
end
