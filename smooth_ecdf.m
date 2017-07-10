function [Fx, F, X] = smooth_ecdf(x, x_p, inverse)
% The code is adapted from 
% http://de.mathworks.com/help/stats/examples/nonparametric-estimates-of-
% cumulative-distribution-functions-and-their-inverses.html

if nargin < 3, inverse = 0; end
if nargin < 2, x_p = x; end

% Estimate the ecdf
[F, X] = ecdf(x);

% Calculate the length of the cdf-vector
n = length(F);

% Smooth the empirical cdf
F = (F(1:end-1)+F(2:end))/2;

% Throw away the first element of X
X = X(2:end);
n = length(X);
% Make sure that the start- and endpoints reach the 0 and 1
X  = [X(1)-F(1)*(X(2)-X(1))/((F(2)-F(1)));
      X;
      X(n)+(1-F(n))*((X(n)-X(n-1))/(F(n)-F(n-1)))];
         
F = [0; F; 1];

if inverse == 0
    % Calculate the cdf-values of x
    Fx = interp1(X, F, x_p, 'linear', 'extrap');
elseif inverse == 1
    Fx = interp1(F, X, x_p, 'linear','extrap');
end
