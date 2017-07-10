function tau = copulastat(family,alpha)
%
%   TAU = COPULASTAT(FAMILY,ALPHA)
%
%   Returns Kendall's rank correlation TAU
%
%   INPUT
%       FAMILY: One of 'Clayton', 'Frank', 'Gumbel', 'Gaussian', 't',
%               'AMH', 'FGM', 'Arch12', 'Arch14'.
%       ALPHA: Copula parameter.
%
%   OUTPUT
%       TAU: Kendall's rank correlation.
%
%   Example:
%       Determine the theoretical rank correlation coefficient for a
%       bivariate Gaussian copula with linear correlation parameter -0.7071
%       rho = -.7071
%       tau = copulastat('gaussian',rho)

%   Original version by Peter Perkins, The MathWorks, Inc.
%   Revision: 1.0  Date: 2003/09/05
%
%   Extended and modified by G.Evin, 2005.
%   Modified by D. Huard, Nov. 2006.
%   Requires MATLAB R13.

% -------------------------------------------------------------------------
% Update: Jun 2016: Added support for Joe Copula (C. Lorenz)
% -------------------------------------------------------------------------

if nargin ~= 2
    error('Requires two input arguments.');
end

% pass = check_alpha(family, alpha, 0);
% 
% if any(~pass)
%     error('Invalid parameters.')
% end

tau = zeros(size(alpha));

switch lower(family)
    case 'gaussian' 
        tau = 2.*asin(alpha)./pi;
        
    case 't'
        tau = 2.*asin(alpha)./pi;
        
    case 'clayton'
        tau = alpha ./ (2 + alpha);
        
    case 'frank'    
        tau = real(1 - (4*(1 - (-3*alpha.^2 - pi^2 + 6*alpha.*log(1 - exp(alpha)) + 6*dilog(exp(alpha)))./(6*alpha)))./alpha);
  %      fun = @(t) t./(exp(t)-1);
  %      D   = 1./alpha*integral(fun, 0, alpha);
  %      tau = 1-4./alpha.*(1 - D);
        
    case 'gumbel'
        tau = 1 - 1./alpha;
        
    case 'amh'
        t0 = alpha.^2;
        t1 = log(1-alpha);
        t2 =  t0 .* t1;
        t3 = 2*alpha .* t1;
        tau = 1- 2/3 * (t2 - t3 + alpha + t1) ./ t0;
        
    case 'fgm'
        tau = (2/9).*alpha;
        
    case 'arch12'
        tau = 1-(2/3)./alpha;
        
    case 'arch14'
        tau = 1-2./(1+2.*alpha);
    case 'joe'
%         fun = @(t) (log(1-t.^alpha).*(1-t.^alpha))./(t.^(alpha-1));
%         D   = integral(fun, 0, 1);
%         tau = 1+4./alpha.*D;

        summ  = 0;
        summP = Inf;
        k     = 1;
        while abs(summP - summ) > 1e-16
            summP = summ;
            summ  = summ + 1./(k.*(alpha*k+2).*(alpha.*(k-1)+2));
            k = k + 1;
        end
        tau = 1 - 4*summ;
            
            
    otherwise
        error('Unrecognized copula family: ''%s''',family);
end
