function p = posterior(alpha, family, U, prior_tau)
%
%   P = POSTERIOR(ALPHA, FAMILY, U, [,PRIOR_TAU])
%
%   Return the posterior density computed at ALPHA.
%   \prod c(u,v|ALPHA) * prior(ALPHA|TAU) * prior(TAU) 
%
%   INPUTS
%       FAMILY   : one of { 'arch12', 'arch14' 'ind', 'fgm' 'gaussian', 'gumbel' 'clayton' 'frank' 'amh'}       
%       U        : Nx2 matrix of quantiles (u,v) in [0,1]^2.
%       ALPHA    : 1xM vector of copula parameters.
%       PRIOR_TAU: Function of TAU returning the normalized prior for TAU.
%                  
%
%   OUTPUT
%       L: Likelihood at each parameter (1xM). 
%

%   G. Evin & D. Huard, 2006

%   Compute density
% if size(alpha, 1) > 1
%     alpha = alpha';
% end

c = copulapdf(family, U, alpha);
if min(min(c)) < 0 && strcmp(family, 'frank')
    warning('Negative values in Frank-PDF. Replaced with abs(c(c<0))!ARSCH')
    keyboard
    c(c < 0) = abs(c(c < 0));   
end

if any(isinf(c))
    error('Inf in copulapdf.')
end
% Compute data likelihood
likelihood = sum(log(c), 1);

% Compute prior for the parameter assuming an uniform prior on TAU.
prior_alpha = log(taujacobian(family, alpha));

% Prior for Tau
pr_tau = log(prior_tau(copulastat(family, alpha)));

% Combine likelihood and priors
p                 = exp(likelihood + prior_alpha + pr_tau);
p(abs(p) < 1e-50) = 0;

p(isinf(p))       = realmax;

% if strcmp(family, 'frank')
%     keyboard
% end
if any(isnan(p))
    keyboard
    error('Nan')
end


