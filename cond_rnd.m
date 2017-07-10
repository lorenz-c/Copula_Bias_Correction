function R_c = cond_rnd(uv, dC, n)
% Due to the precision of doubles in Matlab, it "might" be that the
% conditional Copula is not strictly monotonic increasing. Therefore, we'll
% make sure that there are only "unique" (in terms of Matlab's precision)
% numbers, so that we can run interp1.
[dC_un, ia, ic] = unique(dC);
uv_un           = uv(ia);

% Draw n uniformly distributed random samples between 0 and 1
R_u = rand(n, 1);

% Transform the random samples so that they match the conditional copula
% CDF:
R_c = interp1(dC_un, uv_un, R_u, 'linear', 'extrap');

R_c(R_c >= 1) = 1 - eps;
R_c(R_c <= 0) = 0 + eps;
    
