function xi = icdf_transform(p, x, est_method, cdf_samples, varargin)
    

% For the back-transformation from the rank- to the data-space, we also
% % need an estimate of y's CDF at equally spaced points. In order to
% % further deal with values below and above min(y) and max(y),
% % respectively, 1/2*std(y) is added to the limits.
% if min(y) == 0
%     LB = nanmin(y);
% else
%     LB = nanmin(y) - 0.5*nanstd(y);
% end
% 
% if max(y) == 0
%     UB = nanmax(y);
% else
%     UB = nanmax(y) + 0.5*nanstd(y);
% end
% 
% x_d = linspace(LB, UB, cdf_smples);


if strcmp(est_method, 'ksd')
    xi  = ksdensity(x, p, 'function', 'icdf');
elseif strcmp(est_method, 'param')
    xi = param_transform(p, varargin{1}, 1);
elseif strcmp(est_method, 'ecdf')
    xi = smooth_ecdf(x, p, 1);
elseif strcmp(est_method, 'ranks')
    xi = rank_transform(x, o, 1);
end

    





