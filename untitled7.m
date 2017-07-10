function [r, s, r_p, s_d, s_d_u, y_d_u] = rank_transform(x, y, x_pred, est_method, cdf_smples, estcdf, predcdf);

if nargin < 7, predcdf    = 'new'; end
if nargin < 6, calcdf     = 'separate'; end
if nargin < 5, cdf_smples = 1000; end
if nargin < 4, method     = 'ranks'; end



% For the back-transformation from the rank- to the data-space, we also
% need an estimate of y's CDF at equally spaced points. In order to
% further deal with values below and above min(y) and max(y),
% respectively, 1/2*std(y) is added to the limits.
if min(y) == 0
    LB = nanmin(y);
else
    LB = nanmin(y) - 0.5*nanstd(y);
end

if max(y) == 0
    UB = nanmax(y);
else
    UB = nanmax(y) + 0.5*nanstd(y);
end

y_d = linspace(LB, UB, cdf_smples);
    

if strcmp(est_method, 'ksd')
    % Transform data to the rank space using ksdensity
    if strcmp(calcdf, 'separate')
        % Estimate the CDF for x from the values in x only
        r  = ksdensity(x, x, 'function', 'cdf');
    elseif strcmp(calcdf, 'merge')
        % Estimate the CDF for x from both x and x_pred (assuming that
        % both come from the same underlying distribution)
        r  = ksdensity([x; x_pred], x, 'function', 'cdf'); 
    end
    
    if strcmp(predcdf, 'new')
        % Compute a new marginal distribution for the predictor variables
        r_p = ksdensity(x_pred, x_pred, 'function', 'cdf');
    elseif strcmp(predcdf, 'old')
        % Use the previously derived marginal distribution for transforming
        % the x_pred to the rank space
        r_p = ksdensity(x, x_pred, 'function', 'cdf');
    elseif strcmp(predcdf, 'merge')
        % Combine the values from x and x_pred for estimating a combined
        % marginal distribution
        r_p = ksdensity([x; x_pred], x_pred, 'function', 'cdf');
    end

    % Transform the y-data to the rank-space
    s   = ksdensity(y, y, 'function', 'cdf');
    % Transform the data for the back-transformation ot the rank space
    s_d = ksdensity(y, y_d, 'function', 'cdf');

    % Check if s_d is strictly monotonic increasing
    [s_d_u, ia, ic] = unique(s_d);
    y_d_u           = y_d(ia);
    
    % Clean up
    clear ia ic LB UB 
    
elseif strcmp(est_method, 'param')
    % Transform data to the rank space using ksdensity
    if strcmp(calcdf, 'separate')
        % Estimate the CDF for x from the values in x only
        [D_x PD_x] = allfitdist(x);
        r          = param_transform(x, D_x);
    elseif strcmp(calcdf, 'merge')
        % Estimate the CDF for x from both x and x_pred (assuming that
        % both come from the same underlying distribution)
        [D_x PD_x] = allfitdist([x; x_pred]);
        r          = param_transform([x; x_pred], D_x);
    end
    
    
    if strcmp(predcdf, 'new')
        % Compute a new marginal distribution for the predictor variables
        [D_x_p PD_x_p] = allfitdist(x_pred);
        r_p            = param_transform(x_pred, D_x_p);
    elseif strcmp(predcdf, 'old')
        % Use the previously derived marginal distribution for transforming
        % the x_pred to the rank space
        [D_x_p PD_x_p] = allfitdist(x);
        r_p            = param_transform(x_pred, D_x_p);
    elseif strcmp(predcdf, 'merge')
        % Combine the values from x and x_pred for estimating a combined
        % marginal distribution
        [D_x_p PD_x_p] = allfitdist([x; x_pred]);
        r_p            = param_transform([x; x_pred], D_x_p);
    end

    % Fit a distribution to the y-data
    [D_y PD_y] = allfitdist(y);
    % Transform the y-data to the rank-space
    s          = param_transform(y, D_y);
    % Transform the data for the back-transformation ot the rank space
    s_d        = param_transform(y_d, D_y);

    % Check if s_d is strictly monotonic increasing
    [s_d_u, ia, ic] = unique(s_d);
    y_d_u           = y_d(ia);
    
    % Clean up
    clear ia ic LB UB
    
elseif strcmp(est_method, 'ecdf')
    
    if strcmp(calcdf, 'separate')
        % Estimate the CDF for x from the values in x only
        [r_tmp, x_tmp]  = ecdf(x);
        r               = interp1(x_tmp(2:end), r_tmp(2:end), x);
    elseif strcmp(calcdf, 'merge')
        % Estimate the CDF for x from both x and x_pred (assuming that
        % both come from the same underlying distribution)
        r = ecdf([x; x_pred], x)
        r = interp1([x; x_pred], r, x);
    end
    [s_tmp, y_tmp]  = ecdf(y);
    s               = interp1(y_tmp(2:end), s_tmp(2:end), y);

    if strcmp(predcdf, 'new')
        % Compute a new marginal distribution for the predictor variables
        [r_p_tmp, x_pred_tmp] = ecdf(x_pred);
        r_p                   = interp1(x_pred_tmp(2:end), ...
                              r_p_tmp(2:end), x_pred,  'linear', 'extrap');   
    elseif strcmp(predcdf, 'old')
        % Use the previously derived marginal distribution for transforming
        % the x_pred to the rank space
        [r_p_tmp, x_tmp] = ecdf(x);
        r_p              = interp1(x_tmp(2:end), r_p_tmp(2:end), ...
                                              x_pred,  'linear', 'extrap');   
    elseif strcmp(predcdf, 'merge')
        % Combine the values from x and x_pred for estimating a combined
        % marginal distribution
        [r_p_tmp, x_tmp] = ecdf([x; x_pred]);
        r_p              = interp1(x_tmp(2:end), r_p_tmp(2:end), ...
                                               x_pred, 'linear', 'extrap');   
    end
    
    [s_d_tmp, y_tmp] = ecdf(y);
    s_d              = interp1(y_tmp(2:end), s_d_tmp(2:end), y_d, ...
                                                       'linear', 'extrap');

    % Check if s_d is strictly monotonic increasing
    [s_d_u, ia, ic] = unique(s_d);
    y_d_u           = y_d(ia);
    
    % Clean up
    clear ia ic LB UB
    
elseif strcmp(est_method, 'ranks')
    
    if strcmp(calcdf, 'separate')
        % Estimate the CDF for x from the values in x only
        for i = 1:length(x)
            tmp  = find(x <= x(i));
            r(i) = length(tmp)/(length(x) + 1);
        end
    elseif strcmp(calcdf, 'merge')
        % Estimate the CDF for x from both x and x_pred (assuming that
        % both come from the same underlying distribution)
        for i = 1:length(x)
            tmp  = find([x; x_pred] <= x(i));
            r(i) = length(tmp)/(length([x; x_pred]) + 1);
        end
    end
    
    
    if strcmp(predcdf, 'new')
        % Compute a new marginal distribution for the predictor variables
        for i = 1:length(x_pred)
            tmp  = find(x_pred <= x_pred(i));
            r_p(i) = length(tmp)/(length(x_pred) + 1);
        end
    elseif strcmp(predcdf, 'old')
        % Use the previously derived marginal distribution for transforming
        % the x_pred to the rank space
        for i = 1:length(x_pred)
            tmp  = find(x <= x_pred(i));
            r_p(i) = length(tmp)/(length(x) + 1);
        end
    elseif strcmp(predcdf, 'merge')
        % Combine the values from x and x_pred for estimating a combined
        % marginal distribution
        for i = 1:length(x_pred)
            tmp  = find([x; x_pred] <= x_pred(i));
            r_p(i) = length(tmp)/(length([x; x_pred]) + 1);
        end
     end
    
        
    for i = 1:length(x)
        tmp = find(y <= y(i));
        s(i) = length(tmp)/(length(y) + 1);
    end
    
    for i = 1:length(y_d)
        tmp = find(y,<= y_d(i))
        s_d(i) = length(tmp)/(length(y) + 1);
    end
                                          
    % Check if s_d is strictly monotonic increasing
    [s_d_u, ia, ic] = unique(s_d);
    y_d_u           = y_d(ia);
    
    % Clean up
    clear ia ic LB UB
    
    
end



if plotflag == 1
    figure
    subplot(1,2,1)
    scatterhist(x, y, 'nbins', floor(length(x)/10));
    subplot(1,2,2)
    scatterhist(r, s, 'nbins', floor(length(r)/10));
end