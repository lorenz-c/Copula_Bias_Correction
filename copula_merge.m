function [s_dC, std_s_dC, theta_rl, id, varargout] = copula_merge(x, y, x_pred, settings);

if nargin < 4
    settings.cdf_est   = 'pareto_ksd'; 
    settings.calcdf    = 'separate';
    settings.predcdf   = 'old';
    settings.copnms    = {'gumbel'; 'clayton'; 'frank'; 'amh'; 'fgm'};
    settings.cdfsmpls  = 500;
    settings.condvals  = 500;
    settings.nsamples  = 500;
    settings.paramest  = 'tau';
    settings.checkAC   = false;
    settings.copulaest = 'selected';
    settings.removenan = 1;
    settings.copout    = 'maxp';
    settings.xp_output = 'mean';
end

%% --------------------- REMOVE NANS IN THE DATA --------------------------
if settings.removenan == 1
    % Keep the old x and y values for the cdf estimation
    x_old = x;
    y_old = y;
    x_pred_old = x_pred;
    
    % Set missing values in both x and y to NaN
    x(isnan(y)) = NaN;
    y(isnan(x)) = NaN;
    % Remove the NaNs from the data
    x(isnan(x)) = [];
    y(isnan(y)) = [];
    n_x_pred    = length(x_pred);
    % Save the locations with data in x_pred
    dta_ids     = find(~isnan(x_pred));
    % Remove the NaNs from x_pred
    x_pred(isnan(x_pred)) = [];
else
    x_old = x;
    y_old = y;
end

%% ------------------- FIT MARGINAL DISTRIBUTIONS -------------------------
% Generate a vector for the x-data
if strcmp(settings.calcdf, 'separate')
    % Estimate the CDF for x from the values in x only
    x_cdf = x_old;
elseif strcmp(settings.calcdf, 'merge')
    % Estimate the CDF for x from both x and x_pred (assuming that
    % both come from the same underlying distribution)
    if length(x_old) == length(x_pred_old)
        if abs(sum(x_old - x_pred_old)) == 0
            x_cdf = x_old;
        end
    else
        x_cdf = [x_old; x_pred_old];
    end
end


if strcmp(settings.predcdf, 'old')
    % Use the previously derived marginal distribution for transforming
    % the x_pred to the rank space (e.g. if x_pred is not "long" enough)
    x_p_cdf = x_old;
elseif strcmp(settings.predcdf, 'new')
    % Compute a new marginal distribution for the predictor variables
    x_p_cdf = x_pred_old;
elseif strcmp(settings.predcdf, 'merge')
    % Combine the values from x and x_pred for estimating a combined
    % marginal distribution
    if length(x_old) == length(x_pred_old)
        if abs(sum(x_old - x_pred_old)) == 0
            x_p_cdf = x_old;
        end
    else
        x_p_cdf = [x_old; x_pred_old];
    end
end

if strcmp(settings.cdf_est, 'param')
    [r, r_tmp, x_tmp, D_x]         = cdf_transform(x_cdf, x, settings.cdf_est);
    [r_p, r_p_tmp, x_p_tmp, D_x_p] = cdf_transform(x_p_cdf, x_pred, settings.cdf_est);
    [s, s_tmp, y_tmp, D_y]         = cdf_transform(y_old, y, settings.cdf_est);
    
    varargout{1}.D_x   = D_x;
    varargout{1}.D_y   = D_y;
    varargout{1}.D_x_p = D_x_p;
else
    [r, r_tmp, x_tmp]       = cdf_transform(x_cdf, x, settings.cdf_est);
    [r_p, r_p_tmp, x_p_tmp] = cdf_transform(x_p_cdf, x_pred, settings.cdf_est);
    [s, s_tmp, y_tmp]       = cdf_transform(y_old, y, settings.cdf_est);
end

% Check for 0 and 1 in the transformed data
r_p(r_p == 0) = 1e-16;
r(r == 0)     = 1e-16;
s(s == 0)     = 1e-16;

r_p(r_p == 1) = 1 - 1e-16;
r(r == 1)     = 1 - 1e-16;
s(s == 1)     = 1 - 1e-16;

%% ------------------- CHECK FOR ANTI-CORRELATION -------------------------
if settings.checkAC == true
    CC = corr(r, s, 'type', 'Kendall');
    
    if CC < 0
        r   = 1 - r;
        r_p = 1 - r_p;
        AC = true;
    else
        AC = false;
    end
end

%% ------------------------------------------------------------------------
%         Perform an initial test which Copulas might fit to the
%                    data, based on Kendall's tau
%  ------------------------------------------------------------------------
tau      = kendall(x, y, 0);
families = settings.copnms;

for i = 1:length(settings.copnms)
    boolean(i) = check_tau(families{i}, tau);
end

% Remove the Copula-families, which did not pass the alpha-test
families(boolean == 0) = [];

%% ---      Find a suitable Copula, based on the Bayesian approach      ---
p = bcs(families, [r s], [-0.86 0.86]);
p = p./sum(p);

[p_max id] = max(p);
sel_family = families{id};



%% ------------------------------------------------------------------------
%                     ESTIMATE THE COPULA PARAMETERS
%  ------------------------------------------------------------------------
% The Copula-parameter can be estimated via Canonical Maximum Likelihood
% (default) and through Kendall's tau
if strcmp(settings.paramest, 'cml')
    theta_rl = CML(sel_family, r, s, 1);
elseif strcmp(settings.paramest, 'tau')
    theta_rl = copulaparam(sel_family, tau);
end


%% ------------------------------------------------------------------------
%          COMPUTE THE CONDITIONAL COPULA AND TRANSFORM RANK DATA
%  ------------------------------------------------------------------------

% Create a vector with values between 0 and 1
v   = (linspace(0 + eps, 1 - eps, settings.condvals))';

% Using matrix-calculus speeds up the code:
r_p = repmat(r_p', size(v, 1), 1);
v   = repmat(v, 1, size(r_p, 2));

    

% for i = 1:length(families)
% 
%     % Compute the conditional Copulas
%     dC{i}   = conditionalcdf(families{i}, r_p, v, theta_rl(i, 1));
% 
%     for j = 1:size(r_p, 2)
% %        if j == 21
% %            keyboard
% %        end
%         % Draw some random samples from the conditional Copula  
%         R_dC(:, j) = cond_rnd(v, dC{i}(:, j), settings.nsamples);
%     end
% 
%     if settings.checkAC == 1 
%         if AC == true
%             R_dC  = 1 - R_dC;
%         end
%     end
%     
%     if strcmp(settings.cdf_est, 'param')
%         s_dC{i} = cdf_transform(y, R_dC, settings.cdf_est, 1, D_y);
%         
%     else
%         s_dC{i} = cdf_transform(y, R_dC, settings.cdf_est, 1);
%     end
%     
%     if settings.removenan == 1
%         s_dC_out             = NaN(settings.nsamples,  n_x_pred);
%         s_dC_out(:, dta_ids) = s_dC{i};
%         s_dC{i} = s_dC_out;
%     end
% end

% 
dC = conditionalcdf(sel_family, r_p, v, theta_rl);

for j = 1:size(r_p, 2)  
    try
        R_dC(:, j) = cond_rnd(v, dC(:, j), settings.nsamples);
    catch
        keyboard
    end
end

if settings.checkAC == 1 
    if AC == true
        R_dC  = 1 - R_dC;
    end
end

% R_dC     = mean(R_dC, 1);
% std_R_dC = std(R_dC, 1);
if strcmp(settings.cdf_est, 'param')
    s_dC = cdf_transform(y_old, R_dC, settings.cdf_est, 1, D_y);
else
    s_dC = cdf_transform(y_old, R_dC(:), settings.cdf_est, 1);
end


s_dC = reshape(s_dC, size(R_dC, 1), size(R_dC, 2));

std_s_dC = std(s_dC, 0, 1);
s_dC     = mean(s_dC);

if settings.removenan == 1
	s_dC_out          = NaN(n_x_pred, 1);
    std_s_dC_out      = NaN(n_x_pred, 1);
    
    s_dC_out(dta_ids) = s_dC;
    std_s_dC_out(dta_ids) = std_s_dC;
    
    s_dC              = s_dC_out;
    std_s_dC          = std_s_dC_out;
end
    
    
    
    
    
 
 
 
 







