function [Fxp, Fx, X, varargout] = cdf_transform(x, x_p, est_method, inverse, varargin);

if nargin < 4, inverse = 0; end
if nargin < 3, est_method = 'ecdf'; end

warning('off', 'all');

if inverse == 0
    if strcmp(est_method, 'ksd')
        Fxp = ksdensity(x, x_p, 'function', 'cdf', 'support', 'positive');
        Fx  = ksdensity(x, x, 'function', 'cdf', 'support', 'positive');
        X   = x;
    elseif strcmp(est_method, 'param')
        [D PD]       = allfitdist(x);
        varargout{1} = D;
        Fxp          = param_transform(x_p, D);
        Fx           = param_transform(x, D);
        X            = x;       
    elseif strcmp(est_method, 'ecdf')
        [Fxp, Fx, X] = smooth_ecdf(x, x_p);
    elseif strcmp(est_method, 'ranks')
        [Fxp, Fx, X] = rank_transform(x, x_p);
    elseif strcmp(est_method, 'pareto_ecdf')
        if isempty(varargin)
            LB = 0.1; UB = 0.9;
        else
            LB = varargin{1}; UB = varargin{2};
        end
        tmp = paretotails(x, LB, UB, 'ecdf');
        Fxp = cdf(tmp, x_p);
        Fx  = cdf(tmp, x);
        X   = x;
    elseif strcmp(est_method, 'pareto_ksd')
        if isempty(varargin)
            LB = 0.1; UB = 0.9;
        else
            LB = varargin{1}; UB = varargin{2};
        end
        tmp = paretotails(x, LB, UB, 'kernel');
        Fxp = cdf(tmp, x_p);
        Fx  = cdf(tmp, x);
        X   = x;
    end
    
    [Fx, ids] = sort(Fx, 'ascend');
    X         = X(ids);
                
elseif inverse == 1
    if strcmp(est_method, 'ksd')
        Fxp  = ksdensity(x, x_p(:), 'function', 'icdf');
        Fxp  = reshape(Fxp, size(x_p));
    elseif strcmp(est_method, 'param')
        [D PD]       = allfitdist(x);
        varargout{1} = varargin{1};
        Fxp          = param_transform(x_p, varargin{1}, 1);
    elseif strcmp(est_method, 'ecdf')
        Fxp = smooth_ecdf(x, x_p, 1);
    elseif strcmp(est_method, 'ranks')
        Fxp = rank_transform(x, x_p, 1);
        elseif strcmp(est_method, 'pareto_ecdf')
        if isempty(varargin)
            LB = 0.1; UB = 0.9;
        else
            LB = varargin{1}; UB = varargin{2};
        end
        tmp = paretotails(x, LB, UB, 'ecdf');
        Fxp = icdf(tmp, x_p);
        Fx  = icdf(tmp, x);
        X   = x;
    elseif strcmp(est_method, 'pareto_ksd')
        if isempty(varargin)
            LB = 0.1; UB = 0.9;
        else
            LB = varargin{1}; UB = varargin{2};
        end
        tmp = paretotails(x, LB, UB, 'kernel');
        Fxp = icdf(tmp, x_p);
        Fx  = icdf(tmp, x);
        X   = x;
    end   
       

    Fx = NaN;
    X  = NaN;
end
    
    
