function [pi, varargout] = param_transform(x, D, inverse)

if nargin < 3, inverse = 0; end

DistName = D(1).DistName;
Params   = D(1).Params;
    
if inverse == 0
    
    if length(Params) == 1
        % One-parameter distributions
        pi = cdf(DistName, x, Params(1));
    elseif length(Params) == 2
        % Two-parameter distributions
        pi = cdf(DistName, x, Params(1), Params(2));
    elseif length(Params) == 3
        % Three-parameter distributions
        pi = cdf(DistName, x, Params(1), Params(2), Params(3));    
    end
    
elseif inverse == 1
    
    if length(Params) == 1
        % One-parameter distributions
        pi = icdf(DistName, x, Params(1));
    elseif length(Params) == 2
        % Two-parameter distributions
        pi = icdf(DistName, x, Params(1), Params(2));
    elseif length(Params) == 3
        % Three-parameter distributions
        pi = icdf(DistName, x, Params(1), Params(2), Params(3));    
    end
    
end

varargout{1} = DistName;
varargout{2} = Params;
    
    
