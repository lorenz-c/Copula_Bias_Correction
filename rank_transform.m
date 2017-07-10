function [Fx, F, X] = rank_transform(x, x_p, inverse)

if nargin < 3, inverse = 0; end
if nargin < 2, x_p = x; end

n = length(x);

for i = 1:n
    tmp  = find(x <= x(i));
    F(i) = length(tmp)/(n + 1);
    X(i) = x(i);
end
    

[F_unique, un_ids] = unique(F);
X_unique           = X(un_ids);

[F, ids] = sort(F_unique, 'ascend');
X        = X_unique(ids);

if inverse == 0
    Fx = interp1(X, F, x_p, 'linear', 'extrap');
elseif inverse == 1
    Fx = interp1(F, X, x_p, 'linear', 'extrap');
end


    
    

