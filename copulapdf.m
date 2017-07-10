function c = copulapdf(family, varargin)
%
%   FUNCTION C = COPULAPDF(FAMILY, U, ALPHA)
%
%   Return c(u,v|ALPHA), the copula density. 
%
%   INPUTS
%       FAMILY: one of {'ind', 'gaussian', 'gumbel' 'clayton' 'frank' 'amh' 'joe' 'fgm' 'arch12' 'arch14'}       
%       U: Nx2 vector (u,v) in [0,1]^2.
%       ALPHA: 1xM vector of copula parameters
%           or
%   COPULAPDF(FAMILY, U1, U2, ALPHA)
%       U1, U2:  Matrices or row vectors.
%                If U1 is (1xN) and U2 is (1xM), c is (NxM)
%                If U1 is (NxM), U2 must be (N,M).
%       ALPHA:   Scalar copula parameter
%
%   OUTPUT
%       C:       Copula density c(u,v|ALPHA) (NxM).

%   Guillaume EVIN, 13 May, 2004.
%   D. Huard, Nov. 2006
%  
if nargin == 3
    U = varargin{1};
    u = U(:,1);
    v = U(:,2);
    alpha = varargin{2};
elseif nargin == 4
    u = varargin{1};
    v = varargin{2};
    alpha = varargin{3};
end

% Check alpha is in the domain covered by the family.
pass = check_alpha(family, alpha);
if ~all(pass)
    error('Some parameters are not valid.\n%f', alpha(~pass))
end

% Check u,v are in [0,1]^2
if any( (u < 0) | (u > 1)) | any((v < 0) | (v > 1) )
    error('Some quantiles are outside the unit hypercube.')
end

if nargin == 3
    % Shape checking
    [NU, MU] = size(U);
    [NA, MA] = size(alpha);
    
    if MU ~= 2
        error('Bad shape. U is not Nx2, but rather %s.', mat2str(size(U)))
    end
    
    % Reshape ALPHA
    if NA == 1 
        alpha = repmat(alpha, NU, 1);
    elseif NA ~= NU && NU ~= 1
        error('Number of parameters must be 1, identical to number of couples in U, or a row vector.')
    end
    
    % Reshape u,v
    if NU == 1
        u = repmat(u, NA, MA);
        v = repmat(v, NA, MA);
    else
        u = repmat(u, 1, MA);
        v = repmat(v, 1, MA);
    end
elseif nargin == 4
    if ~all(size(alpha)==1)
        error('Alpha must be a scalar.')
    end
    su = size(u);
    sv = size(v);
    if ~any(su==1)
        if all(su ==sv)
            alpha = repmat(alpha, size(u));
        else
            error('If U1 and U2 are matrices, they must have the same size.')
        end
    else
        [u,v] = meshgrid(u,v);
        alpha = repmat(alpha, size(u));
    end
end

switch lower(family)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ellipitical copulas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 'gaussian'
        v1 = norminv(u);
        v2 = norminv(v);
        % Octave ----
        %v1 = normal_inv(u);
        %v2 = normal_inv(v);
        % -----------
        c = (1./sqrt(1-alpha.^2)).*exp(-(v1.^2+v2.^2-(2.*alpha).*v1.*v2)./(2*(1-alpha.^2)) + (v1.^2+v2.^2)./2);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Archimedean copulas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'ind'
        c = ones(size(u)); % the density for the independant copula is one
        
    case 'gumbel'
        % the gumbel copula : C(u,v) = exp(-((-ln u)^alpha + (-ln v)^alpha)^(1/alpha))
        C = copulacdf('gumbel',[u(:,1),v(:,1)], alpha(1,:));
        logC = (-log(C)).^alpha;
        v1 = ((-log(u)).^(alpha-1))./u ;
        v2 = ((-log(v)).^(alpha-1))./v ;
        c = (v1.*v2.*C).*(logC.^((2-2.*alpha)./alpha)-(1-alpha).*(logC.^((1-2.*alpha)./alpha)));
   
    case 'clayton'
        % (u^(-a)+v^(-a)-1)^(-(1+2*a)/a)*v^(-a-1)*u^(-a-1)+(u^(-a)+v^(-a)-1)^(-(1+2*a)/a)*u^(-a-1)*v^(-a-1)*a 
        %C1 = copulacdf('clayton',[u(:,1),v(:,1)], alpha(1,:));
        %c= (u.^(-alpha-1)).*(v.^(-alpha-1)).*(alpha+1).*C1.^(1+2.*alpha);
        t1 = u .^ (-alpha);
        t2 = v .^ (-alpha);
        t8 = (t1 + t2 - 1) .^ (-(1 + 2 * alpha) ./ alpha);
        t9 = -alpha - 1;
        t10 = v .^ t9;
        t12 = u .^ t9;
        t13 = t8 .* t10 .* t12;
         c = t13 + t13 .* alpha;
        
    case 'frank'
        % the frank copula : C(u,v) = (-1/alpha)*(1 + (exp(-alpha*u)-1)(exp(-alpha*v)-1)/(exp(-alpha)-1))
        t3 = exp(-alpha .* (u + v));
        t5 = exp(-alpha);
        t6 = alpha .* u;
        t7 = alpha .* v;
        t9 = exp(-t6 - t7);
        t11 = exp(-t6);
        t12 = exp(-t7);
        t14 = (t5 + t3 - t11 - t12) .^ 2;
        c   = -t3 .* alpha .* (-1 + t5 + t3 - t9) ./ t14;
           
    case 'frank_genest'
        % The frank copula from Genest (1987))
        % The alphaameterization is different
        v1 = alpha.^u ;
        v2 = alpha.^v ;
        v3 = alpha.^(u +v);
        c = ((alpha-1).*log(alpha).*v3)./(((alpha-1)+(v1-1).*(v2-1)).^2);
        
    case 'joe'
        % the joe copula : C(u,v) = 1 - [(1-u)^alpha+(1-v)^alpha-
        % (1-u)^alpha(1-v)^alpha]^(1/alpha)
        C = 1 - copulacdf('joe',[u(:,1),v(:,1)], alpha(1,:));
        v1=1-u;
        v2=1-v;
        c = (v1.^(alpha-1)).*(v2.^(alpha-1)).*alpha.*(C.^(1-alpha)) + (alpha-1).*((v1.^(alpha-1)).*(v2.^(alpha-1))...
            .*(1-v1.^alpha).*(1-v2.^alpha)).*(C.^(1-2*alpha));
        
    case 'arch12'
        % the 12th archimedean copula in Nelsen.
        t1 = -1 + v  ;
        t2 = 1 ./ v  ;
        t3 = t1 .* t2 ;
        t4 = (-t3) .^ alpha ;
        t5 = -1 + u  ;
        t6 = 1 ./ u  ;
        t7 = t5 .* t6 ;
        t8 = (-t7) .^ alpha ;
        t9 = t8 + t4 ;
        t11 = 1 ./ alpha ;
        t14 = t9 .^ (-2 .* (-1 + alpha) .* t11) ;
        t15 = t4 .* t14 ;
        t16 = 3 .* alpha ;
        t17 = (-t7) .^ t16 ;
        t19 = 2 .* alpha ;
        t20 = (-t7) .^ t19 ;
        t22 = (-t3) .^ t19 ;
        t25 = t8 .* t14 ;
        t26 = (-t3) .^ t16 ;
        t28 = t4 .* t8 ;
        t29 = t9 .^ t11 ;
        t47 = 1 + t29 ;
        t48 = t47 .^ 2 ;
        c = (t15 .* t17 + 2 .* t14 .* t20 .* t22 + t25 .* t26 - t28 .* t29 + t28 .* ...
            alpha .* t29 + t15 .* alpha .* t17 + 2 .* t22 .* t20 .* t14 .* alpha + t25 .* alpha .* ...
        t26) .* t6 ./ t5 .* t2 ./ t1 ./ t48 ./ t47 ./ (t20 + 2 .* t28 + t22) ;
%         theta = alpha;
% 
%         c = (2*(1./u - 1).^(theta - 1).*(1./v - 1).^(theta - 1).*((1./u - 1).^theta + (1./v - 1).^theta).^(2./theta - 2))./(u.^2.*v.^2.*(((1./u - 1).^theta + (1./v - 1).^theta).^(1./theta) + 1).^3) - (theta.*(1./theta - 1).*(1./u - 1).^(theta - 1).*(1./v - 1).^(theta - 1).*((1./u - 1).^theta + (1./v - 1).^theta).^(1./theta - 2))./(u.^2.*v.^2.*(((1./u - 1).^theta + (1./v - 1).^theta).^(1./theta) + 1).^2);
    case 'arch14'
        % the 14th archimedean copula in Nelsen.
%         v1 = (-1 + u.^(-1./alpha)).^alpha;
%         v2 = (-1 + v.^(-1./alpha)).^alpha;
%         c = v1.*v2.*((v1+v2).^(-2+1./alpha)).*(1+(v1+v2).^(1./alpha)).^(-2-alpha).*(-1+alpha+2.*alpha.*...
%         ((v1+v2).^(1./alpha)))./(alpha.*u.*v.*(-1+u.^(1./alpha)).*(-1+v.^(1./alpha)));
        t1 = 1 ./ alpha;
        t2 = u .^ (-t1);
        t3 = t2 - 1;
        t4 = t3 .^ alpha;
        t5 = v .^ (-t1);
        t6 = t5 - 1;
        t7 = t6 .^ alpha;
        t8 = t4 + t7;
        t9 = t8 .^ t1;
        t10 = 1 + t9;
        t11 = t10 .^ (-alpha);
        t12 = t9 .^ 2;
        t13 = t11 .* t12;
        t15 = 1 ./ v;
        t17 = 1 ./ t6;
        t17 = replaceInf(t17);
        t18 = t5 .* t15 .* t17;
        t20 = t8 .^ 2;
        t21 = replaceInf(1 ./ t20);
        t22 = t10 .^ 2;
        t24 = t21 ./ t22;
        t27 = t2 ./ u;
        t28 = 1 ./ t3;
        t28 = replaceInf(t28);
        t29 = t27 .* t28;
        t32 = t11 .* t9;
        t34 = t7 .* t5;
        t39 = 1 ./ t10;
        c = t13 .* t7 .* t18 .* t24 .* t4 .* t29 - t32 .* t1 .* t34 .* t15 .* t17 .* t21 .* t4 .* t27 .* t28 .* t39 + t32 .* t4 .* t29 .* t21 .* t39 .* t7 .* t18 + t13 .* t4 .* t29 .* t24 .* t1 .* t34 .* t15 .* t17;
%           c = ((1./u.^(1./alpha) - 1).^alpha.*(1./v.^(1./alpha) - 1).^alpha.*((1./u.^(1./alpha) - 1).^alpha + (1./v.^(1./alpha) - 1).^alpha).^(1./alpha - 2).*(alpha + 2*alpha.*((1./u.^(1./alpha) - 1).^alpha + (1./v.^(1./alpha) - 1).^alpha).^(1./alpha) - 1))./(alpha.*u.*v.*(u.^(1./alpha) - 1).*(v.^(1./alpha) - 1).*(((1./u.^(1./alpha) - 1).^alpha + (1./v.^(1./alpha) - 1).^alpha).^(1./alpha) + 1).^(alpha + 2));
    case 'fgm'
        % the Farlie-Gumbel-Morgenstern Copula
        % = ...
        t3 = alpha .* u ;
        c = 1 + alpha - (2 .* alpha) .* v  - 2 .* t3 + 4 .* t3 .* v ;
        
    case 'amh'
        % the Ali-Mikhail-Haq Copula
        % = (1-2*a+a^2+a*u*v+a*v+a*u-a^2*v-a^2*u+a^2*u*v)/(-1+a-a*v-a*u+a*u*v)^3
        t2 = alpha .* u ;
        t3 = t2 .* v ;
        t4 = alpha.^2;
        t5 = t4.*u ;
        t7 = alpha.*v ;
        t10 = -1 + alpha - t7 - t2 + t3;
        t11 = t10.^2;
        c = -(1 - 2 .* alpha + t3 + t5 .* v  + t2 + t7 - t5 - t4 .* v  + t4) ./ t11 ./ t10;
        
    otherwise
        error('Copula family ''%s'' not recognized.', family)
end

if any(any(isnan(c)))
    warning('NaN in %s at!', family)
end

function x = replace0(x)
% 
% Return x with zeros replaced by eps.

i = x == 0;
x(i) = eps;

function x = replaceInf(x)
%
% Return realmax instead of Inf.
%
i = isinf(x);
x(i) = sign(x(i)).*realmax;