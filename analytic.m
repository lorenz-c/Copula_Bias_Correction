function dC = condcopula(family, u, v, theta)

%	    P:  	Conditional cumulative  distribution of v given u!


dC = zeros(size(u)); 
nz = ((u ~= 0) & (v ~= 0));
u  = u(nz);
v  = v(nz);



switch lower(family)
    case 'clayton'
        dC(nz) = 1./(u.^(theta + 1).*(1./u.^theta + 1./v.^theta - 1).^(1./theta + 1));
     case 'frank'
        et  = exp(-theta);
        etu = exp(-theta.*u);
        etv = exp(-theta.*v);
        dC(nz)  = (etu.*(etv - 1))./(et + etu.*etv - etu - etv);
    case 'gumbel'
        dC(nz) = (exp(-((-log(u)).^theta + (-log(v)).^theta).^(1./theta)).*(-log(u)).^(theta - 1).*((-log(u)).^theta + (-log(v)).^theta).^(1./theta - 1))./u;
	case 'arch12'
        dC(nz) = ((1./u - 1).^(theta - 1).*((1./u - 1).^theta + (1./v - 1).^theta).^(1./theta - 1))./(u.^2.*(((1./u - 1).^theta + (1./v - 1).^theta).^(1../theta) + 1).^2);
    case 'arch14'
        dC(nz) = ((1./u.^(1./theta) - 1).^(theta - 1).*((1./u.^(1./theta) - 1).^theta + (1./v.^(1/theta) - 1).^theta).^(1./theta - 1))./(u.^(1./theta + 1).*(((1./u.^(1./theta) - 1).^theta + (1./v.^(1./theta) - 1).^theta).^(1./theta) + 1).^(theta + 1));
	case 'gb'
        dC(nz) = -v.*exp(-theta.*log(u).*log(v)).*(theta.*log(v) - 1);
	case 'fgm'
        dC(nz) = v.*(theta.*(u - 1).*(v - 1) + 1) + theta.*u.*v.*(v - 1);
    case 'amh'
        dC(nz) = (theta.*u.*v.*(v - 1))./(theta.*(u - 1).*(v - 1) - 1).^2 - v./(theta.*(u - 1).*(v - 1) - 1);
	case 'joe'
        dC(nz) = (((1 - v).^theta - 1).*(1 - u).^theta)./((u - 1).*((1 - u).^theta + (1 - v).^theta - (1 - u).^theta.*(1 - v).^theta).^((theta - 1)./theta));               
end




