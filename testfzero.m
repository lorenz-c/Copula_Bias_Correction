function alpha = testfzero(tau)

for i = 1:length(tau)
    if tau(i) == 0
        alpha(i) = 1;
    else
        %alpha(i) = %fzero(@invcopulastat,[1+eps, 54],[],'joe', tau(i));
        fun      = @(alpha) invcopulastat(alpha, 'joe', tau(i));
        keyboard
        alpha(i) = fzero(fun, [1+eps, 100]);
    end
    
end

function err = invcopulastat(alpha, family, target_tau)
% Return difference between target_tau and tau computed from copulastat
% with a guess for alpha.
guess_tau = copulastat(family, alpha);
err       = guess_tau - target_tau;