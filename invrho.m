function DeltMin = invrho(bar_ns)

% As robustness is achieved if the minimum step size PhiMin is chosen
% bigger than the value satisfying the equation bar_ns = rho(PhiMinBar)/2, the
% function invrho computes an approximation (from above) of the value
% PhiMinBar, hence guaranteeing robustness to any noise whose norm is
% bounded by bar_ns.

global lambda_s

temp = 0;
while(rho(lambda_s*temp)<=2*bar_ns+0.0001)
   temp = temp + 0.01;
end

DeltMin = temp;
end