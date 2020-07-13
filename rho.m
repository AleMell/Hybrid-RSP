function smalloDelta = rho(Delta)

% The map rho used to impose the sufficient decrese condition. It is
% defined as a strictly increasing positive definite o(Delta) function for 
% Delta->0

if (Delta<=exp(1))
    smalloDelta = Delta^(1/Delta);
else
    smalloDelta = Delta+(exp(1/exp(1))-exp(1));

end

