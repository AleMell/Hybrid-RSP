function smalloDelta = rho(Delta)
%RHO A positive definite o(Delta) function for Delta->0

%smalloDelta = Delta^3;

if (Delta<=exp(1))
    smalloDelta = Delta^(1/Delta);
else
    smalloDelta = Delta+(exp(1/exp(1))-exp(1));

end

