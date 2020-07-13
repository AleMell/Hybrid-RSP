function xdot = f(x)

% The flow map F definining the continuous time closed loop dynamics in the
% flow set C

global PInd VSt VEn DeltaInd RegDim a

P = x(PInd);
Delta = x(DeltaInd);
V = x(VSt:VEn);

xdot = [a*P*Delta*V; 1; zeros(RegDim-1,1)];

end