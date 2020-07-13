function [value] = C(x) 

% The flow set C

global TimInd

tau = x(TimInd);

if (tau >= 0) && (tau < 1)  % flow condition
    value = 1;  % report flow
else 
    value = 0;  % do not report flow
end
end