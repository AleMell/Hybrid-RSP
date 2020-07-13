function inside = D(x) 

% The jump set D

global TimInd

tau=x(TimInd);

if (tau>=1)
    inside = 1;
else
    inside = 0;
end



end