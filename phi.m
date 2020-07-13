function d = phi(D_j, Alpha, Lambda, V)
% phi 

global n taudet

if(n>1)
    dM = abs(det([vec2mat(D_j(n+1:n*n), n)' ((Alpha+Lambda*V))]))
    if (1000*dM>=taudet)
        Alpha
        Lambda
        V
        Alpha+Lambda*V
        temp= (Alpha+Lambda*V)/(norm(Alpha+Lambda*V))

    else 
        Alpha
        Lambda*V
        Alpha+Lambda*V
        vec2mat(D_j(n+1:n*n), n)'
        temp = D_j(1:n);
    end
elseif (n==1)
    if (Alpha+Lambda*V~=0)
    temp= (Alpha+Lambda*V)/(norm(Alpha+Lambda*V));

    else 
        temp = D_j(1:n);
    end
else 
    error('Error: n<=0 is not supported');
end

d=temp;

end

