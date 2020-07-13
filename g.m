function xplus = g(x)
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: g_ex1_2.m
%--------------------------------------------------------------------------
% Project: Simulation of a hybrid system (bouncing ball)
% Description: Jump map
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00


global n XSt XEn LamInd ZInd QInd MInd PInd KInd AlphaSt AlphaEn VSt VEn DeltaInd D_jSt D_jEn Delta_jEn Delta_jSt gamma theta finalvalfx Posvector lambda_s lambda_t mu PhiInd Alpha_barInd PhiMin Delta_j_temp


%state

X = x(XSt:XEn);
Lambda = x(LamInd);
Z = x(ZInd);
P = x(PInd);
Q = x(QInd);
M = x(MInd);
K = x(KInd);
Alpha = x(AlphaSt:AlphaEn);
V = x(VSt:VEn);
Delta = x(DeltaInd);
D_j = x(D_jSt:D_jEn);
Delta_j = x(Delta_jSt:Delta_jEn);
Phi = x(PhiInd);
Alpha_bar = x(Alpha_barInd);


% timer and position discrete dynamics

Timplus = 0;
Xplus = X;

% 1) Continue positive line search

if ((fx(X,Delta,Z, Delta_j)<=Z-rho(Delta)) && (P==1) && ((Q==1) || (Q==0)) && (M==0))

    Zplus = fx(X,Delta,Z, Delta_j);
    Qplus = 1;
    Operation_in_g=1
    if (K==0)
        
        Lambdaplus = Lambda+Delta_j(n)*P;
        
        if (gamma*Delta_j(n) <= lambda_t*Phi && gamma*Delta_j(n)>= lambda_s*Phi)
            Deltaplus = gamma*Delta_j(n);
            Delta_jplus = [Delta_j(1:n-1); gamma*Delta_j(n)];
        elseif (gamma*Delta_j(n) >= lambda_t*Phi)
            Deltaplus = lambda_t*Phi;
            Delta_jplus = [Delta_j(1:n-1); lambda_t*Phi];
        else
            Deltaplus = lambda_s*Phi;
            Delta_jplus = [Delta_j(1:n-1); lambda_s*Phi];
        end
    else
        
        Lambdaplus = Lambda+Delta_j(K)*P;

        if (gamma*Delta_j(K) <= lambda_t*Phi && gamma*Delta_j(K) >= lambda_s*Phi)
            Deltaplus = gamma*Delta_j(K);
            Delta_jplus = [Delta_j(1:K-1); gamma*Delta_j(K); Delta_j(K+1:n)];
        elseif (gamma*Delta_j(K) >= lambda_t*Phi)
            Deltaplus = lambda_t*Phi;
            Delta_jplus = [Delta_j(1:K-1); lambda_t*Phi; Delta_j(K+1:n)];
        else
            Deltaplus = lambda_s*Phi;
            Delta_jplus = [Delta_j(1:K-1); lambda_s*Phi; Delta_j(K+1:n)];
        end
    end
    Pplus = P; Mplus = M; Kplus = K; Alphaplus = Alpha; Vplus = V;
    D_jplus = D_j; Phiplus = Phi; Alpha_barplus = Alpha_bar;
    
% 2) Correct overshoot

elseif ((fx(X,Delta,Z, Delta_j) >=Z-rho(Delta)) && ((Q==1) || (Q==0)) && (M==0))

    Operation_in_g=2
    Pplus = -P;
    Qplus = Q+1;
    Mplus = 1;
        
    Lambdaplus = Lambda; Zplus = Z; Kplus = K; 
    Alphaplus = Alpha; Vplus = V; Deltaplus = Delta; D_jplus = D_j;
    Delta_jplus = Delta_j; Phiplus = Phi; Alpha_barplus = Alpha_bar;

% 3) Starting negative line search
elseif ((P==-1) && (Q==1) && (M==1))
    Operation_in_g=3
    Zplus = fx(X,Delta,Z, Delta_j);
    Mplus = 0;
    Lambdaplus = 0;

    Pplus = P; Qplus = Q; Kplus = K; 
    Alphaplus = Alpha; Vplus = V; Deltaplus = Delta; D_jplus = D_j;
    Delta_jplus = Delta_j; Phiplus = Phi; Alpha_barplus = Alpha_bar;
   
% 4) Continue a negative line search    
elseif ((fx(X,Delta,Z, Delta_j)<=Z-rho(Delta)) && (P==-1) && (Q==1) && (M==0))

    Operation_in_g=4
    Zplus = fx(X,Delta,Z, Delta_j);
    % option 1: put all Delta_j into plus and then only change one
    % corresponding
    % option 2: use an array like below
    if (K==0)
        Lambdaplus = Lambda+Delta_j(n)*P;
                        
        if (gamma*Delta_j(n) <= lambda_t*Phi && gamma*Delta_j(n)>= lambda_s*Phi)
            Deltaplus = gamma*Delta_j(n);
            Delta_jplus = [Delta_j(1:n-1); gamma*Delta_j(n)];
        elseif (gamma*Delta_j(n) >= lambda_t*Phi)
            Deltaplus = lambda_t*Phi;
            Delta_jplus = [Delta_j(1:n-1); lambda_t*Phi];
        else
            Deltaplus = lambda_s*Phi;
            Delta_jplus = [Delta_j(1:n-1); lambda_s*Phi];
        end
        
    else
        Lambdaplus = Lambda+Delta_j(K)*P;
        
        
        if (gamma*Delta_j(K) <= lambda_t*Phi && gamma*Delta_j(K) >= lambda_s*Phi)
            Deltaplus = gamma*Delta_j(K);
            Delta_jplus = [Delta_j(1:K-1); gamma*Delta_j(K); Delta_j(K+1:n)];
        elseif (gamma*Delta_j(K) >= lambda_t*Phi)
            Deltaplus = lambda_t*Phi;
            Delta_jplus = [Delta_j(1:K-1); lambda_t*Phi; Delta_j(K+1:n)];
        else
            Deltaplus = lambda_s*Phi;
            Delta_jplus = [Delta_j(1:K-1); lambda_s*Phi; Delta_j(K+1:n)];
        end
        
    end
    
    Pplus = P; Qplus = Q; Mplus = M; Kplus = K; Alphaplus = Alpha; Vplus = V;
    D_jplus = D_j; Phiplus = Phi; Alpha_barplus = Alpha_bar;

% 5) Update direction and start positive line search
elseif ((Q == 2))
    
    Operation_in_g=5
    Qplus = 0;
    Pplus = 1;
    Lambdaplus = 0;
    Mplus = 0;
    Kplus = mod(K+1,n+1);
    if (K==0)
        
        Alphaplus = Alpha + Lambda*V; %zeros(n,1)
        Alpha_barplus = Alpha_bar + norm(Lambda*V);
        Vplus = D_j(1:n);
        D_jplus = D_j;
        Deltaplus = Delta_j(mod(K+1,n+1));
        
        if (abs(Lambda)<=Delta_j(n)/2)
            
            if (theta*Delta_j(n)<= lambda_s*Phi)
                
                if (theta*Delta_j(n)<= lambda_s*PhiMin)
                    
                    Delta_j_temp(n) = lambda_s*PhiMin;
                else
                    
                    Delta_j_temp(n) = lambda_s*Phi;
                end
            else
                Delta_j_temp(n) = theta*Delta_j(n);
            end
        else
            Delta_j_temp(n) = Delta_j(n);
        end
        Delta_jplus = [Delta_j(1:n-1); Delta_j_temp(n)];
        Phiplus = Phi;
        
    elseif (K == n)
        
        Alphaplus = zeros(n,1)
        Alpha_barplus = 0
        Vplus = phi(D_j, Alpha, Lambda, V)
        D_jplus = [D_j(n+1:n*n); phi(D_j, Alpha, Lambda, V)]

        if (Alpha_bar + norm(Lambda*V)<=min(Delta_j)/2)
                       
            if (mu*Phi>= PhiMin)
                
                Phiplus = mu*Phi;
                for i=(1:length(Delta_j)-1)
                    if (Delta_j(i) >= mu*lambda_t*Phi)
                        Delta_j_temp(i) = mu*lambda_t*Phi;
                    elseif (Delta_j(i) <= mu*lambda_s*Phi)
                        Delta_j_temp(i) = mu*lambda_s*Phi; 
                    end
                end                
            else
                
                Phiplus = PhiMin
                for i=1:(length(Delta_j)-1)
                    if (Delta_j(i) >= lambda_t*PhiMin)
                        
                        Delta_j_temp(i) = lambda_t*PhiMin
                    elseif (Delta_j(i) <= lambda_s*PhiMin)
                        
                        Delta_j_temp(i) = lambda_s*PhiMin
                    end
                end
            end
            
            if (theta*Delta_j(K)<= mu*lambda_s*Phi)
                if (theta*Delta_j(K)<= lambda_s*PhiMin)
                    
                    Delta_j_temp(K) = lambda_s*PhiMin
                else
                    
                    Delta_j_temp(K) = lambda_s*Phi
                end
            elseif (theta*Delta_j_temp(K) >= mu*lambda_t*Phi)
                if (theta*Delta_j(K)<= lambda_s*PhiMin)
                    
                    Delta_j_temp(K) = lambda_s*PhiMin
                else
                    
                    Delta_j_temp(K) = mu*lambda_t*Phi
                end
            else
                if (theta*Delta_j(K)<= lambda_s*PhiMin)
                    
                    Delta_j_temp(K) = lambda_s*PhiMin
                else
                    
                    Delta_j_temp(K) = theta*Delta_j(K)
                end
            end
        else
            
            Phiplus = Phi;
            
            if (abs(Lambda)<=Delta_j(K)/2)
                if (theta*Delta_j(K)<= lambda_s*Phi)
                    
                    if (theta*Delta_j(K)<= lambda_s*PhiMin)
                        
                        Delta_j_temp(K) = lambda_s*PhiMin;
                    else
                        
                        Delta_j_temp(K) = lambda_s*Phi;
                    end
                else
                    
                    Delta_j_temp(K) = theta*Delta_j(K);
                end
            else
             
                Delta_j(K)
                Delta_j_temp(K) = Delta_j(K)
            end
            
            Delta_j_temp=[Delta_j(1:n-1);Delta_j_temp(K)]
        end
        
        Delta_jplus = [Delta_j_temp(2:n);Delta_j_temp(1)]               
        
        Deltaplus = Delta_j(1);
    else

        Alphaplus = Alpha + Lambda*V;
        Alpha_barplus = Alpha_bar + norm(Lambda*V);
        Vplus = D_j(K*n+1:(K+1)*n);
        D_jplus = D_j;
        Deltaplus = Delta_j(mod(K+1,n+1));
        
        if (abs(Lambda)<=Delta_j(K)/2)
            if (theta*Delta_j(K)<= lambda_s*Phi)

                if (theta*Delta_j(K)<= lambda_s*PhiMin)
                    Delta_j_temp(K) = lambda_s*PhiMin;
                else
                    Delta_j_temp(K) = lambda_s*Phi;
                end
            else

                Delta_j_temp(K) = theta*Delta_j(K);
            end
        else

             Delta_j_temp(K) = Delta_j(K);
        end
                    
        Delta_jplus = [Delta_j(1:K-1); Delta_j_temp(K) ;Delta_j(K+1:n)];
        Phiplus = Phi;
    end
    
    Zplus = fx(X,Delta,Z, Delta_j); 
    
% 6) All the other possible conditions:
else
    Operation_in_g=6
    
    Lambdaplus = Lambda; Zplus = Z; Pplus = P; Qplus = 2;
    Mplus = M; Kplus = K; Alphaplus = Alpha; Vplus = V;
    Deltaplus = Delta; D_jplus = D_j; Delta_jplus = Delta_j; Phiplus = Phi;
    Alpha_barplus = Alpha_bar;
end
   
% Print for debugging purposes 

finalvalfx=fx(X,Delta,Z, Delta_j)
X
Posvector = [Posvector, X];
P;
Q;
M;
Z
K
rhoD=rho(Delta)
Alpha
Lambda*V
V
Lambda
Alphaplus
Delta_jplus
Delta_j
D_j
xplus = [Xplus; Timplus; Lambdaplus; Zplus; Pplus; Qplus;...
    Mplus; Kplus; Alphaplus; Vplus; Deltaplus; D_jplus; Delta_jplus; Phiplus; Alpha_barplus];

end