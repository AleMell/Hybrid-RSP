%clc
%clear all
%close all

% INITIALIZATION

global n XSt XEn TimInd LamInd ZInd QInd MInd PInd KInd AlphaSt AlphaEn VSt VEn DeltaInd D_jSt D_jEn Delta_jEn Delta_jSt RegDim a gamma theta taudet lambda_s lambda_t mu PhiInd Alpha_barInd PhiMin bar_ns Mode noise PhiMin_theoretical finalvalfx

n=2;             % space dimension
rand_var = 0;    % logic variable defining if the initial condition have to be randomly defined (1) or not (0)
Delta_jDim = n;    % dimension of the controller state Delta (vector of step-sizes, one for each direction)
AlphaDim = n;    % dimension of the controller state alpha (total displacement per cycle)
D_jDim = n*n;    % dimension of the vector of directions 
VDim = n;        % dimension of the vector v (current direction getting explored)
% dimension of timer tau = 1
% dimension of lambda = 1       (current displacement along the current direction)  *TO DO(could be substituted with the corresponding alpha_j)
% dimension of z = 1            (previous measurement)
% dimension of p = 1            ( {1,-1} direction of travel along the current direction)
% dimension of q = 1            ( {0,1,2} states of the line minimization)
% dimension of m = 1            ( {0,1} first and second (in the opposite direction) search along a direction)
% dimension of k = 1            ( {0,1,..,n} state of a cycle)
% dimension of Delta = 1        ( current step-size) *TO DO (can be substituted with the actual Delta_j)
% dimension of Phi = 1          ( global step size )
% dimension of alpha_bar = 1    ( norm of the total distance travelled in cycle k )


a = 1;            % speed constant

bar_ns = 0.1;     % upper bound on the norm of the noise signal

gamma = 1;        % increment of Delta constant
theta = 1;        % reduction of Delta constant1
taudet = 0.001;     % minimum determinant
lambda_s = 1;   % constant defining the lower bound on Delta based on Phi, namely Delta>=lambda_s*Phi (same for all Delta_j)
lambda_t = 1;       % constant defining the upper bound on Delta based on Phi, namely Delta<=lambda_t*Phi (same for all Delta_j)
mu = 0.15;          % decrement of Phi constant. Has to be in (0, 1/lambda_t)

% indexes of the state in the whole state vector x
XSt = 1;                                          % position start 1
XEn = n;                                          % position end n
TimInd= XEn + 1;                                  % timer n+1
LamInd = TimInd + 1;                              % lambda n+2
ZInd = LamInd + 1;
PInd = ZInd + 1;
QInd = PInd + 1;
MInd = QInd + 1;
KInd = MInd + 1;
AlphaSt = KInd + 1;
AlphaEn = AlphaSt + (AlphaDim - 1);
VSt = AlphaEn + 1;
VEn = VSt + (VDim -1);
DeltaInd = VEn + 1;
D_jSt = DeltaInd + 1;
D_jEn = D_jSt + (D_jDim -1);
Delta_jSt = D_jEn +1;
Delta_jEn = Delta_jSt + (Delta_jDim -1);
PhiInd = Delta_jEn +1;
Alpha_barInd = PhiInd + 1;

RegDim= 10 + Delta_jDim + AlphaDim + D_jDim + VDim;    % number of states of the controller

% Choose mode of operation:

%Mode = 1; % for the noiseless case
Mode = 2; % for the noisy case without lower bound > 0 on the minimum step size. In this case, the noise is
          % tailored in order to make the state escape to the horizon
%Mode = 3; % for the noisy case with lower bound > 0 on the minimum step size. Knowledge of an upper bound on the
          % infinite norm of the noise signal is required.
          

% By default, in Mode 3, the lower bound on the minimum step size is
% defined by the minimum possible step size able to guarantee robustness,
% namely such that namely any PhiMin>PhiMinBar, with PhiMinBar such that 
% bar_ns = rho(lambda_s*PhiMinBar)/2

PhiMin_theoretical = invrho(bar_ns);

if (Mode == 3)
    PhiMin = PhiMin_theoretical;
else PhiMin = 0;
end

% noise signal initilization

noise = 0;

%% Definition of the matrix with column the directions of exploration
%Different examples:

% 1) Initial set of directions used in the examples on the paper

DirMat=[cos(pi/8) sin(pi/8); -sin(pi/8) cos(pi/8)];

% 2) The identity matrix

%DirMat=eye(n);

% 3) A random generated matrix

% Temp1 = (rand(n,1)-ones(n,1)/2);
% Temp = Temp1/norm(Temp1);
% i=2;
% while (i<=n)
%     Tempr = (rand(n,1)-ones(n,1)/2);
%     Tempn = Tempr/norm(Tempr)
%     if(rank([Temp,Tempn])==i)
%         Temp = [Temp,Tempn]
%         i=i+1;
%     end
% end
% DirMat = Temp;

%% Initial conditions
if(rand_var ==0)
    
    Alpha_bar0 = 0;
    X0 = [1.5;0];
    Tim0 = 0;
    Lambda0 = 0;
    %Z0 = fx(X0,Delta0,Z0);
    Z0=0;
    P0 = 1;
    Q0 = 0;
    M0 = 0;
    K0 = 0;
    Alpha0 = zeros(n,1);
    V0 = DirMat(:,n);
    %Delta0 = PhiMin;
    Delta0 = 0.01;
    D_j0 = DirMat(:);
    Delta_j0 = ones(n,1)*Delta0;
    %Phi0 = PhiMin;
    Phi0 = 0.01;
    
    if(Mode == 3 && Phi0 < PhiMin)
        Phi0 = PhiMin;
    end
    
elseif (rand_var == 1)
    
    minInCond = -10;
    maxInCond = 10;
    LambdaRange = 5;
    Alpha_barRange = 5;
    PRange = [-1;1];
    QRange = [0;1;2];
    MRange = [0;1];
    AlphaRange = 5;
    DeltaRange = 5;
    PhiRange = 5;
    
    X0= minInCond + (maxInCond-minInCond).*rand(n,1);
   
    Tim0 = rand; %it is not considering the case Tim0 = 0 or = 1
    
    Lambda0 = -LambdaRange + (2*LambdaRange).*rand;
    Alpha_bar0 = -Alpha_jRange + (2*Alpha_jRange).*rand;
    %Z0 = fx(X0,Delta0,Z0);
    Z0=0;
    
    P0 = PRange(randi(length(PRange)));
    Q0 = QRange(randi(length(QRange)));
    M0 = MRange(randi(length(MRange)));
    K0 = randi([0 n]);
    Alpha0 = -AlphaRange + (2*AlphaRange).*rand(n,1);
    
    %Vtemp = rand(n,1);
    %V0 = Vtemp/norm(Vtemp);
    V0 = DirMat(:,n);
    
    Phi0 = PhiRange*rand;
    
    if (Mode ==3 && Phi0 <= PhiMin)
        Phi0 = PhiMin;
    end
        
    Delta0 = DeltaRange*rand;
    
    D_j0 = DirMat(:);
    %Delta_j0 = DeltaRange.*rand(n,1);
    Delta_j0 = ones(n,1)*Delta0;

end