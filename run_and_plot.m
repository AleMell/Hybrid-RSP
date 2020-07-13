%--------------------------------------------------------------------------
clc

global  n


% initial conditions

x0 = [X0; Tim0; Lambda0; Z0; P0; Q0; M0; K0; Alpha0; V0; Delta0; D_j0; Delta_j0; Phi0; Alpha_bar0];


% simulation horizon
TSPAN=[0 100];
JSPAN = [0 100];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);

% simulate
[t,j,x] = HyEQsolver( @f,@g,@C,@D,x0,TSPAN,JSPAN,rule,options,'ode23t');

%%
%plot solution
if (n==3)
    figure(1) % position
    clf
    subplot(4,1,1), plotHarc(t,j,x(:,1));

    grid on
    Ylab1=ylabel('$x_1$');
    set(Ylab1,'Interpreter','latex');
    set(Ylab1,'FontSize',15);
    subplot(4,1,2), plotHarc(t,j,x(:,2));

    grid on
    Ylab2=ylabel('$x_2$');
    set(Ylab2,'Interpreter','latex');
    set(Ylab2,'FontSize',15);
    
    subplot(4,1,3), plotHarc(t,j,x(:,3));
    
    grid on
    Ylab3=ylabel('$x_3$');
    set(Ylab3,'Interpreter','latex');
    set(Ylab3,'FontSize',15);

    %subplot(3,1,3), plotHarc(t,j,(1-x(:,1)).^2+10.*(x(:,2)-x(:,1).^2).^2); % rosenbrock 
    grid on
    Ylab4=ylabel('$f(x)$');
    set(Ylab4,'Interpreter','latex');
    set(Ylab4,'FontSize',15);
    Xlab1=xlabel('$t$');
    set(Xlab1,'Interpreter','latex');
    set(Xlab1,'FontSize',15);

    
    
elseif (n==2)
    figure(1) % position
    clf
    subplot(3,1,1), plotHarc(t,j,x(:,1));

    grid on
    Ylab1=ylabel('$x_1$');
    set(Ylab1,'Interpreter','latex');
    set(Ylab1,'FontSize',15);

    subplot(3,1,2), plotHarc(t,j,x(:,2));
    
    grid on
    Ylab2=ylabel('$x_2$');
    set(Ylab2,'Interpreter','latex');
    set(Ylab2,'FontSize',15);
    
    % QUADRATIC
    
    subplot(3,1,3), plotHarc(t,j,1.*x(:,1).^2+5.*x(:,2).^2); 

    % DROP WAVE
    % %- (1+cos(12*(x(1)^2+x(2)^2)^(1/2)))/(0.5*(x(1)^2+x(2)^2)+2)
    %subplot(3,1,3), plotHarc(t,j, - (1+cos(12.*(x(:,1).^2+x(:,2).^2).^(1/2)))./(0.5.*(x(:,1).^2+x(:,2).^2)+2));
    
    % ROSENBROCK
    
    %subplot(3,1,3), plotHarc(t,j,(1-x(:,1)).^2+10.*(x(:,2)-x(:,1).^2).^2);
    
    grid on
    Ylab4=ylabel('$f(x)$');
    set(Ylab4,'Interpreter','latex');
    set(Ylab4,'FontSize',15);
    Xlab1=xlabel('$t$');
    set(Xlab1,'Interpreter','latex');
    set(Xlab1,'FontSize',15);
end

%%  PLOT OF PHASE PLANE
figure(2) % position
clf

lab1=xlabel('$x_1$');
set(lab1,'Interpreter','latex');
lab2=ylabel('$x_2$');
set(lab2,'Interpreter','latex');
lab3=zlabel('$f(x)$');
set(lab3,'Interpreter','latex');

grid on
hold on

points = 500;

X_ax1 = linspace(-3,3,points);
Y_ax1 = linspace(-3,3,points); 

[X_ax,Y_ax]=meshgrid(X_ax1,Y_ax1);

%%

%Z_ax = (1-X_ax).^2+10.*(Y_ax-X_ax.^2).^2;  % rosenbrock

%Z_ax=X_ax.^3-3*X_ax.*Y_ax.^2; %monkey saddle

%Z_ax=X_ax.^2-4.*Y_ax.^2; % classic saddle

Z_ax=X_ax.^2+5.*Y_ax.^2; % quadratic in fx

%Z_ax = -(1+cos(12.*(X_ax.^2+Y_ax.^2).^(1/2)))./(0.5.*(X_ax.^2+Y_ax.^2)+2);  % drop wave

%Z_ax=1/2*((X_ax.^4-16.*X_ax.^2+5.*X_ax)+(Y_ax.^4-16.*Y_ax.^2+5.*Y_ax)); %STYBLINSKI-TANG FUNCTION 

%Z_ax=(cos(X_ax) .* sin(Y_ax)) - (X_ax ./ ((Y_ax .^ 2) + 1));  %ADJ

contour3(X_ax,Y_ax,Z_ax,300);

%%
hold on

Minimum_for_quadratic = plot3(0,0,0,'r*','MarkerSize',5);

Initial_position = plot3(x(1,1), x(1,2), x(1,1).^2+5.*x(1,2).^2,'g*','MarkerSize',10);

% QUADRATIC

plotHarcColor3D(t,j,[x(:,1) x(:,2), x(:,1).^2+5.*x(:,2).^2],[],[],[],[]);

% DROP WAVE

%plotHarcColor3D(t,j,[x(:,1) x(:,2), - (1+cos(12.*(x(:,1).^2+x(:,2).^2).^(1/2)))./(0.5.*(x(:,1).^2+x(:,2).^2)+2)],[],[],[],[]);

% ROSENBROCK

%plotHarcColor3D(t,j,[x(:,1) x(:,2), (1-x(:,1)).^2+10.*(x(:,2)-x(:,1).^2).^2],[],[],[],[]); 

