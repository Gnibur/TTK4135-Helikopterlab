% ********************************************************
% *                                                      *
% *       Optimization and Control                       *
% *                                                      *
% *       Helikopterlab                                  *
% *                                                      *
% * MAIN1NY.m                                            *
% *                                                      *
% *                                                      *
% * Updated on 04/2009 by Agus Ismail Hasan              *
% *                                                      *
% ********************************************************

clear; clc;

init;
delta_t	  = 0.25;	                    % sampling time
h = delta_t;
sek_forst = 5;

% System model. x=[lambda r p p_dot]'

A = [ 1 h  0           0              0            0; 
      0 1 -h*K_2       0              0            0;
      0 0  1           h              0            0;
      0 0 -h*K_1*K_pp (1-h*K_1*K_pd)  0            0;
      0 0  0           0              1            h;
      0 0  0           0              -h*K_3*K_ep -h*K_3*K_ed+1];
  
B1 = [0; 0; 0; h*K_1*K_pp];
  
B = [0           0; 
     0           0; 
     0           0; 
     h*K_1*K_pp  0; 
     0           0;
     0           h*K_3*K_ep];


% Number of states and inputs

mx = size(A,2);                        % Number of states (number of columns in A)
mu = size(B,2);                        % Number of inputs(number of columns in B)

% Initial values

x1_0 = pi;                              % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x5_0 = 0;                               % e
x6_0 = 0;                               % e_dot
x0   = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]';          % Initial values

% Time horizon and initialization

N  = 100;                                % Time horizon for states (40)
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;                                 % Initial value for optimization

% Bounds

p_lb  = -30*pi/180;                     % Lower bound on control -- u1, p_k
p_ub  = 30*pi/180;                      % Upper bound on control -- u1, p_k


lambda_t = 2*pi/3;
alpha = 1;
beta = 1;

xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = p_lb;                         % Lower bound on state x3
xu(3)   = p_ub;                         % Upper bound on state x3

xl(2) = -0.1;
xu(2) = 0.1;

% Generate constraints on measurements and inputs

vlb       = [repmat(xl,N,1);repmat([p_lb; -Inf],M,1)];
vub       = [repmat(xu,N,1); repmat([p_ub; Inf],M,1)];
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem) 

Q1 = zeros(mx,mx);
Q1(1,1) = 1;                             % Weight on state x1 (1)
Q1(2,2) = 0;                             % Weight on state x2
Q1(3,3) = 0;                             % Weight on state x3
Q1(4,4) = 0;                             % Weight on state x4
Q1(5,5) = 0;                             % Weight on state x4
Q1(6,6) = 0;                             % Weight on state x4


P1 = zeros(mu, mu);                         
P1(1,1) = 1;                             % Weight on input pc
P1(2,2) = 1;                             % Weight on input ec

Q = blkdiag(kron(eye(N), Q1), kron(eye(N), P1));
c = zeros(N*mx+M*mu,1);                  % Generate c

% Generate system matrices for linear model

Aeq = [kron(diag(ones(N-1,1), -1), -A) + eye(N * mx), kron(eye(N),-B)];
beq = [A*x0; zeros((N-1)*mx,1)]; 	     

objfun = @(x) x'*Q*x+c'*x;
%options = optimset('Algorithm', 'sqp');

% Solve Qp problem with non-linear model
tic
[z,lambda] = fmincon(objfun,zeros(N*(mx+mu),1),[],[],Aeq,beq,vlb,vub,@nonlinearconstraints);
t1=toc;

% Calculate objective value

phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+Q(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end

% Extract control inputs and states

u  = [z(N*mx+1:N*mx+M*mu)]; % Control input from solution
u_pc = [u(1:mu:M*mu); 0];
u_ec = [u(2:mu:M*mu); 0];

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(5:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);z(6:mx:N*mx)];              % State x6 from solution

Antall = 5/delta_t;
Nuller = zeros(Antall,1);
Enere  = ones(Antall,1);

u_pc = [Nuller; u_pc; Nuller];
u_ec = [Nuller; u_ec; Nuller];
x1   = [pi*Enere; x1; Nuller];
x2   = [Nuller; x2; Nuller];
x3   = [Nuller; x3; Nuller];
x4   = [Nuller; x4; Nuller];
x5   = [Nuller; x5; Nuller];
x6   = [Nuller; x6; Nuller];

%save trajektory

% figure
t = 0:delta_t:delta_t*(length(u_pc) - 1);                % real time
t = t';

figure(2)
subplot(411)
stairs(t,u_pc),grid
ylabel('u_pc')
subplot(412)
stairs(t,u_ec),grid
ylabel('u_ec')

subplot(413)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(414)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')

figure(3);
subplot(411)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(412)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')
subplot(413)
plot(t,x5,'m',t,x5','mo'),grid
xlabel('tid (s)'),ylabel('e')
subplot(414)
plot(t,x6,'m',t,x6','mo'),grid
xlabel('tid (s)'),ylabel('edot')
