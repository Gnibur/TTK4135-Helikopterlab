clc;

init;


h = 0.25; %% Discretization time step

%% Define system parameters

A = [1 h 0 0; 0 1 -h*K_2 0; 0 0 1 h; 0 0 -h*K_1*K_pp (1-h*K_1*K_pd)];

B = [0; 0; 0; h*K_1*K_pp];


%% Define objective function

lambda_0 = pi;
lambda_f = 0;

q = 1; % penalty
N = 100; % number of time steps / time horizont
t_sample = 0.25; % sampling time in seconds

H_i = blkdiag(1,0,0,0);

H = H_i;

for i = 1:(N-1)
    H = blkdiag(H,H_i);
end

for i = 1:N
   H = blkdiag(H, q); 
end
   

f = [-2*lambda_f 0 0 0]';
f = repmat(f,N,1);
f = [f' zeros(1,N)];

%% Define constraints

x0 = [lambda_0 0 0 0]';
xf = [lambda_f 0 0 0]';





Venstre = -A;
for i = 1:(N-2)
    Venstre = blkdiag(Venstre, -A);
end
Venstre = [zeros(4,4*N-4); Venstre]; % sleng på øverste linje
Venstre = [Venstre zeros(4*N, 4)]; % sleng på høyre kolonne
Venstre = Venstre + eye(4*N);



Hoyre = -B;
for i = 1:(N-1)
    Hoyre = blkdiag(Hoyre, -B);
end

A_eq = [ Venstre Hoyre];
b_eq = [A*x0; zeros(4*N-4,1)];

p_k_constraint = 30*pi/180;

x_k_constraint = [ Inf Inf p_k_constraint Inf]';

lb_vector = [repmat(-x_k_constraint,N,1); repmat(-p_k_constraint,N, 1)];
ub_vector = [repmat(x_k_constraint,N,1); repmat(p_k_constraint,N, 1)];


startingpoint = [x0; zeros(5*N-4,1)];

x = quadprog(H, f, [], [], A_eq, b_eq, lb_vector, ub_vector) %,startingpoint


u = x(4*N+1:end);


u = [zeros(20, 1); u; zeros(20,1)]; %% add 10 seconds before and after control input
 
t = 0.25:0.25:35;
t = t';

disp('plot')
figure(22)
stairs(t,u),grid













