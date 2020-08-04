clear
format long

%% Finite time horizon
T = 2;

%% System dynamics
x = mpvar('x',[6,1]);
x1 = x(1); 
x2 = x(2); 
x3 = x(3); 
x4 = x(4);
x5 = x(5);
x6 = x(6);

%% create time t
pvar t

%% Define the system
% parameters
g = 9.81;
K = 0.89/1.4;
d0 = 70;
d1 = 17;
n0 = 55;
f = [x3; x4; 0; -g; x6; -d0*x5-d1*x6];
g1 = [0; 0; K*(-0.166*x5^3+x5); K*(-0.498*x5^2+1); 0; 0];
g2 = [0; 0; 0; 0; 0; n0];

dynamics.f = f;
dynamics.g1 = g1;
dynamics.g2 = g2;

%% target function
rt.rt1 = (1.7 - x1)*(x1 + 1.7);
rt.rt2 = (0.85 - x2)*(x2 + 0.85);
rt.rt3 = (0.8 - x3)*(x3 + 0.8);
rt.rt4 = (1 - x4)*(x4 + 1);
rt.rt5 = (pi/12 - x5)*(x5 + pi/12);
rt.rt6 = (pi/2 - x6)*(x6 + pi/2);

%% Initialize
timer = tic;
P = Vinitialize;
V0val = x'*P*x;
% load('V2u2.mat','Vval')
% V0val = Vval;

%%
max_iter = 30;
uM = [1.5+g/K; pi/12];
um = [-1.5+g/K; -pi/12];
epsilon = 1e-4;
[gamma_list, Vval, u1val, u2val] = ...
    CLF(dynamics, T, x, t, rt, V0val, max_iter, uM, um, epsilon);
compute_time = toc(timer)

%% Save data
filename = ['Reach',datestr(now, 'dd-mmm-yyyy-HH-MM-PM'),'.mat'];

save(filename)