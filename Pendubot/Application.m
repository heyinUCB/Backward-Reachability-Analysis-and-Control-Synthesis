clear
format long

%% Finite time horizon
T = 4;

%% System dynamics
x = mpvar('x',[4,1]);
x1 = x(1); 
x2 = x(2); 
x3 = x(3); 
x4 = x(4); 

%% create time t
pvar t

%% Define the system
f2 = - 10.6560*x1^3 + 11.5309*x1^2*x3 + 7.8850*x1*x3^2 + 0.7972*x2^2*x3 ...
  + 0.8408*x2*x3*x4 + 21.0492*x3^3 + 0.4204*x3*x4^2 + 66.5225*x1 - 24.5110*x3;
f4 = 10.9955*x1^3 - 48.9151*x1^2*x3 - 6.4044*x1*x3^2 - 2.3955*x2^2*x3 ...
  - 1.5943*x2*x3*x4 - 51.9088*x3^3 - 0.7971*x3*x4^2 - 68.6419*x1 + 103.9783*x3;
f = [x2; f2; x4; f4];
g2 = -10.0959*x3^2 + 44.2521;
g4 =  37.8015*x3^2 - 83.9120;
g = [0; g2; 0; g4];

dynamics.f = f;
dynamics.g = g;

%% target function
rT = x'*blkdiag(1/0.1^2, 1/0.35^2, 1/0.1^2, 1/0.35^2)*x - 1;

%% Initialize
timer = tic;

% V0val = 182.2102*x1^2 + 67.8981*x1*x2 + 314.3265*x1*x3 + 37.2705*x1*x4 + ...
%     6.4123*x2^2 + 59.0528*x2*x3 + 7.0314*x2*x4 + 138.9343*x3^2 + ...
%     32.5944*x3*x4 + 1.9521*x4^2;

load('u1T4_step4.mat','Vval')
V0val = Vval;

%%
max_iter = 10;
uM = 1;
um = -1;
[gamma_list, Vval, uval] = ...
    CLF(dynamics, T, x, t, rT, V0val, max_iter, uM, um);
compute_time = toc(timer)

%% Save data
filename = ['Reach',datestr(now, 'dd-mmm-yyyy-HH-MM-PM'),'.mat'];
save(filename)