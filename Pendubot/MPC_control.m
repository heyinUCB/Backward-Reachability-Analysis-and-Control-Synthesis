% This script is aming at using control Lyapunov function to design
% controller
clear

%% System dynamics
% time horizon for MPC
N = 10;
% number of states
nX = 4;

% define decision variables x and u
x = sdpvar(nX, N+1);
u = sdpvar(1,N);

% symbolic x
xsym = sdpvar(nX, 1);
x1 = xsym(1);
x2 = xsym(2);
x3 = xsym(3);
x4 = xsym(4);

% symbolic t
t = sdpvar(1,1);

% Define the system
f2 = - 10.6560*x1^3 + 11.5309*x1^2*x3 + 7.8850*x1*x3^2 + 0.7972*x2^2*x3 ...
    + 0.8408*x2*x3*x4 + 21.0492*x3^3 + 0.4204*x3*x4^2 + 66.5225*x1 - 24.5110*x3;
f4 = 10.9955*x1^3 - 48.9151*x1^2*x3 - 6.4044*x1*x3^2 - 2.3955*x2^2*x3 ...
    - 1.5943*x2*x3*x4 - 51.9088*x3^3 - 0.7971*x3*x4^2 - 68.6419*x1 + 103.9783*x3;
f = [x2; f2; x4; f4];
g2 = -10.0959*x3^2 + 44.2521;
g4 =  37.8015*x3^2 - 83.9120;
g = [0; g2; 0; g4];

% discrete time: for simulation
dt = 0.02;

% sample time: for MPC controller
Ts = 0.02;

% time interval [0, T]
T = 2;

% initial condition
x0 = [-0.35; 2.6; 0.35; -4];

%% load storage function
file = 'u1T4_step5.mat';
Vcoeff = loadVcoeff(file);

Vdeg = 4;
Vmonom_yalmip = monolist([xsym; t], Vdeg);
Vval_yalmip = Vcoeff'*Vmonom_yalmip;
DVvalx = jacobian(Vval_yalmip, xsym);
DVvalt = jacobian(Vval_yalmip, t);

%% control constraints
umin = -1;
umax = 1;

%% Main algorithm
xsim = [x0];
usim = [];
time_line = 0:dt:T;
for i = 1:length(time_line)
    % initialize the cost and constraints
    constr = [];
    
    tk = time_line(i);
    for k = 1:N
        % plug in the measured state and time
        DVx = replace(replace(DVvalx, xsym, x(:,k)),t,tk);
        DVt = replace(replace(DVvalt, xsym, x(:,k)),t,tk);
        fk = replace(f,xsym,x(:,k));
        gk = replace(g,xsym,x(:,k));
        constr = [constr; ...
            DVt + DVx*(fk + gk*u(:,k)) <= 0;...
            x(:,k+1) == x(:,k) + (fk + gk*u(:,k))*Ts;...
            umin <= u(:,k), u(:,k) <= umax...
            ];
        tk = tk+Ts;
    end
    obj = u*u';
    options = sdpsettings('verbose',0,'solver','ipopt');
    optimize([constr, x(:,1)==xsim(:,i)],obj,options)
    % open-loop input sequence
    uval = value(u);
    % closed-loop input data
    usim = [usim, uval(1)]
    
    % simulating the system with uval
    fval = replace(f,xsym,xsim(:,i));
    gval = replace(g,xsym,xsim(:,i));
    x_next = xsim(:,i) + (fval + gval*usim(end))*dt;
    
    xsim = [xsim, x_next];
end

%%
subplot(2,1,1)
plot(0:dt:T,usim,'linewidth',2)
xlabel('$t$ (sec)','interpreter','latex')
ylabel('$u$ (N)','interpreter','latex')
grid on

subplot(2,1,2)
plot(0:dt:T,xsim(1,2:end),'linewidth',2)
hold on
plot(0:dt:T,xsim(2,2:end),'linewidth',2)
hold on
plot(0:dt:T,xsim(3,2:end),'linewidth',2)
hold on
plot(0:dt:T,xsim(4,2:end),'linewidth',2)
leg = legend('$d$ (m)','$\theta$ (rad)','$v$ (m/s)','$\dot{\theta}$ (rad/s)');
set(leg,'interpreter','latex')
grid on

xlabel('$t$ (sec)','interpreter','latex')
ylabel('states value','interpreter','latex')

