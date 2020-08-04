% This script is aming at using control Lyapunov function to design
% controller

%% System dynamics
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);

x = [x1; x2; x3; x4];
% create time t
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

% define control input
u = sdpvar(1,1);

% discrete time
dt = 0.005;

% time horizon
T = 1;

% initial condition
x0 = [0.18; -1.6; -0.15; 2.3];

%% load storage function
file = 'result.mat';
Vcoeff = loadVcoeff(file);

Vdeg = 4;
Vmonom_yalmip = monolist([x; t], Vdeg);
Vval_yalmip = Vcoeff'*Vmonom_yalmip;
DVvalx = jacobian(Vval_yalmip, x);
DVvalt = jacobian(Vval_yalmip, t);

%% Main algorithm
xsim = [];
xsim = [xsim, x0];
uopt = [];
time_line = 0:dt:T;
for i = 1:length(time_line)
    % plug in the measured state and time
    DVx = replace(replace(DVvalx, x, xsim(:,i)),t,time_line(i));
    DVt = replace(replace(DVvalt, x, xsim(:,i)),t,time_line(i));
    
    options = sdpsettings('verbose',1,'solver','mosek');
    fval = replace(f,x,xsim(:,i));
    gval = replace(g,x,xsim(:,i));
    constr = DVt + DVx*(fval + gval*u) <= 0; 
    obj = u'*u;
    optimize(constr,obj,options)
    % optimal value of u
    uval = value(u);
    uopt = [uopt, uval];
    
    % simulating the system with uval
    x_next = [xsim(1,i) + dt*xsim(2,i);...
        
              xsim(2,i) + dt*(- 10.6560*xsim(1,i)^3 + 11.5309*xsim(1,i)^2*xsim(3,i) + 7.8850*xsim(1,i)*xsim(3,i)^2 + 0.7972*xsim(2,i)^2*xsim(3,i) ...
              + 0.8408*xsim(2,i)*xsim(3,i)*xsim(4,i) + 21.0492*xsim(3,i)^3 + 0.4204*xsim(3,i)*xsim(4,i)^2 + 66.5225*xsim(1,i) - 24.5110*xsim(3,i)+ ...
              uval*(-10.0959*xsim(3,i)^2 + 44.2521));...
              
              xsim(3,i) + dt*xsim(4,i);...
              
              xsim(4,i) + dt*(10.9955*xsim(1,i)^3 - 48.9151*xsim(1,i)^2*xsim(3,i) - 6.4044*xsim(1,i)*xsim(3,i)^2 - 2.3955*xsim(2,i)^2*xsim(3,i) ...
              - 1.5943*xsim(2,i)*xsim(3,i)*xsim(4,i) - 51.9088*xsim(3,i)^3 - 0.7971*xsim(3,i)*xsim(4,i)^2 - 68.6419*xsim(1,i) + 103.9783*xsim(3,i)+ ...
              uval*(37.8015*xsim(3,i)^2 - 83.9120))];
          
    xsim = [xsim, x_next];
end
%%
save simulation.mat uopt xsim dt T

%%
subplot(2,1,1)
plot(0:dt:T,uopt,'linewidth',2)
xlabel('$t$ (sec)','interpreter','latex')
ylabel('$u$ (Nm)','interpreter','latex')
grid on

subplot(2,1,2)
plot(0:dt:T,xsim(1,2:end),'linewidth',2)
hold on
plot(0:dt:T,xsim(2,2:end),'linewidth',2)
hold on
plot(0:dt:T,xsim(3,2:end),'linewidth',2)
hold on
plot(0:dt:T,xsim(4,2:end),'linewidth',2)
leg = legend('$\theta_1$ (rad/s)','$\dot{\theta}_1$ (rad/s)',...
    '$\theta_2$ (rad/s)','$\dot{\theta}_2$ (rad/s)');
set(leg,'interpreter','latex')

xlabel('$t$ (sec)','interpreter','latex')
ylabel('states value','interpreter','latex')
grid on
