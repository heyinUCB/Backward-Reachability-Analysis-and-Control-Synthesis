% This script is aming at using control Lyapunov function to design
% controller

% discrete time
dt = 0.005;

% initial condition
x0 = [-0.35; 2.6; 0.35; -4];

%% load storage function
load u1T4_step5.mat

%% Main algorithm
xsim = [];
xsim = [xsim, x0];
uopt = [];
Vtx = [];
time_line = 0:dt:T;
for i = 1:length(time_line)
    % plug in the measured state and time    
    fval = subs(f,x,xsim(:,i));
    gval = subs(g,x,xsim(:,i));
    Vtx = [Vtx; subs(subs(Vval,x,xsim(:,i)),t,time_line(i))];
    ucmd = subs(subs(uval, x, xsim(:,i)),t,time_line(i));
    uopt = [uopt, ucmd];
    
    % simulating the system with uval
%     x_next = [xsim(1,i) + dt*xsim(2,i);...
%         
%               xsim(2,i) + dt*(- 10.6560*xsim(1,i)^3 + 11.5309*xsim(1,i)^2*xsim(3,i) + 7.8850*xsim(1,i)*xsim(3,i)^2 + 0.7972*xsim(2,i)^2*xsim(3,i) ...
%               + 0.8408*xsim(2,i)*xsim(3,i)*xsim(4,i) + 21.0492*xsim(3,i)^3 + 0.4204*xsim(3,i)*xsim(4,i)^2 + 66.5225*xsim(1,i) - 24.5110*xsim(3,i)+ ...
%               uval*(-10.0959*xsim(3,i)^2 + 44.2521));...
%               
%               xsim(3,i) + dt*xsim(4,i);...
%               
%               xsim(4,i) + dt*(10.9955*xsim(1,i)^3 - 48.9151*xsim(1,i)^2*xsim(3,i) - 6.4044*xsim(1,i)*xsim(3,i)^2 - 2.3955*xsim(2,i)^2*xsim(3,i) ...
%               - 1.5943*xsim(2,i)*xsim(3,i)*xsim(4,i) - 51.9088*xsim(3,i)^3 - 0.7971*xsim(3,i)*xsim(4,i)^2 - 68.6419*xsim(1,i) + 103.9783*xsim(3,i)+ ...
%               uval*(37.8015*xsim(3,i)^2 - 83.9120))];
    x_next = xsim(:,i)+dt*(fval+gval*ucmd);
    xsim = [xsim, x_next];
end

%%
subplot(3,1,1)
plot(0:dt:T,uopt,'r','linewidth',2)
hold on
xlabel('$t$ (sec)','interpreter','latex')
ylabel('$u$ (Nm)','interpreter','latex')
grid on

subplot(3,1,2)
plot(0:dt:T,xsim(2,2:end),'b','linewidth',2)
hold on
plot(0:dt:T,xsim(4,2:end),'m','linewidth',2)
hold on
leg = legend('$\dot{\theta}_1$ (rad/s)','$\dot{\theta}_2$ (rad/s)');
set(leg,'interpreter','latex')

xlabel('$t$ (sec)','interpreter','latex')
ylabel('$\dot{\theta_1}$, $\dot{\theta_2}$','interpreter','latex')
grid on

subplot(3,1,3)
plot(0:dt:T,xsim(1,2:end),'r','linewidth',2)
hold on
plot(0:dt:T,xsim(3,2:end),'g','linewidth',2)
hold on
leg = legend('$\theta_1$ (rad)','$\theta_2$ (rad)');
set(leg,'interpreter','latex')

xlabel('$t$ (sec)','interpreter','latex')
ylabel('$\theta_1$, $\theta_2$','interpreter','latex')
grid on
