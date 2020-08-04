% This code is trying to scale the units of state from the ones used for
% computation to [m/s, deg, deg/s, deg];
Dmax = diag([20 20 50 20]);
Vval = Vval_list(end);
gamma_use = gamma_list(end);

d2r = pi/180;
r2d = 180/pi;

xeq = [45.0; 0.04924*r2d; 0.0; 0.04924*r2d];

plotaxes1 = [-40 40 -280 220];
plotaxes2 = [35 55 -28 35];

V0val = subs(Vval,t,0);
V0val = subs(V0val,x,inv(Dmax)*x);
V0val = subs(V0val,x,x-xeq);

VTval = subs(Vval,t,T);
VTval = subs(VTval,x,inv(Dmax)*x);
VTval = subs(VTval,x,x-xeq);

rT = subs(rT,x,inv(Dmax)*x);
rT = subs(rT,x,x-xeq);

%% x2 vs. x3
subplot(1,2,1)
x1x4 = [45; 0.04924*r2d];
V023 = subs(V0val,[x1;x4],x1x4);
pcontour(V023, gamma_use, plotaxes1,'g',[500 500]);
hold on
VT23 = subs(VTval,[x1;x4],x1x4);
pcontour(VT23, gamma_use, plotaxes1,'m',[500 500]);
hold on
rT23 = subs(rT,[x1;x4],x1x4);
pcontour(rT23, 0, plotaxes1,'b',[500 500]);
xlabel('$x_2$ (deg)','interpreter','latex')
ylabel('$x_3$ (deg/s)','interpreter','latex')
title('$x_2 - x_3$ plane','interpreter','latex')
grid on

%% x1 vs. x4
x2x3 = [0.04924*r2d; 0];
subplot(1,2,2)
V014 = subs(V0val,[x2;x3],x2x3);
pcontour(V014, gamma_use, plotaxes2,'g',[500,500]);
hold on
VT14 = subs(VTval,[x2;x3],x2x3);
pcontour(VT14, gamma_use, plotaxes2,'m',[500,500]);
hold on
rT14 = subs(rT,[x2;x3],x2x3);
pcontour(rT14, 0, plotaxes2,'b',[500,500]);
xlabel('$x_1$ (m/s)','interpreter','latex')
ylabel('$x_4$ (deg)','interpreter','latex')
title('$x_1 - x_4$ plane','interpreter','latex')
grid on
leg = legend('$\Omega_{t_0,\gamma}^V$','$\Omega_{T,\gamma}^V$',...
    '$\Omega_0^{r_T}$');
set(leg,'interpreter','latex')
