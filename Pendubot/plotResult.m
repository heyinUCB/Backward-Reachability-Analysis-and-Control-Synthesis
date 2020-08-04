plotaxes1 = [-3 3 -4.5 4.5];
plotaxes2 = [-0.5 0.5 -0.5 0.5];
% alpha_use = alpha_list(end);
gamma_use = gamma_list(end);

V0val = subs(Vval,t,0);
VTval = subs(Vval,t,T);

%%
subplot(1,2,1)
x2x4 = [0;0];
V013 = subs(V0val,[x2;x4],x2x4);
pcontour(V013, gamma_use, plotaxes2,'r',[500, 500]);
hold on

VT13 = subs(VTval,[x2;x4],x2x4);
pcontour(VT13, gamma_use, plotaxes2,'m',[500, 500]);
hold on

rT13 = subs(rT,[x2;x4],x2x4);
pcontour(rT13, 0, plotaxes2,'b',[500, 500]);
hold on

xlabel('$\theta_1$ (rad)','interpreter','latex')
ylabel('$\theta_2$ (rad)','interpreter','latex')
axis equal
grid on
legend('\Omega_{t_0,\gamma}^V','\Omega_{T,\gamma}^V','\Omega_{0}^{r_T}')
title('$\theta_1 - \theta_2$ plane','interpreter','latex')

%% 
subplot(1,2,2)
x1x3 = [0;0];
V024 = subs(V0val,[x1;x3],x1x3);
pcontour(V024, gamma_use, plotaxes1,'r',[500, 500]);
hold on

VT24 = subs(VTval,[x1;x3],x1x3);
pcontour(VT24, gamma_use, plotaxes1,'m',[500, 500]);
hold on

rT24 = subs(rT,[x1;x3],x1x3);
pcontour(rT24, 0, plotaxes1,'b',[500, 500]);
hold on

xlabel('$\dot{\theta}_1$ (rad/s)','interpreter','latex')
ylabel('$\dot{\theta}_2$ (rad/s)','interpreter','latex')
axis equal
grid on
title('$\dot{\theta}_1 - \dot{\theta}_2$ plane','interpreter','latex')
