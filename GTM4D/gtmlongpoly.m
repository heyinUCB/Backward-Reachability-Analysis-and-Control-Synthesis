% This function creates a 4-state polynomial model of the
% longitudinal GTM dynamics.  Least squares is used to 
% compute polynomial fits to trigonometric functions, 1/V,
% aerodynamic tables, and the engine model.
%
% States:
%  V = x(1);                   % True Air Speed, V (m/sec); 
%  alpha = x(2);               % Angle of Attack, alpha (rad); 
%  q = x(3);                   % Pitch Rate, q (rad/sec);
%  theta = x(4);               % Pitch Euler Angle, theta (rad)
%
% Inputs:
%  dthr = u(1);                % Throttle Position, dthr (Percent, 0 to 100)
%  de = u(2);                  % Elevator Deflection, de (rad)
%
% 10/5/2009 PJS

%--------------------------------------------------------------
% Polynomial Fits
%--------------------------------------------------------------
gtmpolyfits;

%--------------------------------------------------------------
% Engine Model
%--------------------------------------------------------------

% L and R engines are assumed to have equal setting.
cg_pos0 = [-4.7750  0  -0.9401];
engl_ang = [0 2.1460 1.6830];
engl_pos = [-4.32525  -1.1833  -0.6425];
engr_ang = [0 2.1460 -1.6830];
engr_pos = [-4.32525  1.1833  -0.6425];
deltaLEng = engl_pos - cg_pos0;

% Convert deltaLEng to SI units
ft2m      = 0.3048;         % ft to meters       [ m/ft ]
deltaLEng = deltaLEng*ft2m;

% Direction cosine vector (DCV) to transform thrust axis 
% to airplane body axis. 
DCV_LEng = [cosd(engl_ang(2))*cosd(engl_ang(3));...
           -sind(engl_ang(3));...
            sind(engl_ang(2))*cosd(engl_ang(3))];

% Engine Thrust: Body Axis Forces Forces and  Moments
neng = 2;
Tx = neng*DCV_LEng(1)*Tdthr;
Tz = neng*DCV_LEng(3)*Tdthr;
Tm = (Tx*deltaLEng(3)-Tz*deltaLEng(1));

%--------------------------------------------------------------
% Form polynomial Model
%--------------------------------------------------------------

% Combine Aero coef terms
CL = CLa + CLde + CLq;
CD = CDa + CDde + CDq;
Cm = Cma + Cmde + Cmq;

% Forces / Moments
qbar = 0.5*rho0*V^2;
L = qbar*Sref*CL;
D = qbar*Sref*CD;
My = qbar*Sref*cbar*Cm + Tm;

% State Equations
Vdot = ( -D - mass*g*singam + Tx*cosalp + Tz*sinalp  )/mass;
alphadot = ( -L + mass*g*cosgam -Tx*sinalp + Tz*cosalp )*Vinv/mass + q;
qdot = My/Iyy;
thetadot = q;

% Stack state, input, and state derivatives
flong = [Vdot; alphadot;  qdot; thetadot];
x = [V; alpha; q; theta];
u = [dthr; de];

