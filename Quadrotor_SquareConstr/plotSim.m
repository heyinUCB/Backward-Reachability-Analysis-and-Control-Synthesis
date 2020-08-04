%% Parameters
% number of initial points
npts = 200;

% Simulation time
T = 1;
% number of steps
% num_step = length(V0_cell);
num_step = 1;
pvar x1 x2

%% Define epsilon0
% Initial set epsilon0: [1,x1,x2]*R0*[1;x1;x2] <= 0
% Partition R0 into [a, b'; b, C].
C = [1/(0.1)^2, 0; 0, 1/(0.1)^2];

x0Center = [0.0; 0.0]; % center of the set epsilon0
b = -C*x0Center;
a = b'*inv(C)*b - 1;
R0 = [a, b'; b, C];

epsilon0 = [1, x1, x2]*R0*[1; x1; x2];

% Generate random points inside a ellipsoid as initial points
pointInit = PointInEllipsoid(R0, npts);

% Initialize
alpha_use_pre = 1.0;
ti = 0;
tf = num_step*T;

% Trajectory
xT = [];
yT = [];

% Final point
xF = [];
yF = [];

% figure('rend','painters','pos',[10 10 900 500])
for i = 1:npts
    Tspan = [ti, tf];
    dyn =@(t, x) sysTunnelDiode(t, x);
    
    [TOUT,YOUT] = ode45(dyn,Tspan,[pointInit(i,1); pointInit(i,2)]);
    xF = [xF; YOUT(end,1)];
    yF = [yF; YOUT(end,2)];
%     plot(YOUT(:,1), YOUT(:,2), 'b')
%     hold on
end


for j = 1:length(xF)
    plot(xF(j), yF(j),'.g')
    hold on
    axis equal
end

for i = 1:size(pointInit,1)
    plot(pointInit(i,1), pointInit(i,2),'.r')
    hold on
    axis equal
end