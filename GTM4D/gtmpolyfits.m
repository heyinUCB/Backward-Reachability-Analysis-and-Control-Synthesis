% gtmpolyfits.m
%
% This function uses least squares to compute polynomial 
% fits needed for the 4-state longitudinal GTM model. 
% Polynomials fits are computed for trigonometric functions, 
% 1/V, aerodynamic tables (Wind Axes), and the engine model.
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
% 12/20/2009 PJS -- Update to convert to SI units

%--------------------------------------------------------------
% Setup
%--------------------------------------------------------------

% PlotsIdx are indices of figures to plot.  Set PlotsIdx=[] for no plots
%PlotsIdx = [4 5];
PlotsIdx = [];

% Load aerodynamic tables from GTM model
load aerodat;

% Conversion factors
d2r       = pi/180;         % Degrees to Radians [ rad/deg ]
r2d       = 180/pi;         % Radians to Degrees [ deg/rad ]

% Physical Parameters / Constant for EoM
g         = 9.8067;       % Gravitational accel [m/sec^2]
rho0      = 1.2240;       % Air density [ kg/m^3 ]
Sref      = 0.5483;       % Wing reference area [ m^2 ]
mass      = 22.498;       % Mass [kg]
cbar      = 0.2790;       % Mean chord [m]
Iyy       = 5.7676;       % Princ. moment of inertia [kg-m^2]
pratio    = 1.0;          % Pressure / Sea Level Pressure [unitless]

% Create polynomial variables
pvar alpha theta V q qhat de dthr

%--------------------------------------------------------------
% Taylor series approximation for trig functions
%--------------------------------------------------------------

% Fits for sin/cos of alpha and theta 
sinalp = alpha - alpha^3/6;
cosalp = 1 - alpha^2/2;

singam = (theta-alpha) - (theta-alpha)^3/6;
cosgam = 1 - (theta-alpha)^2/2;

if ismember(1,PlotsIdx) 
    % Plot poly fit with original curve 
    adata = linspace(-pi,pi);
    figure(1); 
    plot(adata*r2d,sin(adata),'b',...
        adata*r2d,double(subs(sinalp,alpha,adata)),'r--');
    xlabel('theta (deg)'); 
    legend('sin(theta)','Poly Fit');
end    
if ismember(2,PlotsIdx) 
    % Plot poly fit with original curve 
    figure(2); 
    plot(adata*r2d,cos(adata),'b',...
        adata*r2d,double(subs(cosalp,alpha,adata)),'r--');
    xlabel('theta (deg)'); 
    legend('cos(theta)','Poly Fit');
end

%--------------------------------------------------------------
% Least squares fit of 1/V over the range [Vmin, Vmax] m/sec
%--------------------------------------------------------------

pvar c0 c1 c2;
c2 = 0;  % Use a linear fit
Vinv = c0+c1*V+c2*V^2;
Vmin = 30; 
Vmax = 60; 
Vdata = linspace(Vmin,Vmax,100);
Vinv = pdatafit(Vinv,V,Vdata,1./Vdata);

if ismember(3,PlotsIdx) 
    % Plot poly fit with original curve
    figure(3); 
    plot(Vdata,1./Vdata,'b',Vdata,double(subs(Vinv,V,Vdata)),'r--');
    xlabel('V, ft/sec');
    legend('1/V','Poly Fit');    
end

%--------------------------------------------------------------
% Least squares fits to basic aero tables
%--------------------------------------------------------------

% Get aero data evaluated at beta=0
bidx = find(Aero.C6_bas.beta==0);
adata = Aero.C6_bas.alpha;
Cxadata = Aero.C6_bas.data(:,bidx,1);
Czadata = Aero.C6_bas.data(:,bidx,3);
Cmadata = Aero.C6_bas.data(:,bidx,5);
Na = length(adata);

% Transform from body to wind axes (Cx/Cz to CL/CD)
CLadata = Cxadata.*sin(adata*d2r) - Czadata.*cos(adata*d2r);
CDadata = -Cxadata.*cos(adata*d2r) - Czadata.*sin(adata*d2r);

% Fit CLa(alpha,beta=0) 
W = ones(Na,1);  W( (adata <= 10) & (adata>=0) ) = 4;
pvar c0 c1 c2 c3;
CLa = c0 + c1*alpha + c2*alpha^2 + c3*alpha^3; 
CLa = pdatafit(CLa,alpha,adata,CLadata,W);
if ismember(4,PlotsIdx)  
    % Plot poly fit with original curve
    figure(4); 
    plot(adata,CLadata,'b',adata,double(subs(CLa,alpha,adata)),'r--');
    xlabel('alpha (deg)'); ylabel('CLa');
    legend('CLa','Poly Fit');
end

% Fit CDa(alpha,beta=0)
W = ones(Na,1); W( (adata <= 10) & (adata>=0) ) = 20;
pvar c0 c1 c2 c3;
CDa = c0 + c1*alpha + c2*alpha^2 + c3*alpha^3; 
CDa = pdatafit(CDa,alpha,adata,CDadata,W);
if ismember(5,PlotsIdx)  
    % Plot poly fit with original curve
    figure(5); 
    plot(adata,CDadata,'b',adata,double(subs(CDa,alpha,adata)),'r--');
    xlabel('alpha (deg)'); ylabel('CDa');
    legend('CDa','Poly Fit','Location','Best');
end

% Fit Cma(alpha,beta=0)
W = ones(Na,1);  W( (adata <= 10) & (adata>=0) ) = 30;
pvar c0 c1 c2 c3;
Cma = c0 + c1*alpha + c2*alpha^2 + c3*alpha^3; 
Cma = pdatafit(Cma,alpha,adata,Cmadata,W);
if ismember(6,PlotsIdx)  
    % Plot poly fit with original curve
    figure(6); 
    plot(adata,Cmadata,'b',adata,double(subs(Cma,alpha,adata)),'r--');
    xlabel('alpha (deg)'); ylabel('Cma');
    legend('Cma','Poly Fit');
end

% Data/fits have alpha in degs but the model has alpha in rads
CLa = subs(CLa,alpha,alpha*r2d);
CDa = subs(CDa,alpha,alpha*r2d);
Cma = subs(Cma,alpha,alpha*r2d);

%--------------------------------------------------------------
% Least squares fits to elevator aero tables
%--------------------------------------------------------------

% Get aero data evaluated at beta=0 and stab = 0
bidx = find(Aero.dC3_ele.beta==0);
sidx = find(Aero.dC3_ele.stab==0);
adata = Aero.dC3_ele.alpha;
edata = Aero.dC3_ele.elev;

[agrid,egrid] = meshgrid(adata,edata);
Xdata = [agrid(:) egrid(:)]';

Cxdedata = squeeze(Aero.dC3_ele.data(:,bidx,sidx,:,1))';
Czdedata = squeeze(Aero.dC3_ele.data(:,bidx,sidx,:,2))';
Cmdedata = squeeze(Aero.dC3_ele.data(:,bidx,sidx,:,3))';

% Transform from body to wind axes (Cx/Cz to CL/CD)
sinadata = repmat( sin(adata(:)*d2r)' , [size(Cxdedata,1) 1] );
cosadata = repmat( cos(adata(:)*d2r)' , [size(Cxdedata,1) 1] );

CLdedata = Cxdedata.*sinadata - Czdedata.*cosadata;
CDdedata = -Cxdedata.*cosadata - Czdedata.*sinadata;

% Fit CLde(alpha,beta=0,stab=0,de)
W = ones(size(agrid)); 
W( (agrid>=0) & (agrid <= 10) & (abs(egrid)<=10) ) = 4;
W( (agrid>=0) & (agrid <= 10) & (egrid==0) ) = 10;
pvar c0 c1 c2 c3 c4 c5;
CLde = c0 + c1*de + c2*alpha + c3*de^2 + c4*de*alpha + c5*alpha^2; 
CLde = pdatafit(CLde,[alpha;de],Xdata,CLdedata(:),W(:));
if ismember(7,PlotsIdx)  
    % Plot poly fit with original curve
    figure(7); 
    fit = double(subs(CLde,[alpha;de],Xdata));
    fit = reshape(fit,size(CLdedata));
    mesh(adata,edata,CLdedata); hold on;
    mh = surf(adata,edata,fit); hold off;
    set(mh,'CDataMapping','direct')
    xlabel('alpha (deg)'); ylabel('de (deg)'); zlabel('CLde');
end

% Fit CDde(alpha,beta=0,stab=0,de)
W = ones(size(agrid)); 
W( (agrid>=0) & (agrid <= 10) & (abs(egrid)<=10) ) = 4;
W( (agrid>=0) & (agrid <= 10) & (egrid==0) ) = 10;
pvar c0 c1 c2 c3 c4 c5;
CDde = c0 + c1*de + c2*alpha + c3*de^2 + c4*de*alpha + c5*alpha^2; 
CDde = pdatafit(CDde,[alpha;de],Xdata,CDdedata(:),W(:));
if ismember(8,PlotsIdx) 
    % Plot poly fit with original curve
    figure(8); 
    fit = double(subs(CDde,[alpha;de],Xdata));
    fit = reshape(fit,size(CDdedata));
    mesh(adata,edata,CDdedata); hold on;
    mh = surf(adata,edata,fit); hold off;
    set(mh,'CDataMapping','direct')
    xlabel('alpha (deg)'); ylabel('de (deg)'); zlabel('CDde');
end

% Fit Cmde(alpha,beta=0,stab=0,de)
% Polynomial form is such that Cmqde(alpha,0,0,de=0) = 0
W = ones(size(agrid)); 
W( (agrid>=0) & (agrid <= 10) & (abs(egrid)<=10) ) = 10;
W( (agrid>=0) & (agrid <= 10) & (egrid==0) ) = 20;
pvar c0 c1 c2;
Cmde = c1*de+c2*de*alpha; 
Cmde = pdatafit(Cmde,[alpha;de],Xdata,Cmdedata(:),W(:));
if ismember(9,PlotsIdx)  
    % Plot poly fit with original curve
    figure(9); 
    fit = double(subs(Cmde,[alpha;de],Xdata));
    fit = reshape(fit,size(Cmdedata));
    mesh(adata,edata,Cmdedata); hold on;
    mh = surf(adata,edata,fit); hold off;
    set(mh,'CDataMapping','direct')
    xlabel('alpha (deg)'); ylabel('de (deg)'); zlabel('Cmde');     
end

% Data has alpha in degs but the model has alpha in rads
CLde = subs(CLde,alpha,alpha*r2d);
CDde = subs(CDde,alpha,alpha*r2d);
Cmde = subs(Cmde,alpha,alpha*r2d);

% Data has de in degs but the model has de in rads
CLde = subs(CLde,de,de*r2d);
CDde = subs(CDde,de,de*r2d);
Cmde = subs(Cmde,de,de*r2d);

%--------------------------------------------------------------
% Least squares fits to pitch rate aero tables
%--------------------------------------------------------------

% Get aero data  (qhat = q*cbar/2/V)
adata = Aero.dC3_q.alpha;
qdata = Aero.dC3_q.qhat;

[agrid,qgrid] = meshgrid(adata,qdata);
Xdata = [agrid(:) qgrid(:)]';

Cxqdata = squeeze(Aero.dC3_q.data(:,:,1))';
Czqdata = squeeze(Aero.dC3_q.data(:,:,2))';
Cmqdata = squeeze(Aero.dC3_q.data(:,:,3))';

% Transform from body to wind axes (Cx/Cz to CL/CD)
sinadata = repmat( sin(adata(:)*d2r)' , [size(Cxqdata,1) 1] );
cosadata = repmat( cos(adata(:)*d2r)' , [size(Cxqdata,1) 1] );

CLqdata = Cxqdata.*sinadata - Czqdata.*cosadata;
CDqdata = -Cxqdata.*cosadata - Czqdata.*sinadata;

% Fit CLq(alpha,qhat)
W = ones(size(agrid)); 
W( (agrid>=0) & (agrid <= 10) & (abs(qgrid)<=0.002) ) = 4;
W( (agrid>=0) & (agrid <= 10) & (qgrid==0) ) = 10;
pvar c0 c1 c2 c3 c4 c5;
CLq = c0 + c1*qhat + c2*alpha + c3*qhat^2 + c4*qhat*alpha + c5*alpha^2; 
CLq = pdatafit(CLq,[alpha;qhat],Xdata,CLqdata(:),W(:));
if ismember(10,PlotsIdx) 
    % Plot poly fit with original curve
    figure(10); 
    fit = double(subs(CLq,[alpha;qhat],Xdata));
    fit = reshape(fit,size(CLqdata));
    mesh(adata,qdata,CLqdata); hold on;
    mh = surf(adata,qdata,fit); hold off;
    set(mh,'CDataMapping','direct')
    xlabel('alpha (deg)'); ylabel('qhat (unitless)'); zlabel('CLq');
end

% Fit CDq(alpha,qhat)
W = ones(size(agrid)); 
W( (agrid>=0) & (agrid <= 10) & (abs(qgrid)<=0.002) ) = 4;
W( (agrid>=0) & (agrid <= 10) & (qgrid==0) ) = 10;
pvar c0 c1 c2 c3 c4 c5;
CDq = c0 + c1*qhat + c2*alpha + c3*qhat^2 + c4*qhat*alpha + c5*alpha^2; 
CDq = pdatafit(CDq,[alpha;qhat],Xdata,CDqdata(:),W(:));
if ismember(11,PlotsIdx) 
    % Plot poly fit with original curve
    figure(11); 
    fit = double(subs(CDq,[alpha;qhat],Xdata));
    fit = reshape(fit,size(CDqdata));
    mesh(adata,qdata,CDqdata); hold on;
    mh = surf(adata,qdata,fit); hold off;
    set(mh,'CDataMapping','direct')
    xlabel('alpha (deg)'); ylabel('qhat (unitless)'); zlabel('CDq');
end

% Fit Cmq(alpha,qhat)
% Polynomial form is such that Cmq(alpha,0) = 0
W = ones(size(agrid)); 
W( (agrid>=0) & (agrid <= 10) & (abs(qgrid)<=0.002) ) = 4;
W( (agrid>=0) & (agrid <= 10) & (qgrid==0) ) = 10;
pvar c0 c1;
Cmq = c1*qhat; 
Cmq = pdatafit(Cmq,[alpha;qhat],Xdata,Cmqdata(:),W(:));
if ismember(12,PlotsIdx) 
    % Plot poly fit with original curve
    figure(12); 
    fit = double(subs(Cmq,[alpha;qhat],Xdata));
    fit = reshape(fit,size(Cmqdata));
    mesh(adata,qdata,Cmqdata); hold on;
    mh = surf(adata,qdata,fit); hold off;
    set(mh,'CDataMapping','direct')
    xlabel('alpha (deg)'); ylabel('qhat (unitless)'); zlabel('Cmq');
end

% Data has alpha in degs but the model has alpha in rads
CLq = subs(CLq,alpha,alpha*r2d);
CDq = subs(CDq,alpha,alpha*r2d);
Cmq = subs(Cmq,alpha,alpha*r2d);

% Data uses normalized pitch rate. 
% Sub in qhat = cbar*q/2/V to get dependence on q and V.
CLq = subs(CLq,qhat,cbar*q/2*Vinv);
CDq = subs(CDq,qhat,cbar*q/2*Vinv);
Cmq = subs(Cmq,qhat,cbar*q/2*Vinv);

%--------------------------------------------------------------
% CEP Paper Modification
% XXX Abhi :  Modifed to get Cubic qhat dependence (for CEP Paper)
if exist('Cmq3','var')
    Cmq = Cmq + (Cmq3*cbar*q/2*Vinv)^3;
end


%--------------------------------------------------------------
% Least squares fits to engine model 
% (Thrust in Newtons)
%--------------------------------------------------------------

% T(dthr): GTM model uses a high order polynomial function
% Coeff_thrust and Coeff_throttle are set in the Jetcat P70 mask 
% initialization of the Engines subsystem in the 6DOF GTM model.
Coeff_thrust   = [0.000028864242464  -0.001914835880159 ...
    0.105916926095459  -1.251646340183624 ];
Coeff_throttle = [0.000033035538475  -0.008354385314499 ...
    1.191126728328792 31.518959847101669];
RPMref = Coeff_throttle*[dthr^3; dthr^2; dthr; 1];
Tdthr = Coeff_thrust*[RPMref^3; RPMref^2; RPMref;1];
Tdthr = pratio*Tdthr;

% Thrust model is in lbs. Convert to Newtons.
lbs2N      = 4.44822162;    % pound force to Newton [ N/lbs ]
Tdthr = Tdthr*lbs2N;      
Xdata = linspace(0,100);
Ydata = double(subs(Tdthr,dthr,Xdata));
        
% Fit T(dthr) with a lower order polynomial
pvar c0 c1 c2 c3;
Tdthr =  c0+c1*dthr+c2*dthr^2+c3*dthr^3;
W = ones(length(Xdata),1); W( Xdata<30 ) = 10;
Tdthr = pdatafit(Tdthr,dthr,Xdata,Ydata,W);
if ismember(13,PlotsIdx) 
    % Plot poly fit with original curve
    figure(13); 
    fit = double(subs(Tdthr,dthr,Xdata));
    plot(Xdata,Ydata,'b',Xdata,fit,'r--');    
    xlabel('dthr (percent)'); ylabel('Tdthr (N)'); 
end


