% This file generates the required GTM models for analysis presented in the
% CEP paper. Specifically, the file generates both 2-state short-period
% polynomial models and longitudinal 4-state models including off-nominal
% values of Cmq3 derivative.
%

% ======================================================
% (1) Generate Polynomial Models

d2r       = pi/180;         % Degrees to Radians [ rad/deg ]
r2d       = 180/pi;         % Radians to Degrees [ deg/rad ]


%--------- Cmq3 Values
% XXX gtmpolyfits has been modified Line 352 to add Cmq3
Cmq3Values = [ 0 ; 1.7455*57.29];

% Initialize model /trim storage variable
fgtm = cell(length(Cmq3Values),5);

for i1 = 1 :length(Cmq3Values)
    
    Cmq3 = Cmq3Values(i1);
    
    %--------------------------------------------------------------
    %  Create 4-state longitudinal poly model:
    %  The poly fits are done in gtmpolyfits which is called by
    %  gtmlongpoly. See gtmlongpoly for info on the 4-state model
    %--------------------------------------------------------------
    gtmlongpoly;
    
    %--------------------------------------------------------------
    % Compute level-flight trim condition for longitudinal model
    %--------------------------------------------------------------
    
    Vtrim = 45;   % Air Speed (m/sec) ~ 150 ft/s
    x0 = [45; 0.05;0;-0.05];
    u0 = [16;3*d2r];
    [xtp,utp]=ptrim(flong,x,u,x0,u0,[V-Vtrim; alpha-theta]);
    
    %--------------------------------------------------------------
    % Compute short period trim condition
    %--------------------------------------------------------------
    xtsp = xtp(2:3);
    utsp = utp(2);
    
    %--------------------------------------------------------------
    %  Nonlinear Short-period Model
    %  Sub in trim velocity, pitch, and thrust to obtain 2-state
    %  1-input short-period model
    %--------------------------------------------------------------
    fsp = subs(flong,[V;theta;dthr],[xtp(1);xtp(4);utp(1)]);
    fsp = fsp(2:3);
    xsp = [alpha;q];
    usp = de;
    
    %--------------------------------------------------------------
    % Shift the equilibrium
    % Model is in x = m/s , rad ;  u = rad , %
    %--------------------------------------------------------------
    f4 = subs(flong, [x;u] , [x+xtp;u+utp]);
    f2 = subs(fsp,[alpha;q;de], [alpha;q;de]+[xtp(2);xtp(3);utp(2)]);
    
    % store model and trim values
    fgtm{i1,1} = f4;    % 4-state model
    fgtm{i1,2} = f2;    % 2-state model
    fgtm{i1,3} = xtp;   % trimmed state values
    fgtm{i1,4} = utp;   % trimmed input values
    fgtm{i1,5} = [x;u]; % state / input variable
    
end


%--------------------------------------------------------------------------
% Extract short-period model with NO cubic pitch damping nonlineairty
xtrim2 = fgtm{1,3}; utrim2 = fgtm{1,3};
f2model = fgtm{1,2};
f2sim = subs(f2model, de,0);
f2nom = cleanpoly(f2sim,1e-08);

%--------------------------------------------------------------------------
% Extract short-period model with cubic nonlineairty in pitch damping
xtrimC = fgtm{2,3}; utrimC = fgtm{2,3};
f2modelCubic = fgtm{2,2};
f2simCubic = subs(f2modelCubic, de,0);
f2Cubic = cleanpoly(f2simCubic,1e-08);
xtrimC = fgtm{2,3}; utrimC = fgtm{2,3};

%--------------------------------------------------------------------------
% Extract 4-state longitudinal model and form Closed Loop
xtrim4 = fgtm{1,3}; utrim4 = fgtm{1,3};
f4model = fgtm{1,1};
f4sim = subs(f4model,u,[0;0]);
f4nom = cleanpoly(f4sim,1e-14);

% 4-state Polynomial Closed Loop Model
Kq = -4/r2d;
f4nomClp = subs(f4model,[de;dthr],[-Kq*q;0]);
f4Clp   = cleanpoly(f4nomClp,10^-12);

