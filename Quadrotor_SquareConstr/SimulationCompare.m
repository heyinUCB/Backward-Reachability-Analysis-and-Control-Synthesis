%% System dynamics
x = mpvar('x',[4,1]);
x1 = x(1); % distance of the cart
x2 = x(2); % angle of the rod
x3 = x(3); % speed of the cart
x4 = x(4); % angular velocity of the rod

%% create time t
pvar t

%% Define the system
f3 = 0.11707*x2^5 + 0.035908*x2^3*x4^2 - 1.6032*x2^3 - 0.17201*x2*x4^2 + 3.0313*x2;
f4 = 0.24902*x2^5 + 0.13049*x2^3*x4^2 - 5.6188*x2^3 - 0.29147*x2*x4^2 + 23.9892*x2;
f = [x3; x4; f3; f4];
g3 = 0.02905*x2^4 - 0.11289*x2^2 + 0.3955;
g4 = 0.096371*x2^4 - 0.54277*x2^2 + 0.7831;
g = [0; 0; g3; g4];

cleanpoly(fcubicR,10^-5)
[xtraj]=psim(fscaleR,x,pointInit(i,:)',T);