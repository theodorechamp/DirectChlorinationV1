clear

%%%%%%%%%%
% Inputs %
%%%%%%%%%%
P = 101.325; % kPa
T = 333; % K
mc2h40 = 0.0; % kg/s, from MatBal
mcl20 = 0.0;
mc2h4cl20 = 0.0;
mhcl0 = 0.0;
mc2h3cl30 = 0.0;


%%%%%%%%%%%%%
% Constants %
%%%%%%%%%%%%%
R = 8.314; % [L kPa / mol K]
mwc2h4 = 28.05;
mwcl2 = 70.9;
mwc2h4cl2 = 98.96;
mwhcl = 36.46;
mwc2h3cl3 = 133.4;

%%%%%%%%%
% Logic %
%%%%%%%%%

% Convert to molar flows [kmols/s]
nc2h40 = mwc2h4 * mc2h40;
ncl20 = mwcl2 * mcl20;
nc2h4cl20 = mwc2h4cl2 * mc2h4cl20;
nhcl0 = mwhcl * mhcl0;
nc2h3cl30 = mwc2h3cl3 * mwc2h3cl30;

% Composition of inlet
ntot = nc2h4 + ncl20 + nc2h4cl20 + nhcl0 + nc2h3cl3;
xc2h40 = nc2h4/ntot;
xcl20 = ncl20/ntot;
xc2h4cl20 = nc2h4cl20/ntot;
xhcl0 = nhcl0/ntot;
xc2h3cl30 = nc2h3cl3/ntot;

% calculate concentrations in [kmols/L]
% Ci0 refers to inlet concentration, Ci refers to reactor and outlet
Ctot = P/(R*T);
Cc2h40 = Ctot*xc2h40;
Ccl2 = Ctot*xcl2;
Cc2h4cl20 = Ctot*xc2h4cl20;
Chcl0 = Ctot*xhcl0;
Cc2h3cl30 = Ctot*xc2h3cl30;

Vin = ntot*

%%%%%%%%%%%%%
% Equations %
%%%%%%%%%%%%%

syms x y;
Soln = vpasolve([y == x + 2, y == x^2], [x,y], [0 Inf; 0 Inf]);
disp(Soln.x)