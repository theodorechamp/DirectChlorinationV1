clear

%%%%%%%%%%
% Inputs %
%%%%%%%%%%
P = 101.325; % kPa
T = 333; % K
mc2h40 = 7.1; % kg/s, from MatBal
mcl20 = 20.28;
mc2h4cl20 = 1000;
mhcl0 = 0.0;
mc2h3cl30 = 0.0;
Vr = 100; %m^3


%%%%%%%%%%%%%
% Constants %
%%%%%%%%%%%%%
R = 8.314; % [L kPa / mol K]
mwc2h4 = 28.05; % [kg / kmol]
mwcl2 = 70.9;
mwc2h4cl2 = 98.96;
mwhcl = 36.46;
mwc2h3cl3 = 133.4;
rhoc2h4cl2 = 1.25; %g/cm^3

k1 = .132; %m3 mol-1 s-1
k2 = 0.0239; %m3 mol-2 s-1
k3 = 6.12*10^-9; %m3 mol-1 s-1

%%%%%%%%%
% Logic %
%%%%%%%%%

% Convert to molar flows [mols/s]
nc2h40 = mc2h40/mwc2h4*1000;
ncl20 = mcl20/mwcl2*1000;
nc2h4cl20 = mc2h4cl20/mwc2h4cl2*1000;
nhcl0 = mhcl0/mwhcl*1000;
nc2h3cl30 = mc2h3cl30/mwc2h3cl3*1000;

% Composition of inlet
ntot = nc2h40 + ncl20 + nc2h4cl20 + nhcl0 + nc2h3cl30;
xc2h40 = nc2h40/ntot;
xcl20 = ncl20/ntot;
xc2h4cl20 = nc2h4cl20/ntot;
xhcl0 = nhcl0/ntot;
xc2h3cl30 = nc2h3cl30/ntot;

% Volume of EDC

Vi= mc2h4cl20 / rhoc2h4cl2 / 1000; %m^3/s
% calculate concentrations in [mol/m3]
% Ci0 refers to inlet concentration, Ci refers to reactor and outlet
Cc2h40 = nc2h40/Vi;
Ccl20 = ncl20/Vi;
Cc2h4cl20 = nc2h4cl20/Vi;
Chcl0 = nhcl0/Vi;
Cc2h3cl30 = nc2h3cl30/Vi;

disp(Cc2h40)


%%%%%%%%%%%%%
% Equations %
%%%%%%%%%%%%%
%define rho and mw of C2H4Cl2 so the units when divided are mols/m3
mw = mwc2h4cl2; %g/mols
rho = rhoc2h4cl2*100^3; %g/m^3

syms Cc2h4 Ccl2 Vo Chcl Cc2h3cl3;
[Cc2h4, Ccl2, Vo, Chcl, Cc2h3cl3] = vpasolve([...
    Vi*Cc2h40 - k1*Cc2h4*Ccl2*Vr - k2*Cc2h4*Ccl2^2*Vr == Vo*Cc2h4...
    Vi*Ccl20 - k1*Cc2h4*Ccl2*Vr - 2*k2*Cc2h4*Ccl2^2*Vr - k3*rho/mw*Ccl2*Vr == Vo*Ccl2...
    Vi + k1*Cc2h4*Ccl2*Vr - k3*rho/mw*Ccl2*Vr == Vo...
    Vi*Chcl0 + k2*Cc2h4*Ccl2^2*Vr +  k3*rho/mw*Ccl2*Vr == Vo*Chcl...
    Vi*Cc2h3cl30 + k2*Cc2h4*Ccl2^2*Vr + k3*rho/mw*Ccl2*Vr == Vo*Cc2h3cl3],...
    [Cc2h4,Ccl2,Vo,Chcl, Cc2h3cl3], [300; 300; 1; 10; 10]);
disp(Cc2h4[1] + 

