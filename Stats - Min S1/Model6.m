function [M] = Model6(t,y,vDNAin0)
load('Model6_mpars_min.mat');
mpar = mpars;

% Calculate vDNAtot at current t
[vDNAtot] = Calc_vDNA(t,vDNAin0);

% Model parameters
ks1 = mpar(1);
Km1 = mpar(2);
kd1 = mpar(3);
Km3 = mpar(4);
Km5 = mpar(5);

ks2 = mpar(6);
Km2 = mpar(7);
kd2 = mpar(8);
Km4 = mpar(9);
Km6 = mpar(10);

ks_C = mpar(11);
ks_P = mpar(12);
kex = mpar(13);

% State variables
Protein1 = y(1);
Protein2 = y(2);
Capsid = y(3);
Particle = y(4);
Virus = y(5);

dProtein1 = ks1*(vDNAtot/(Km1+vDNAtot))*(Km3/(Km3+Protein1)) - kd1*...
    (vDNAtot/(Km5+vDNAtot))*Protein1 - ks_C*vDNAtot*Protein1;
 
dProtein2 = ks2*(vDNAtot/(Km2+vDNAtot))*(Km4/(Km4+Protein2)) - kd2*...
    (vDNAtot/(Km6+vDNAtot))*Protein2 - ks_P*Capsid*Protein2;
 
dCapsid = ks_C*vDNAtot*Protein1 - ks_P*Capsid*Protein2;
 
dParticle = ks_P*Capsid*Protein2 - kex*Particle;
 
dVirus = kex*(1/500)*Particle;

M = [dProtein1,dProtein2,dCapsid,dParticle,dVirus]';
end