function [M] = Model1(t,y,vDNAin0)
load('Model1_mpars_min.mat');
mpar = mpars;

% Calculate vDNAtot at current t
[vDNAtot] = Calc_vDNA(t,vDNAin0);

% Model parameters
ks1 = mpar(1);
Km1 = mpar(2);
kd1 = mpar(3);

ks2 = mpar(4);
Km2 = mpar(5);
kd2 = mpar(6);

ks_C = mpar(7);
ks_P = mpar(8);
kex = mpar(9);

% State variables
Protein1 = y(1);
Protein2 = y(2);
Capsid = y(3);
Particle = y(4);
Virus = y(5);

dProtein1 = ks1*(vDNAtot/(Km1+vDNAtot)) - kd1*Protein1 ...
    - ks_C*vDNAtot*Protein1;
 
dProtein2 = ks2*(vDNAtot/(Km2+vDNAtot)) - kd2*Protein2 ...
    - ks_P*Capsid*Protein2;
 
dCapsid = ks_C*vDNAtot*Protein1 - ks_P*Capsid*Protein2;
 
dParticle = ks_P*Capsid*Protein2 - kex*Particle;
 
dVirus = kex*(1/500)*Particle;

M = [dProtein1,dProtein2,dCapsid,dParticle,dVirus]';
end