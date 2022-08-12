function [c,ceq] = mycon(x,~,~,~,~)

Protein1Max = 5000;
Protein2Max = 2500;
CapsidMax = 500;
ParticleMax = 25;
VirusMax = 1;

vDNAmax = 130.7684499;
tmax = 96;

% Initial conditions
Protein10 = 0; Protein20 = 0; Capsid0 = 0;
Particle0 = 0; Virus0 = 1E-15;
y0 = [Protein10,Protein20,Capsid0,Particle0,Virus0];

tspan = [0:1:96];
ODE_FH = @(t,y) Model(t,y,x,vDNAmax);
sols1 = ode15s(ODE_FH,tspan,y0,[]);
y = deval(tmax,sols1);

Protein1 = y(1,:);
Protein2 = y(2,:);
Capsid = y(3,:);
Particle = y(4,:);
Virus = y(5,:);

c = 0;

ceq = [Protein1-Protein1Max,Protein2-Protein2Max,Capsid-CapsidMax,...
    Particle-ParticleMax,Virus-VirusMax];
end