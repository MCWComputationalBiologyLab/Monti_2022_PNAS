function [GMaxes] = DetermineGMaxes(mpar,vDNAin0s) %Output: SP1,SP2,Virus
% Initial conditions
Protein10 = 0; Protein20 = 0; Capsid0 = 0;
Particle0 = 0; Virus0 = 1E-15;
y0 = [Protein10,Protein20,Capsid0,Particle0,Virus0];

options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'InitialStep',1e-2,...
    'NonNegative',(1:5), 'MaxOrder',5, 'BDF','on', 'Stats','off');

%Place to store models after they are calculated
Protein_1_Models = cell(2,length(vDNAin0s));
Protein_2_Models = cell(2,length(vDNAin0s));
Capsid_Models = cell(2,length(vDNAin0s));
Particle_Models = cell(2,length(vDNAin0s));
Virus_Models = cell(2,length(vDNAin0s));

Tp5_1tdata = [2,8,24,48,72,96];
Virustdata = [48,72,96];

% Simulate and plot model results
for i = 1:length(vDNAin0s)
    ODE_FH = @(t,y) Model(t,y,mpar,vDNAin0s(i));
    sols1 = ode15s(ODE_FH,Tp5_1tdata,y0,options);
    sols2 = ode15s(ODE_FH,Virustdata,y0,options);
    y1 = deval(Tp5_1tdata,sols1);
    y2 = deval(Virustdata,sols2);
    
  % Assigning models
    tout = Tp5_1tdata;
    toutVirus = Virustdata;
    Protein1 = y1(1,:);
    Protein2 = y1(2,:);
    Capsid = y1(3,:);
    Particle = y1(4,:);
    Virus = y2(5,:);
    
    % Storing Time Points
    Protein_1_Models{1,i} = tout;
    Protein_2_Models{1,i} = tout;
    Capsid_Models{1,i} = tout;
    Particle_Models{1,i} = tout;
    Virus_Models{1,i} = toutVirus;
    
    % Storing Models
    Protein_1_Models{2,i} = Protein1;
    Protein_2_Models{2,i} = Protein2;
    Capsid_Models{2,i} = Capsid;
    Particle_Models{2,i} = Particle;
    Virus_Models{2,i} = Virus;
end

% Define maxima
Protein1Max = Protein_1_Models{2,length(vDNAin0s)}(end);
Protein2Max = Protein_2_Models{2,length(vDNAin0s)}(end);
CapsidMax = Capsid_Models{2,length(vDNAin0s)}(end);
ParticleMax = Particle_Models{2,length(vDNAin0s)}(end);
VirusMax = Virus_Models{2,length(vDNAin0s)}(end);

% Find the max of each model
GMax_SP1 = Protein1Max+CapsidMax+ParticleMax;
GMax_SP2 = Protein2Max+ParticleMax;
GMax_Virus = VirusMax;

GMaxes = [GMax_SP1,GMax_SP2,GMax_Virus];

end