function Err = Calc_SSE(Model,Data1,Data2,Data3,vDNAin0s,Maxes)
Tp5_1tdata = Data1(:,1); %xdata
Tp5_1Data = Data1(:,2:length(Data1(1,:))); %ydata
Tp5_2tdata = Data2(:,1); %xdata
Tp5_2Data = Data2(:,2:length(Data2(1,:))); %ydata
Virustdata = Data3(:,1); %xdata
VirusData = Data3(:,2:length(Data3(1,:))); %ydata

options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'InitialStep',1e-2,...
    'NonNegative',(1:5), 'MaxOrder',5, 'BDF','on', 'Stats','off');

% Initial conditions
Protein10 = 0; Protein20 = 0; Capsid0 = 0;
Particle0 = 0; Virus0 = 1e-20;
y0 = [Protein10,Protein20,Capsid0,Particle0,Virus0];

%Place to store models after they are calculated
Protein_1_Models = cell(2,length(Tp5_1Data(1,:)));
Protein_2_Models = cell(2,length(Tp5_2Data(1,:)));
Capsid_Models = cell(2,length(Tp5_1Data(1,:)));
Particle_Models = cell(2,length(Tp5_1Data(1,:)));
Virus_Models = cell(2,length(VirusData(1,:)));

% Simulate and plot model results
for i = 1:length(vDNAin0s)
    % Calculate model solution for each vDNAin0
    ODE_FH = @(t,y) Model(t,y,vDNAin0s(i));
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
    Virus =y2(5,:);
    
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

% Incorporating capsid and particle
Protein1Max = Maxes(1);
Protein2Max = Maxes(2);
CapsidMax = Maxes(3);
ParticleMax = Maxes(4);
VirusMax = Maxes(5);

SumP1s = cell(2,length(vDNAin0s));
SumP2s = cell(2,length(vDNAin0s));

for i = 1:length(vDNAin0s)
    SumP1s{1,i} = Protein_1_Models{1,i};
    currP1 = Protein_1_Models{2,i};
    SumP2s{1,i} = Protein_2_Models{1,i};
    currP2 = Protein_2_Models{2,i};
    
    currCapsid = Capsid_Models{2,i};
    currParticle = Particle_Models{2,i};
    
    SumP1s{2,i} = currP1 + currCapsid + currParticle;
    SumP2s{2,i} = currP2 + currParticle;
end

% Find the max of each model
GMax_SP1 = Protein1Max+CapsidMax+ParticleMax;
GMax_SP2 = Protein2Max+ParticleMax;
GMax_Virus = VirusMax;

Err = 0;

%% Calculate SSE: Data norm to its GMax and model norm to its GMax
for l = 1:length(vDNAin0s)
    % Isolating normalized data (Data/GMax_Data)
    Tp5_1Data_curr = Tp5_1Data(:,l)';
    Tp5_2Data_curr = Tp5_2Data(:,l)';
    VirusData_curr = VirusData(:,l)';
    
    % Isolating normalized model (Model/GMax_Model)
    Protein_1_Model_curr = SumP1s{2,l}/GMax_SP1;
    Protein_2_Model_curr = SumP2s{2,l}/GMax_SP2;
    VirusModel_curr = Virus_Models{2,l}/GMax_Virus;
    
    % Calculate SSE
    Err_curr_Tp5_1 = sum((Tp5_1Data_curr - Protein_1_Model_curr).^2);
    
    Err_curr_Tp5_2 = sum((Tp5_2Data_curr - Protein_2_Model_curr).^2);
    
    Err_curr_Virus = sum((VirusData_curr - VirusModel_curr).^2);
    
    % Calculate total error
    Err = Err + Err_curr_Tp5_1 + Err_curr_Tp5_2 + Err_curr_Virus;
end

end