%% Figure 4 Other Images
clear; clc; clf; close all;
format long g;

load('Model_mpars_avg.mat'); 
%mpars

%% Plot Figure 4E Left (Figure 1)
vDNAin0 = [1:4:400];

% Initial conditions
Tp5_10 = 0; Tp5_20 = 0; Capsid0 = 0; Particle0 = 0; Virus0 = 1E-15;
y0 = [Tp5_10,Tp5_20,Capsid0,Particle0,Virus0]; %IC

% Place to store models after they are calculated
%Virus_Models0 = [];
Virus_Models24 = [];
Virus_Models96 = [];

tspan = [0,24,96];
fact = 1E7;

options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'InitialStep',1e-2,...
    'NonNegative',(1:5), 'MaxOrder',5, 'BDF','on', 'Stats','off');

% Simulate and plot model results
for i = 1:length(vDNAin0)
    % Calculate model solution for each vDNAin0 at 96 hpi
    ODE_FH = @(t,y) Model(t,y,mpars,vDNAin0(i));
    sols1 = ode15s(ODE_FH,tspan,y0,options);
    y = deval(sols1,tspan);
    Virus96 = y(5,3)*fact;
    Virus24 = y(5,2)*fact;

    % Storing Models
    Virus_Models24(i) = Virus24;
    Virus_Models96(i) = Virus96;
end

% Find dg/dvDNAin0
g = Virus_Models96./Virus_Models24;
dg = gradient(g);

figure(1);
plot(vDNAin0,g,'-k','LineWidth',1.5);
set(gcf,'color','white');
xlim([0 100]);
set(gca,'FontName','Arial','FontSize',16,'LineWidth',1,'FontName','Arial');
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 5.5 3.75]);
set(gcf,'Units','inches','Position',[0.5 0.5 5.5 3.75]);
box off;

% Find change from + to - (indicating a max) in dg
maxIdx = 1;
idx = 1;
plusMinusChecker = false;

% Logic: if dg at the current index is + and dg at the next index is - 
% stop the loop
while(~plusMinusChecker)
    if(dg(idx) > 0 && dg(idx+1) < 0)
        maxIdx = idx
        plusMinusChecker = true;
    end
    idx = idx + 1;
end

lb = vDNAin0(maxIdx); %Last positive
ub = vDNAin0(maxIdx+1); %First negative

maxvDNAin0 = [lb,ub]

lb = lb/26.17 * 1E6 * (1/500) * (1/(1E3));
ub = ub/26.17 * 1E6 * (1/500) * (1/(1E3));

maxMOI = [lb,ub]

%% Plot Figure 4E Right (Figure 2)

% Place to store models after they are calculated
vDNA_Models0 = vDNAin0;
vDNA_Models96 = [];

tspan = [0,96];

% Simulate and plot model results
for i = 1:length(vDNAin0)
    vDNAin0_curr = vDNAin0(i);

    % Calculate and store model solution for each vDNAin0 at 96 hpi
    vDNA_Models96(i) = Calc_vDNA(96,vDNAin0_curr);
end

% Find dg/dvDNAin0
g = vDNA_Models96./vDNA_Models0;
dg = gradient(g);

figure(2);
plot(vDNAin0,g,'-k','LineWidth',1.5);
set(gcf,'color','white');
xlim([0 100]);
set(gca,'FontName','Arial','FontSize',16,'LineWidth',1,'FontName','Arial');
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 5.5 3.75]);
set(gcf,'Units','inches','Position',[0.5 0.5 5.5 3.75]);
box off;

% Find change from + to - (indicating a max) in dg
maxIdx = 1;
idx = 1;
plusMinusChecker = false;

% Logic: if dg at the current index is + and dg at the next index is - 
% stop the loop
while(~plusMinusChecker)
    if(dg(idx) > 0 && dg(idx+1) < 0)
        maxIdx = idx
        plusMinusChecker = true;
    end
    idx = idx + 1;
end

lb = vDNAin0(maxIdx); %Last positive
ub = vDNAin0(maxIdx+1); %First negative

maxvDNAin0 = [lb,ub]

lb = lb/26.17 * 1E6 * (1/500) * (1/(1E3));
ub = ub/26.17 * 1E6 * (1/500) * (1/(1E3));

maxMOI = [lb,ub]


%% Plot IncuCyte Figures (Figure 3)
vDNAData = load('vDNAData.txt');
vDNAin0 = vDNAData(1,2:length(vDNAData(1,:)));
Data4 = load('IC Data_2.txt');
ICtdata = Data4(:,1);
ICData = Data4(:,2:length(Data4(1,:))); %ydata

% Place to store models after they are calculated
Protein_1_Models = cell(2,length(ICData(1,:)));
Protein_2_Models = cell(2,length(ICData(1,:)));
Capsid_Models = cell(2,length(ICData(1,:)));
Particle_Models = cell(2,length(ICData(1,:)));
Virus_Models = cell(2,length(ICData(1,:)));
Virus_Models_Raw = cell(2,length(ICData(1,:)));

% Initial conditions
Protein10 = 0; Protein20 = 0; Capsid0 = 0; Particle0 = 0; Virus0 = 1E-15;
y0 = [Protein10,Protein20,Capsid0,Particle0,Virus0];

tspan = [0:1:96];
options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'InitialStep',1e-2,...
    'NonNegative',(1:5), 'MaxOrder',5, 'BDF','on', 'Stats','off');

% Simulate and plot model results
for i = 1:length(vDNAin0)
    % Calculate model solution for each vDNAin0
    ODE_FH = @(t,y) Model(t,y,mpars,vDNAin0(i));
    sols1 = ode15s(ODE_FH,tspan,y0,[]);
    
    % Assigning models
    tout = sols1.x;
    Protein1 = sols1.y(1,:);
    Protein2 = sols1.y(2,:);
    Capsid = sols1.y(3,:);
    Particle = sols1.y(4,:);
    Virus = sols1.y(5,:);
    
    % Storing Time Points
    Protein_1_Models{1,i} = tout;
    Protein_2_Models{1,i} = tout;
    Capsid_Models{1,i} = tout;
    Particle_Models{1,i} = tout;
    Virus_Models{1,i} = tout;
    Virus_Models_Raw{1,i} = tout;
    
    % Storing Models
    Protein_1_Models{2,i} = Protein1;
    Protein_2_Models{2,i} = Protein2;
    Capsid_Models{2,i} = Capsid;
    Particle_Models{2,i} = Particle;
    Virus_Models{2,i} = Virus;
    Virus_Models_Raw{2,i} = Virus;
end

% Define maxima
Protein1Max = Protein_1_Models{2,length(vDNAin0)}(end);
Protein2Max = Protein_2_Models{2,length(vDNAin0)}(end);
CapsidMax = Capsid_Models{2,length(vDNAin0)}(end);
ParticleMax = Particle_Models{2,length(vDNAin0)}(end);
VirusMax = Virus_Models{2,length(vDNAin0)}(end);

SumP1s = cell(2,length(vDNAin0));
SumP2s = cell(2,length(vDNAin0));

for i = 1:length(vDNAin0)
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

figure(3);
colors = ['b','g','r'];
line_colors = ["-b","-g","-r"];
shapes = ['o','^','s'];
% Plot
for i = 1:length(vDNAin0)
    vDNAin0_curr = vDNAin0(i);
    
    % Isolating normalized data
    ICData_curr = ICData(:,i);
    
    % Tp5_2 Figure
    figure(3);
    % Plot Tp5_2 model
    lc = line_colors(i);
    curr_line = plot(SumP2s{1,i},SumP2s{2,i}/GMax_SP2,lc,'LineWidth',2); 
    curr_line.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on; 
end

figure(3);
set(gca,'FontName','Arial','FontSize',16,'XTick',(0:24:96),'LineWidth',1)
axis([0 96 0 1.2]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 5.5 3.75]);
set(gcf,'Units','inches','Position',[0.5 0.5 5.5 3.75]);
box off
set(gcf,'color','white');

figure(4);
colors = ['b','g','r'];
line_colors = ["-b","-g","-r"];
shapes = ['o','^','s'];
% Plot
for i = 1:length(vDNAin0)
    vDNAin0_curr = vDNAin0(i);
    
    % Isolating normalized data
    ICData_curr = ICData(:,i);
    
    % Tp5_2 Figure
    figure(4);
    % Plot Tp5_2 model
    c = colors(i);
    plot(ICtdata,ICData_curr,shapes(i),'MarkerSize',8,'MarkerFaceColor',...
        c,'MarkerEdgeColor',c,'LineWidth',1.5,'Color',c); hold on;
end

figure(4);
set(gca,'FontName','Arial','FontSize',16,'XTick',(0:24:96),'LineWidth',1)
axis([0 96 0 1.2]);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 5.5 3.75]);
set(gcf,'Units','inches','Position',[0.5 0.5 5.5 3.75]);
box off
set(gcf,'color','white');