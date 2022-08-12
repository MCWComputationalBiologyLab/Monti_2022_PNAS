%% Driver_3D_Plot
clear; clc; clf; close all;
format long g;

% Load in parameters 
load('Model_mpars_avg.mat'); %Change here
vDNAData = load('vDNAData.txt');
%mpars

%% Simulate full range of MOIs
logvDNAin0v = [-4:0.25:4];
vDNAin0 = 10.^logvDNAin0v;

% Initial conditions
Tp5_10 = 0; Tp5_20 = 0; Capsid0 = 0; Particle0 = 0; Virus0 = 1E-15;
y0 = [Tp5_10,Tp5_20,Capsid0,Particle0,Virus0]; %IC

tspan = [0:1:96];
options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'InitialStep',1e-2,...
    'NonNegative',(1:5), 'MaxOrder',5, 'BDF','on', 'Stats','off');

% Place to store models after they are calculated
Protein_1_Models = cell(2,length(vDNAin0));
Protein_2_Models = cell(2,length(vDNAin0));
Capsid_Models = cell(2,length(vDNAin0));
Particle_Models = cell(2,length(vDNAin0));
Virus_Models = cell(2,length(vDNAin0));

% Simulate and plot model results
for i = 1:length(vDNAin0)
    % Calculate model solution for each vDNAin0
    ODE_FH = @(t,y) Model(t,y,mpars,vDNAin0(i));
    sols1 = ode15s(ODE_FH,tspan,y0,options);
    y = deval(sols1,tspan);
    
    % Assigning models
    tout = tspan;
    Protein1 = y(1,:);
    Protein2 = y(2,:);
    Capsid = y(3,:);
    Particle = y(4,:);
    Virus = y(5,:);
    
    % Storing Time Points
    Protein_1_Models{1,i} = tout;
    Protein_2_Models{1,i} = tout;
    Capsid_Models{1,i} = tout;
    Particle_Models{1,i} = tout;
    Virus_Models{1,i} = tout;
    
    % Storing Models
    Protein_1_Models{2,i} = Protein1;
    Protein_2_Models{2,i} = Protein2;
    Capsid_Models{2,i} = Capsid;
    Particle_Models{2,i} = Particle;
    Virus_Models{2,i} = Virus;
end

% Define maxima
[GMax_SP1,GMax_SP2,GMax_Virus,na] = parEstSims(mpars,vDNAData(1,2:...
    length(vDNAData(1,:))));

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

%% Arrange Data
Tp5_1Matrix = [];
for i = 1:length(vDNAin0)
   Tp5_1Matrix(:,i) = SumP1s{2,i}; 
end

Tp5_2Matrix = [];
for i = 1:length(vDNAin0)
   Tp5_2Matrix(:,i) = SumP2s{2,i}; 
end

VirusMatrix = [];
for i = 1:length(vDNAin0)
   VirusMatrix(:,i) = Virus_Models{2,i}; 
end

Tp5_1Matrix = Tp5_1Matrix';
Tp5_2Matrix = Tp5_2Matrix';
VirusMatrix = VirusMatrix';

%% Tp5_1 3D (Figure 1)
figure(1);
mesh(tspan,vDNAin0,Tp5_1Matrix/GMax_SP1,'LineWidth',1.5); hold on;
colormap jet;
axis([0 96 10^-4 10^4 0 1.2]); box on;
set(gca,'XTick',(0:24:96),'YTick',[10^-4 10^0 10^4],'ZTick',(0:0.5:1.2),...
    'FontSize',16,'LineWidth',1);
set(gca,'YScale','log');
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 5.5 3.5]);
set(gcf,'Units','inches','Position',[0.5 0.5 5.5 3.5]);
set(gcf,'color','white');

%% Tp5_2 3D (Figure 2)
figure(2);
mesh(tspan,vDNAin0,Tp5_2Matrix/GMax_SP2,'LineWidth',1.5); hold on;
colormap jet;
axis([0 96 10^-4 10^4 0 1.2]); box on;
set(gca,'XTick',(0:24:96),'YTick',[10^-4 10^0 10^4],'ZTick',(0:0.5:1.2),...
    'FontSize',16,'LineWidth',1);
set(gca,'YScale','log');
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 5.5 3.5]);
set(gcf,'Units','inches','Position',[0.5 0.5 5.5 3.5]);
set(gcf,'color','white');

%% Un-normalized virus 3D (Figure 3)
foldChanges = [];

Virus24s = [];
fact = 1E7;
for i = 1:length(vDNAin0)
    curr_vDNAin0 = vDNAin0(i);
    %Conv. vDNAin0 to MOI then from MOI (IU/cell) multiply by number of
    %cells (1E6) then divide by volume of viral stock used (500 uL) then
    %convert to mL (1E3 uL/mL)
    temp = Virus_Models{2,i};
    curr_Virus = temp(length(tspan))*fact;
    Virus24 = temp(25)*fact;
    foldChanges(i,1) = curr_vDNAin0;
    foldChanges(i,2) = curr_Virus/Virus24;
end

colorMap = [];
for i = 1:length(foldChanges)
    for j = 1:length(tspan)
        colorMap(i,j) = foldChanges(i,2);
    end
end

figure(3);
mesh(tspan,vDNAin0,VirusMatrix*fact,colorMap,'LineWidth',1.5); hold on;
colormap jet;
axis([0 96 10^-4 10^4 10^-8 10^8]); box on;
set(gca,'XTick',(0:24:96),'YTick',[10^-4,10^0,10^4],'ZTick',10.^...
    (-8:5:8));
set(gca,'LineWidth',1,'FontSize',16);
set(gca,'YScale','log'); set(gca,'ZScale','log','ZMinorTick','on');
set(gca,'ColorScale','log');
set(gca,'MinorGridLineStyle','-','MinorGridColor',[0.15 0.15 0.15],...
    'MinorGridAlpha',0.15);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 5.5 3.5]);
set(gcf,'Units','inches','Position',[0.5 0.5 5.5 3.5]);
set(gcf,'color','white');

foldChanges = [];

%% Normalized virus 3D (Figure 4)
Virus24s = [];
fact = 1;
for i = 1:length(vDNAin0)
    curr_vDNAin0 = vDNAin0(i);
    %Conv. vDNAin0 to MOI then from MOI (IU/cell) multiply by number of
    %cells (1E6) then divide by volume of viral stock used (500 uL) then
    %convert to mL (1E3 uL/mL)
    temp = Virus_Models{2,i};
    curr_Virus = temp(length(tspan))*fact;
    Virus24 = temp(25)*fact;
    foldChanges(i,1) = curr_vDNAin0;
    foldChanges(i,2) = curr_Virus/Virus24;
end

colorMap = [];
for i = 1:length(foldChanges)
    for j = 1:length(tspan)
        colorMap(i,j) = foldChanges(i,2);
    end
end

figure(4);
mesh(tspan,vDNAin0,VirusMatrix*fact,'LineWidth',1.5); hold on;
colormap jet;
axis([0 96 10^-4 10^4 10^-15 10^1]); box on;
set(gca,'XTick',(0:24:96),'YTick',[10^-4,10^0,10^4],'ZTick',...
    10.^(-15:5:0));
set(gca,'LineWidth',1,'FontSize',16);
set(gca,'YScale','log'); set(gca,'ZScale','log','ZMinorTick','on');
set(gca,'ColorScale','log');
set(gca,'MinorGridLineStyle','-','MinorGridColor',[0.15 0.15 0.15],...
    'MinorGridAlpha',0.15);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 5.5 3.5]);
set(gcf,'Units','inches','Position',[0.5 0.5 5.5 3.5]);
set(gcf,'color','white');

%% Figures originally from Fig. 4 in manuscript

% Calc vDNAtot
vDNAtot_Models = cell(2,length(vDNAin0));

for i = 1:length(vDNAin0)
    vDNAin0_curr = vDNAin0(i);
    vDNAtot_curr = Calc_vDNA(tspan,vDNAin0_curr);
    vDNAtot_Models{2,i} = vDNAtot_curr;
end

% Arrange Data
vDNAtotMatrix = [];
for i = 1:length(vDNAin0)
   vDNAtotMatrix(:,i) = vDNAtot_Models{2,i}; 
end
vDNAtotMatrix = vDNAtotMatrix';

VirusMatrix = [];
fact = 1E7;
for i = 1:length(vDNAin0)
   VirusMatrix(:,i) = Virus_Models{2,i}*fact; 
end
VirusMatrix = VirusMatrix';

Tp5_1Matrix = [];
for i = 1:length(vDNAin0)
   Tp5_1Matrix(:,i) = SumP1s{2,i}; 
end
Tp5_1Matrix = Tp5_1Matrix';

Tp5_2Matrix = [];
for i = 1:length(vDNAin0)
   Tp5_2Matrix(:,i) = SumP2s{2,i}; 
end
Tp5_2Matrix = Tp5_2Matrix';

tMatrix = ones(length(vDNAin0),length(tspan));
for i = 1:length(vDNAin0)
    tMatrix(i,:) = tspan.*tMatrix(i,:);
end

% Plot
for i = 1:length(vDNAin0)
    curr_vDNAin0 = vDNAin0(i);
    curr_vDNAEnd = vDNAtotMatrix(i,length(tspan));
    foldChanges(i,1) = curr_vDNAin0;
    foldChanges(i,2) = curr_vDNAEnd/curr_vDNAin0;
end

colorMap = [];
for i = 1:length(foldChanges)
    for j = 1:length(tspan)
        colorMap(i,j) = foldChanges(i,2);
    end
end

figure(5); % Tp5_1 vs vDNAtot and t (New) (Figure 5)
mesh(tMatrix,vDNAtotMatrix,Tp5_1Matrix/GMax_SP1,colorMap,'LineWidth',1.5);
hold on;
colormap jet;
set(gca,'ColorScale','log');
axis([0 96 10^-4 10^4 0 1.2]); box on;
set(gca,'XTick',(0:24:96),'YTick',[10^-4 10^0 10^4],'FontSize',16,...
    'LineWidth',1,'FontName','Arial'); 
set(gca,'YScale','log','MinorGridLineStyle','-',...
    'MinorGridColor',[0.15 0.15 0.15],'MinorGridAlpha',0.15);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 5.5 3.5]);
set(gcf,'Units','inches','Position',[0.5 0.5 5.5 3.5]);
set(gcf,'color','white');

figure(6); % Tp5_2 vs vDNAtot and t (Originally Fig 4C Left) (Figure 6)
mesh(tMatrix,vDNAtotMatrix,Tp5_2Matrix/GMax_SP2,colorMap,'LineWidth',1.5);
hold on;
colormap jet;
set(gca,'ColorScale','log');
axis([0 96 10^-4 10^4 0 1.2]); box on;
set(gca,'XTick',(0:24:96),'YTick',[10^-4 10^0 10^4],'FontSize',16,...
    'LineWidth',1,'FontName','Arial'); 
set(gca,'YScale','log','MinorGridLineStyle','-',...
    'MinorGridColor',[0.15 0.15 0.15],'MinorGridAlpha',0.15);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 5.5 3.5]);
set(gcf,'Units','inches','Position',[0.5 0.5 5.5 3.5]);
set(gcf,'color','white');

figure(7); % Virus vs vDNAtot and t (Originally Fig 4C Right) (Figure 7)
mesh(tMatrix,vDNAtotMatrix,VirusMatrix/fact,colorMap,'LineWidth',1.5); 
hold on;
colormap jet;
set(gca,'ColorScale','log');
axis([0 96 10^-4 10^4 10^-15 10^1]); box on; 
set(gca,'XTick',(0:24:96),'YTick',[10^-4 10^-0 10^4],'ZTick',...
    [10^-15 10^-7 10^1],'FontSize',16,'LineWidth',1,'FontName','Arial'); 
set(gca,'YScale','log','ZScale','log','MinorGridLineStyle','-',...
    'MinorGridColor',[0.15 0.15 0.15],'MinorGridAlpha',0.15);
set(gcf,'Units','inches','PaperPosition',[0.5 0.5 5.5 3.5]);
set(gcf,'Units','inches','Position',[0.5 0.5 5.5 3.5]);
set(gcf,'color','white');
