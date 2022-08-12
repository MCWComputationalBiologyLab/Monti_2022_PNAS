%% Driver
clear; clc; clf; close all;
format long g;

% Load in parameters 
load('Model_mpars_avg.mat');
%mpars

%% Load in data
vDNAData = load('vDNAData.txt');
vDNAin0 = vDNAData(1,2:length(vDNAData(1,:)));
Data1 = load('Tp5_1Data.txt');
Tp5_1tdata = Data1(:,1); %xdata
Tp5_1Data = Data1(:,2:length(Data1(1,:))); %ydata
Data2 = load('Tp5_2Data.txt');
Tp5_2tdata = Data2(:,1); %xdata
Tp5_2Data = Data2(:,2:length(Data2(1,:))); %ydata
Data3 = load('VirusData.txt');
Virustdata = Data3(:,1); %xdata
VirusData = Data3(:,2:length(Data3(1,:))); %ydata

%% Simulate Normalized Data
% Load in SDs
Tp5_1SDs = load('Tp5_1SDs.txt');
Tp5_2SDs = load('Tp5_2SDs.txt');
VirusSDs = load('VirusSDs.txt');

% Place to store models after they are calculated
Protein_1_Models = cell(2,length(Tp5_1Data(1,:)));
Protein_2_Models = cell(2,length(Tp5_2Data(1,:)));
Capsid_Models = cell(2,length(Tp5_1Data(1,:)));
Particle_Models = cell(2,length(Tp5_1Data(1,:)));
Virus_Models = cell(2,length(VirusData(1,:)));
Virus_Models_Raw = cell(2,length(VirusData(1,:)));

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
    sols1 = ode15s(ODE_FH,tspan,y0,options);
    
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

% Find the max of each model
GMax_SP1 = Protein1Max+CapsidMax+ParticleMax;
GMax_SP2 = Protein2Max+ParticleMax;
GMax_Virus = VirusMax;

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


figure(1); set(figure(1),'Units','inches','Position',[0.5 0.5 5.5 3.75]);
figure(2); set(figure(2),'Units','inches','Position',[0.5 0.5 5.5 3.75]);
figure(3); set(figure(3),'Units','inches','Position',[0.5 0.5 5.5 3.75]);
colors = ['b','g','r'];
line_colors = ["-b","-g","-r"];
shapes = ['o','^','s'];
% Plot
for i = 1:length(vDNAin0)
    vDNAin0_curr = vDNAin0(i);
    
    % Isolating normalized data
    Tp5_1Data_curr = Tp5_1Data(:,i);
    Tp5_2Data_curr = Tp5_2Data(:,i);
    VirusData_curr = VirusData(:,i);
   
    % Tp5_1 figure
    figure(1);
    % Plot Tp5_1 model
    lc = line_colors(i);
    curr_line = plot(SumP1s{1,i},SumP1s{2,i}/GMax_SP1,lc,'LineWidth',1.5);
    curr_line.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on; 
    %Plot data
    c = colors(i);
    errorbar(Tp5_1tdata,Tp5_1Data_curr,Tp5_1SDs(:,i),shapes(i),...
        'MarkerSize',8,'MarkerFaceColor',c,'MarkerEdgeColor',c,...
        'LineWidth',1.5,'Color',c,'CapSize',6); hold on;
    
    % Tp5_2 Figure
    figure(2);
    % Plot Tp5_2 model
    lc = line_colors(i);
    curr_line = plot(SumP2s{1,i},SumP2s{2,i}/GMax_SP2,lc,'LineWidth',1.5); 
    curr_line.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on; 
    %Plot data
    c = colors(i);
    errorbar(Tp5_2tdata,Tp5_2Data_curr,Tp5_2SDs(:,i),shapes(i),...
        'MarkerSize',8,'MarkerFaceColor',c,'MarkerEdgeColor',c,...
        'LineWidth',1.5,'Color',c,'CapSize',6); hold on;
    
    % Virus Figure
    figure(3);
    % Plot Virus model
    lc = line_colors(i);
    curr_line = semilogy(Virus_Models{1,i},Virus_Models{2,i},lc,...
        'LineWidth',1.5);
    curr_line.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on; 
    %Plot data
    c = colors(i);
    errorbar(Virustdata,VirusData_curr,VirusSDs(:,i),shapes(i),...
        'MarkerSize',8,'MarkerFaceColor',c,'MarkerEdgeColor',c,...
        'LineWidth',1.5,'Color',c,'CapSize',6); hold on;
end

figure(1);
set(gca,'FontName','Arial','FontSize',16,'XTick',(0:24:96),'LineWidth',1)
% xlabel('Time (hpi)','fontsize',36,'FontName','Arial');
% ylabel('Tp5_{1} (A.U.)','fontsize',36,'FontName','Arial');
% legend('vDNA_{in,0} = 3','vDNA_{in,0} = 14',...
%     'vDNA_{in,0} = 131','Location','NorthWest'); box off; legend boxoff
% [~,icons] = legend('vDNA_{in,0} = 3','vDNA_{in,0} = 14',...
%     'vDNA_{in,0} = 131','Location','NorthWest'); box off; legend boxoff
axis([0 96 0 1.2]);

figure(2);
set(gca,'FontName','Arial','FontSize',16,'XTick',(0:24:96),'LineWidth',1)
% xlabel('Time (hpi)','fontsize',36,'FontName','Arial');
% ylabel('Tp5_{2} (A.U.)','fontsize',36,'FontName','Arial');
% legend('vDNA_{in,0} = 3','vDNA_{in,0} = 14',...
%     'vDNA_{in,0} = 131','Location','NorthWest'); box off; legend boxoff
% [~,icons] = legend('vDNA_{in,0} = 3','vDNA_{in,0} = 14',...
%     'vDNA_{in,0} = 131','Location','NorthWest'); box off; legend boxoff
axis([0 96 0 1.2]);

figure(3);
set(gca,'FontName','Arial','FontSize',16,'XTick',(0:24:96),'LineWidth',1)
% xlabel('Time (hpi)','fontsize',36,'FontName','Arial');
% ylabel('Virus (A.U.)','fontsize',36,'FontName','Arial');
% legend('vDNA_{in,0} = 3','vDNA_{in,0} = 14',...
%     'vDNA_{in,0} = 131','Location','NorthWest'); box off; legend boxoff
% [~,icons] = legend('vDNA_{in,0} = 3','vDNA_{in,0} = 14',...
%     'vDNA_{in,0} = 131','Location','SouthEast'); box off; legend boxoff
axis([0 96 10^-15 10^1]);

for i = 1:3
    figure(i);
    set(gcf,'color','white');
    box off;
end