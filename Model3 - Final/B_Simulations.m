%% Model Simulations
clear; clc; clf; close all;
format long g;

% Load in parameters 
load('Model_mpars_avg.mat'); %Change here
%mpars

%% Simulate Normalized Data
% Load in data
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
Protein10 = 0; Protein20 = 0; Capsid0 = 0;
Particle0 = 0; Virus0 = 1E-15;
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

GMaxes = DetermineGMaxes(mpars,vDNAin0);

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


figure(1); set(figure(1),'Units','inches','Position',[0.5 0.5 16 4])
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
    subplot(1,3,1);
    % Plot Tp5_1 model
    lc = line_colors(i);
    curr_line = plot(SumP1s{1,i},SumP1s{2,i}/GMax_SP1,lc,'LineWidth',2);
    curr_line.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on; 
    %Plot data
    c = colors(i);
    errorbar(Tp5_1tdata,Tp5_1Data_curr,Tp5_1SDs(:,i),shapes(i),...
        'MarkerSize',8,'MarkerFaceColor',c,'MarkerEdgeColor',c,...
        'LineWidth',1.5,'Color',c,'CapSize',6); hold on;
    
    % Tp5_2 Figure
    subplot(1,3,2);
    % Plot Tp5_2 model
    lc = line_colors(i);
    curr_line = plot(SumP2s{1,i},SumP2s{2,i}/GMax_SP2,lc,'LineWidth',2); 
    curr_line.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on; 
    %Plot data
    c = colors(i);
    errorbar(Tp5_2tdata,Tp5_2Data_curr,Tp5_2SDs(:,i),shapes(i),...
        'MarkerSize',8,'MarkerFaceColor',c,'MarkerEdgeColor',c,...
        'LineWidth',1.5,'Color',c,'CapSize',6); hold on;
    
    % Virus Figure
    subplot(1,3,3);
    % Plot Virus model
    lc = line_colors(i);
    curr_line = semilogy(Virus_Models{1,i},Virus_Models{2,i}/GMax_Virus,...
        lc,'LineWidth',2);
    curr_line.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hold on; 
    %Plot data
    c = colors(i);
    errorbar(Virustdata,VirusData_curr,VirusSDs(:,i)/7600000,shapes(i),...
        'MarkerSize',8,'MarkerFaceColor',c,'MarkerEdgeColor',c,...
        'LineWidth',1.5,'Color',c,'CapSize',6); hold on;
end

subplot(1,3,1);
set(gca,'FontName','Arial','FontSize',16,'XTick',(0:24:96),'LineWidth',1);
axis([0 96 0 1.2]);
set(gcf,'color','white');
xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
ylabel('Tp5_{1,tot} (A.U.)','fontsize',18,'FontName','Arial');
legend('vDNA_{in,0} = 3','vDNA_{in,0} = 14','vDNA_{in,0} = 131',...
    'Location','NorthWest'); box off;
legend boxoff;

subplot(1,3,2);
set(gca,'FontName','Arial','FontSize',16,'XTick',(0:24:96),'LineWidth',1);
axis([0 96 0 1.2]);
set(gcf,'color','white');
xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
ylabel('Tp5_{2,tot} (A.U.)','fontsize',18,'FontName','Arial');
legend('vDNA_{in,0} = 3','vDNA_{in,0} = 14','vDNA_{in,0} = 131',...
    'Location','NorthWest'); box off;
legend boxoff;

subplot(1,3,3);
set(gca,'FontName','Arial','FontSize',16,'XTick',(0:24:96),'LineWidth',1);
axis([0 96 10^-15 10^1]);
set(gcf,'color','white');
xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
ylabel('Virus (A.U.)','fontsize',18,'FontName','Arial');
legend('vDNA_{in,0} = 3','vDNA_{in,0} = 14','vDNA_{in,0} = 131',...
    'Location','SouthEast'); box off;
legend boxoff;

%% Simulate unnormalized model
% Load vDNAin0s
logvDNAin0v = [-4:0.5:4];
vDNAin0 = 10.^logvDNAin0v;
%vDNAin0 = [0.3,3,14,131,300,1000];

% Initial conditions
Protein10 = 0; Protein20 = 0; Capsid0 = 0;
Particle0 = 0; Virus0 = 1E-15;
y0 = [Protein10,Protein20,Capsid0,Particle0,Virus0];

tspan = [0:1:96];
options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'InitialStep',1e-2,...
    'NonNegative',(1:5), 'MaxOrder',5, 'BDF','on', 'Stats','off');

% Simulate
for i = 1:length(vDNAin0)
    vDNAin0_curr = vDNAin0(i);
    
    % Calculate model values for each MOI
    ODE_FH = @(t,y) Model(t,y,mpars,vDNAin0_curr);
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
    
    % Storing Models
    Protein_1_Models{2,i} = Protein1;
    Protein_2_Models{2,i} = Protein2;
    Capsid_Models{2,i} = Capsid;
    Particle_Models{2,i} = Particle;
    Virus_Models{2,i} = Virus;
end

% Maxima defined previously (see lines 82-85)

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

% Plot
figure(2); set(figure(2),'Units','inches','Position',[0.5 0.5 16 8])
for i = 1:length(vDNAin0)
    vDNAin0_curr = vDNAin0(i); 
    
    % Tp5_1 Figure
    subplot(2,4,1);
    % Plot Tp5_1
    plot(Protein_1_Models{1,i},Protein_1_Models{2,i},...
        '-','LineWidth',2); 
    hold on;
    
     % Tp5_2 Figure
    subplot(2,4,2);
    % Plot Tp5_1
    plot(Protein_2_Models{1,i},Protein_2_Models{2,i},...
        '-','LineWidth',2); 
    hold on;
    
    % Tp5_1_tot Figure
    subplot(2,4,3);
    % Plot Tp5_1
    plot(SumP1s{1,i},SumP1s{2,i},'-','LineWidth',2); 
    hold on;

    % Tp5_2_tot Figure
    subplot(2,4,4);
    % Plot Tp5_2
    plot(SumP2s{1,i},SumP2s{2,i},'-','LineWidth',2); 
    hold on; 

    % Capsid Figure
    subplot(2,4,5);
    % Plot Capsid
    plot(Capsid_Models{1,i},Capsid_Models{2,i},...
        '-','LineWidth',2); 
    hold on; 

    % Particle Figure
    subplot(2,4,6);
    % Plot Capsid
    plot(Particle_Models{1,i},Particle_Models{2,i},...
        '-','LineWidth',2); 
    hold on;

    % Virus Figures
    subplot(2,4,7);
    % Plot Virus
    semilogy(Virus_Models{1,i},Virus_Models{2,i},'-',...
        'LineWidth',2); 

    hold on; 

    %Un-normalized Virus Figures
    fact = 1E7;
    subplot(2,4,8);
    % Plot Virus
    semilogy(Virus_Models{1,i},Virus_Models{2,i}*fact,'-',...
        'LineWidth',2); 
    hold on;
end

subplot(2,4,1)
xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
ylabel('Tp5_{1,free} (genomes/cell)','fontsize',18,'FontName','Arial');
axis([0 96 0 5000]);
xticks([0:24:96]);
set(gca,'fontsize',16,'FontName','Arial','LineWidth',1);
box off;

subplot(2,4,2)
xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
ylabel('Tp5_{2,free} (genomes/cell)','fontsize',18,'FontName','Arial');
axis([0 96 0 3000]);
xticks([0:24:96]);
set(gca,'fontsize',16,'FontName','Arial');
box off;

subplot(2,4,3)
xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
ylabel('Tp5_{1,tot} (genomes/cell)','fontsize',18,'FontName','Arial');
axis([0 96 0 5000]);
xticks([0:24:96]);
set(gca,'fontsize',16,'FontName','Arial');
box off;

subplot(2,4,4)
xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
ylabel('Tp5_{2,tot} (genomes/cell)','fontsize',18,'FontName','Arial');
axis([0 96 0 3000]);
xticks([0:24:96]);
set(gca,'fontsize',16,'FontName','Arial');
box off;

subplot(2,4,5)
xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
ylabel('Capsid (genomes/cell)','fontsize',18,'FontName','Arial');
axis([0 96 0 600]);
xticks([0:24:96]);
set(gca,'fontsize',16,'FontName','Arial');
box off;

subplot(2,4,6)
xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
ylabel('Particle (genomes/cell)','fontsize',18,'FontName','Arial');
axis([0 96 0 30]);
xticks([0:24:96]);
set(gca,'fontsize',16,'FontName','Arial');
box off;

subplot(2,4,7)
xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
ylabel('Virus (genomes/cell)','fontsize',18,'FontName','Arial');
axis([0 96 0 1]);
xticks([0:24:96]);
set(gca,'fontsize',16,'FontName','Arial');
box off;

subplot(2,4,8)
xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
ylabel('Virus (IU/mL)','fontsize',18,'FontName','Arial');
axis([0 96 10^-8 10^7]);
yticks([10^-8 10^-3 10^2 10^7]);
xticks([0:24:96]);
set(gca,'fontsize',16,'FontName','Arial');
box off;

figure(2);
set(gcf,'color','white');

%% Simulate normalized model
% Load vDNAin0s
logvDNAin0v = [-4:0.5:4];
vDNAin0 = 10.^logvDNAin0v;
vDNAin0(9) = 3;
vDNAin0(10) = 8.5;
vDNAin0(11) = 14;
vDNAin0(13) = 131;
%vDNAin0 = [0.3,3,14,131,300,1000];

% Initial conditions
Protein10 = 0; Protein20 = 0; Capsid0 = 0;
Particle0 = 0; Virus0 = 1E-15;
y0 = [Protein10,Protein20,Capsid0,Particle0,Virus0];

tspan = [0:1:96];
options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'InitialStep',1e-2,...
    'NonNegative',(1:5), 'MaxOrder',5, 'BDF','on', 'Stats','off');

% Simulate
for i = 1:length(vDNAin0)
    vDNAin0_curr = vDNAin0(i);
    
    % Calculate model values for each MOI
    ODE_FH = @(t,y) Model(t,y,mpars,vDNAin0_curr);
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
    
    % Storing Models
    Protein_1_Models{2,i} = Protein1;
    Protein_2_Models{2,i} = Protein2;
    Capsid_Models{2,i} = Capsid;
    Particle_Models{2,i} = Particle;
    Virus_Models{2,i} = Virus;
end

% Maxima defined previously (see lines 82-85)

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

% Plot
figure(3); set(figure(3),'Units','inches','Position',[0.5 0.5 5.5 3.75]);
figure(4); set(figure(4),'Units','inches','Position',[0.5 0.5 5.5 3.75]);
figure(5); set(figure(5),'Units','inches','Position',[0.5 0.5 5.5 3.75]);
for i = 1:length(vDNAin0)
    vDNAin0_curr = vDNAin0(i); 
    
    % Tp5_1_tot Figure
    figure(3);
    % Plot Tp5_1
    plot(SumP1s{1,i},SumP1s{2,i}/GMax_SP1,'-','LineWidth',2);
    hold on;

    % Tp5_2_tot Figure
    figure(4);
    % Plot Tp5_2
    plot(SumP2s{1,i},SumP2s{2,i}/GMax_SP2,'-','LineWidth',2); 
    hold on; 

    % Virus Figures
    figure(5);
    % Plot Virus
    
    if i == 9
        lc = line_colors(1);
        semilogy(Virus_Models{1,i},Virus_Models{2,i}/GMax_Virus,...
            lc,'LineWidth',2);
    elseif i == 11
        lc = line_colors(2);
        semilogy(Virus_Models{1,i},Virus_Models{2,i}/GMax_Virus,...
            lc,'LineWidth',2);
    elseif i == 13
        lc = line_colors(3);
        semilogy(Virus_Models{1,i},Virus_Models{2,i}/GMax_Virus,...
            lc,'LineWidth',2);
    else
        semilogy(Virus_Models{1,i},Virus_Models{2,i}/GMax_Virus,'-k',...
            'LineWidth',2); 
    end
    hold on; 
end

figure(3);
%xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
%ylabel('Tp5_{1,tot} (A.U.)','fontsize',18,'FontName','Arial');
axis([0 96 0 1.2]);
xticks([0:24:96]);
set(gca,'fontsize',16,'FontName','Arial');
box off;

figure(4);
%xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
%ylabel('Tp5_{2,tot} (A.U.)','fontsize',18,'FontName','Arial');
axis([0 96 0 1.2]);
xticks([0:24:96]);
set(gca,'fontsize',16,'FontName','Arial');
box off;

figure(5);
%xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
%ylabel('Virus (A.U.)','fontsize',18,'FontName','Arial');
errorbar(Virustdata,VirusData(:,1),VirusSDs(:,1)/7600000,shapes(1),...
    'MarkerSize',8,'MarkerFaceColor',colors(1),'MarkerEdgeColor',colors(1),...
    'LineWidth',1.5,'Color',colors(1),'CapSize',6); hold on;

errorbar(Virustdata,VirusData(:,2),VirusSDs(:,2)/7600000,shapes(2),...
    'MarkerSize',8,'MarkerFaceColor',colors(2),'MarkerEdgeColor',colors(2),...
    'LineWidth',1.5,'Color',colors(2),'CapSize',6); hold on;

errorbar(Virustdata,VirusData(:,3),VirusSDs(:,3)/7600000,shapes(3),...
    'MarkerSize',8,'MarkerFaceColor',colors(3),'MarkerEdgeColor',colors(3),...
    'LineWidth',1.5,'Color',colors(3),'CapSize',6); hold on;
axis([0 96 0 1]);
xticks([0:24:96]);
set(gca,'fontsize',16,'FontName','Arial');
box off;

figure(3);
set(gcf,'color','white');
figure(4);
set(gcf,'color','white');
figure(5);
set(gcf,'color','white');