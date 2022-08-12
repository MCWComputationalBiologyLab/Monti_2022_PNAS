clear; clc; close all; 
format long g;

load('Model_mpars_avg.mat');
vDNAData = load('vDNAData.txt');
vDNAin0 = vDNAData(1,2:length(vDNAData(1,:)));
%% Obtain values for vDNAtot
% Time span for simulations
tspan = [0:1:96];

vDNAtots = cell(1,length(vDNAin0));
for i = 1:length(vDNAin0)
    vDNAin0_curr = vDNAin0(i);
    [vDNAtot] = Calc_vDNA(tspan,vDNAin0_curr);
    vDNAtots{1,i} = vDNAtot;
end
%% Obtain values for Tp51 and Tp52
% Place to store models after they are calculated
Protein_1_Models = cell(2,length(vDNAin0));
Protein_2_Models = cell(2,length(vDNAin0));
Capsid_Models = cell(2,length(vDNAin0));
Particle_Models = cell(2,length(vDNAin0));
Virus_Models = cell(2,length(vDNAin0));

% Initial conditions
Protein10 = 0; Protein20 = 0; Capsid0 = 0; Particle0 = 0; Virus0 = 1E-15;
y0 = [Protein10,Protein20,Capsid0,Particle0,Virus0];

tspan = [0:1:96];
options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'InitialStep',1e-2,...
    'NonNegative',(1:5), 'MaxOrder',5, 'BDF','on', 'Stats','off');

% Simulate results
for i = 1:length(vDNAin0)
    % Calculate model solution for each vDNAin0
    ODE_FH = @(t,y) Model(t,y,mpars,vDNAin0(i));
    sols1 = ode15s(ODE_FH,tspan,y0,options);
    y1 = deval(tspan,sols1);
    
   % Assigning models
    tout = tspan;
    Protein1 = y1(1,:);
    Protein2 = y1(2,:);
    Capsid = y1(3,:);
    Particle = y1(4,:);
    Virus = y1(5,:);
    
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

%% Plots
figure(1); set(figure(1),'Units','inches','Position',[0.5 0.5 12 4]);
set(gcf,'color','white');
cols = ["-b","-g","-r"];

% Km3 negative regulation
subplot(1,2,1);
Km3 = mpars(4);
for i = 1:length(vDNAin0)
    eq = Km3./(Km3 + Protein_1_Models{2,i});
    disp(num2str(eq(end)));
    plot(tspan,eq,cols(i),'LineWidth',2); hold on;
end
xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
ylabel({'Tp5_{1} Feedback'; 'Inhibition Strength (A.U.)'},'fontsize',18,...
    'FontName','Arial');
title('{\it R_{1,4}(t)} ({\itK_{m,3}})','FontName','Arial','FontSize',18,'FontWeight','normal');
axis([0 96 0 1]);
xticks([0:24:96]);
set(gca,'fontsize',16,'FontName','Arial','LineWidth',1);
box off;


% Km4 negative regulation
subplot(1,2,2);
Km4 = mpars(8);
for i = 1:length(vDNAin0)
    eq = Km4./(Km4 + Protein_2_Models{2,i});
    plot(tspan,eq,cols(i),'LineWidth',2); hold on;
end
xlabel('Time (hpi)','fontsize',18,'FontName','Arial');
ylabel({'Tp5_{2} Feedback'; 'Inhibition Strength (A.U.)'},'fontsize',18,...
    'FontName','Arial');
title('{\it R_{1,4}(t)} ({\itK_{m,4}})','FontName','Arial','FontSize',18,'FontWeight','normal');
axis([0 96 0 1]);
xticks([0:24:96]);
set(gca,'fontsize',16,'FontName','Arial','LineWidth',1);
box off;