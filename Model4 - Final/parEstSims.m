function varargout = parEstSims(varargin)

mpars = varargin{1};
vDNAin0 = varargin{2};

% Place to store models after they are calculated
Protein_1_Models = cell(2,length(vDNAin0));
Protein_2_Models = cell(2,length(vDNAin0));
Capsid_Models = cell(2,length(vDNAin0));
Particle_Models = cell(2,length(vDNAin0));
Virus_Models = cell(2,length(vDNAin0));

% Initial conditions
Protein10 = 0; Protein20 = 0; Capsid0 = 0; Particle0 = 0; Virus0 = 1E-15;
y0 = [Protein10,Protein20,Capsid0,Particle0,Virus0];

% Time span for simulations
tspan = [0:1:96];

tspan = [0:1:96];
options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'InitialStep',1e-2,...
    'NonNegative',(1:5), 'MaxOrder',5, 'BDF','on', 'Stats','off');

% Simulate results
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
    
    % Storing Models
    Protein_1_Models{2,i} = Protein1;
    Protein_2_Models{2,i} = Protein2;
    Capsid_Models{2,i} = Capsid;
    Particle_Models{2,i} = Particle;
    Virus_Models{2,i} = Virus;
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

if (nargout == 3)
    varargout = {SumP1s,SumP2s,Virus_Models}; % Useful for legacy drivers
elseif(nargout == 4)
    varargout = {GMax_SP1,GMax_SP2,GMax_Virus,1};
else
% Assign Data
    vDNAData = varargin{3};
    Data1 = varargin{4};
    Data2 = varargin{5};
    Data3 = varargin{6};
    figNum = varargin{7};
    opts = varargin{8};
    Tp5_1tdata = Data1(:,1); %xdata
    Tp5_1Data = Data1(:,2:length(Data1(1,:))); %ydata
    Tp5_2tdata = Data2(:,1); %xdata
    Tp5_2Data = Data2(:,2:length(Data2(1,:))); %ydata
    Virustdata = Data3(:,1); %xdata
    VirusData = Data3(:,2:length(Data3(1,:))); %ydata

    % Load in SDs
    Tp5_1SDs = load('Tp5_1SDs.txt');
    Tp5_2SDs = load('Tp5_2SDs.txt');
    VirusSDs = load('VirusSDs.txt');
    fittedSDs = load('RM_FittedData_01 05 and 50_SD.txt');


    % Plot model simulations compared to experimental data
    figure(figNum); set(figure(figNum),'Units','inches','Position',...
        [0.5 0.5 16 4])
    
    if strcmp(opts,'FirstIter') == 1
        ls{1} = '--b'; ls{2} = '--g'; ls{3} = '--r';
    else
        ls{1} = '-b'; ls{2} = '-g'; ls{3} = '-r';
    end
    
    c{1} = 'b'; c{2} = 'g'; c{3} = 'r';
    shapes{1} = 'o'; shapes{2} = '^'; shapes{3} = 's';
    %Plot data and model simulations
    for k = 1:length(vDNAin0)
        vDNAin0_curr = vDNAin0(k);
        
        % Isolating normalized data
        Tp5_1Data_curr = Tp5_1Data(:,k);
        Tp5_2Data_curr = Tp5_2Data(:,k);
        VirusData_curr = VirusData(:,k);
        
        % vDNA Figure
%         subplot(2,2,1);
%         [vDNAtot] = Calc_vDNA(tspan,vDNAin0_curr);
%         % Plot vDNA model
%         semilogy(tspan,vDNAtot,ls{k},'LineWidth',2); hold on;
%         %Plot vDNA data
%         errorbar(vDNAData(:,1),vDNAData(:,k+1),fittedSDs(:,k),shapes{k},...
%         'MarkerSize',7,'MarkerFaceColor',c{k},'MarkerEdgeColor',c{k},...
%         'LineWidth',2,'Color',c{k},'CapSize',9); hold on;
        
        % Tp5_1 Figure
        subplot(1,3,1);
        % Plot Tp5_1 model
        plot(SumP1s{1,k},SumP1s{2,k}/GMax_SP1,ls{k},'LineWidth',2); hold on
        %Plot Tp51 data
        errorbar(Tp5_1tdata,Tp5_1Data_curr,Tp5_1SDs(:,k),shapes{k},...
        'MarkerSize',7,'MarkerFaceColor',c{k},'MarkerEdgeColor',c{k},...
        'LineWidth',2,'Color',c{k},'CapSize',9); hold on;
        
        % Tp5_2 Figure
         subplot(1,3,2);
        % Plot Tp5_2 model
        plot(SumP2s{1,k},SumP2s{2,k}/GMax_SP2,ls{k},'LineWidth',2); hold on
        %Plot Tp52 data
        errorbar(Tp5_2tdata,Tp5_2Data_curr,Tp5_2SDs(:,k),shapes{k},...
        'MarkerSize',7,'MarkerFaceColor',c{k},'MarkerEdgeColor',c{k},...
        'LineWidth',2,'Color',c{k},'CapSize',9); hold on;
        
        % Virus Figure
         subplot(1,3,3);
        % Plot Virus model
        semilogy(Virus_Models{1,k},Virus_Models{2,k}/GMax_Virus,ls{k},...
            'LineWidth',2); hold on
        %Plot Virus data
        errorbar(Virustdata,VirusData_curr,VirusSDs(:,k),shapes{k},...
        'MarkerSize',7,'MarkerFaceColor',c{k},'MarkerEdgeColor',...
        c{k},'LineWidth',2,'Color',c{k},'CapSize',9); hold on;
    
    end
    
    %Figure formatting
%      subplot(1,3,1);
%     set(gca,'FontSize',14)
%     xlabel('Time (hpi)','fontsize',18);
%     ylabel('vDNA_{tot}','fontsize',18); box off;
%     axis([0 96 10^-1 10^4]);
%     xticks(0:24:96);
    
    subplot(1,3,1);
    set(gca,'FontSize',14)
    xlabel('Time (hpi)','fontsize',18);
    ylabel('Tp_{5.1} (A.U.)','fontsize',18); box off;
    axis([0 96 0 1.2]);
    xticks(0:24:96);
    
    subplot(1,3,2);
    set(gca,'FontSize',14)
    xlabel('Time (hpi)','fontsize',18);
    ylabel('Tp_{5.2} (A.U.)','fontsize',18); box off;
    axis([0 96 0 1.2]);
    xticks(0:24:96);
    
    subplot(1,3,3);
    set(gca,'FontSize',14)
    xlabel('Time (hpi)','fontsize',18);
    ylabel('Virus (A.U.)','fontsize',18); box off;
    axis([0 96 0 1.2]);
    xticks(0:24:96);
    
    set(gcf,'color','white');
end

end % function end