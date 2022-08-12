%% Stats Driver
clear; clc; close all;
%% Load data
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

Maxes = DetermineMaxes(vDNAin0,1);
SSE1 = Calc_SSE(@Model1,Data1,Data2,Data3,vDNAin0,Maxes);

Maxes = DetermineMaxes(vDNAin0,2);
SSE2 = Calc_SSE(@Model2,Data1,Data2,Data3,vDNAin0,Maxes);

Maxes = DetermineMaxes(vDNAin0,3);
SSE3 = Calc_SSE(@Model3,Data1,Data2,Data3,vDNAin0,Maxes);

Maxes = DetermineMaxes(vDNAin0,4);
SSE4 = Calc_SSE(@Model4,Data1,Data2,Data3,vDNAin0,Maxes);

Maxes = DetermineMaxes(vDNAin0,5);
SSE5 = Calc_SSE(@Model5,Data1,Data2,Data3,vDNAin0,Maxes);

Maxes = DetermineMaxes(vDNAin0,6);
SSE6 = Calc_SSE(@Model6,Data1,Data2,Data3,vDNAin0,Maxes);


%N = 48; % Number of data points.
N = 45;
% Protein data is 6 time points for 3 MOIs and 2 different proteins
% Virus data is 4 time points for 3 MOIs. Total = 48.

% Number of parameters
P1 = 9;
P2 = 11;
P3 = 11;
P4 = 11;
P5 = 13;
P6 = 13;

% Number of parameters + 1
% K1 = 10;
% K2 = 12;
% K3 = 12;
% K4 = 12;
% K5 = 14;
% K6 = 14;

% Two parameters kept constant 
K1 = 8;
K2 = 10;
K3 = 10;
K4 = 10;
K5 = 12;
K6 = 12;

% Degrees of freedom
DF1 = N-P1;
DF2 = N-P2;
DF3 = N-P3;
DF4 = N-P4;
DF5 = N-P5;
DF6 = N-P6;

%% AICs
AIC1 = AIC(SSE1,N,K1);
AIC2 = AIC(SSE2,N,K2);
AIC3 = AIC(SSE3,N,K3);
AIC4 = AIC(SSE4,N,K4);
AIC5 = AIC(SSE5,N,K5);
AIC6 = AIC(SSE6,N,K6);

AICs = [1,2,3,4,5,6;AIC1,AIC2,AIC3,AIC4,AIC5,AIC6]

AIC12 = AIC_Comp(AIC1,AIC2);
AIC13 = AIC_Comp(AIC1,AIC3);
AIC14 = AIC_Comp(AIC1,AIC4);
AIC15 = AIC_Comp(AIC1,AIC5);
AIC16 = AIC_Comp(AIC1,AIC6);

AIC21 = AIC_Comp(AIC2,AIC1);
AIC23 = AIC_Comp(AIC2,AIC3);
AIC24 = AIC_Comp(AIC2,AIC4);
AIC25 = AIC_Comp(AIC2,AIC5);
AIC26 = AIC_Comp(AIC2,AIC6);

AIC31 = AIC_Comp(AIC3,AIC1);
AIC32 = AIC_Comp(AIC3,AIC2);
AIC34 = AIC_Comp(AIC3,AIC4);
AIC35 = AIC_Comp(AIC3,AIC5);
AIC36 = AIC_Comp(AIC3,AIC6);

AIC41 = AIC_Comp(AIC4,AIC1);
AIC42 = AIC_Comp(AIC4,AIC2);
AIC43 = AIC_Comp(AIC4,AIC3);
AIC45 = AIC_Comp(AIC4,AIC5);
AIC46 = AIC_Comp(AIC4,AIC6);

AIC51 = AIC_Comp(AIC5,AIC1);
AIC52 = AIC_Comp(AIC5,AIC2);
AIC53 = AIC_Comp(AIC5,AIC3);
AIC54 = AIC_Comp(AIC5,AIC4);
AIC56 = AIC_Comp(AIC5,AIC6);

AIC61 = AIC_Comp(AIC6,AIC1);
AIC62 = AIC_Comp(AIC6,AIC2);
AIC63 = AIC_Comp(AIC6,AIC3);
AIC64 = AIC_Comp(AIC6,AIC4);
AIC65 = AIC_Comp(AIC6,AIC5);

AIC_Matrix = [0,AIC12,AIC13,AIC14,AIC15,AIC16;AIC21,0,AIC23,AIC24,AIC25,...
    AIC26;AIC31,AIC32,0,AIC34,AIC35,AIC36;AIC41,AIC42,AIC43,0,AIC45,...
    AIC46;AIC51,AIC52,AIC53,AIC54,0,AIC56;AIC61,AIC62,AIC63,AIC64,AIC65,0];

% Plot
figure(1); set(figure(1),'Units','inches','Position',[0.5 0.5 10 8])
imagesc(AIC_Matrix); hold on;
%colormap(flipud(gray)); 
colormap(flipud(pink));

% Add numbers to the boxes
textStrings = num2str(AIC_Matrix(:),'%0.2E');
textStrings = strtrim(cellstr(textStrings));
[x, y] = meshgrid(1:6);
hStrings = text(x(:), y(:), textStrings(:), ...
                'HorizontalAlignment', 'center','FontName','Arial',...
                'FontSize',14);
midValue = mean(get(gca, 'CLim'));
textColors = repmat(AIC_Matrix(:) > midValue, 1, 3);
set(hStrings, {'Color'}, num2cell(textColors, 2));

% Formatting
colorbar;
xtl = ["Model 1";"Model 2";"Model 3";"Model 4";"Model 5";"Model 6"];
xt = [1:1:length(xtl)];
ytl = ["Model 1";"Model 2";"Model 3";"Model 4";"Model 5";"Model 6"];
yt = [1:1:length(xtl)];
set(gca,'XTick',xt,'XTickLabel',xtl,'YTick',yt,'YTickLabel',ytl);
set(gcf,'color','white');
title('Probability of correctness of model n vs. model m (AIC)',...
    'FontName','Arial','FontSize',18,'FontWeight','Normal');
ylabel('Model n','FontName','Arial','FontSize',18);
xlabel('Model m','FontName','Arial','FontSize',18);
set(gca,'FontName','Arial','FontSize',18);

[M,I] = min(AICs(2,:));
disp('The minimum AIC is: ');
disp(num2str(M));
disp('The model with the minimum AIC is model number: ');
modNum = AICs(1,I);
disp(num2str(modNum));

%% F-tests
F12 = FTest(SSE1,SSE2,DF1,DF2);
F13 = FTest(SSE1,SSE3,DF1,DF3);
F14 = FTest(SSE1,SSE4,DF1,DF4);
F15 = FTest(SSE1,SSE5,DF1,DF5);
F16 = FTest(SSE1,SSE6,DF1,DF6);

F25 = FTest(SSE2,SSE5,DF2,DF5);
F26 = FTest(SSE2,SSE6,DF2,DF6);

F35 = FTest(SSE3,SSE5,DF3,DF5);
F36 = FTest(SSE3,SSE6,DF3,DF6);

F45 = FTest(SSE4,SSE5,DF4,DF5);
F46 = FTest(SSE5,SSE6,DF4,DF6);

p12 = fcdf(F12,DF1,DF2);
p13 = fcdf(F13,DF1,DF3);
p14 = fcdf(F14,DF1,DF4);
p15 = fcdf(F15,DF1,DF5);
p16 = fcdf(F16,DF1,DF6);

p25 = fcdf(F25,DF2,DF5);
p26 = fcdf(F26,DF2,DF6);

p35 = fcdf(F35,DF3,DF5);
p36 = fcdf(F36,DF3,DF6);

p45 = fcdf(F45,DF4,DF5);
p46 = fcdf(F46,DF4,DF6);

FMatrix = [0,1-p12,1-p13,1-p14,1-p15,1-p16;1-p12,0,0,0,1-p25,1-p26;...
    1-p13,0,0,0,1-p35,1-p36;1-p14,0,0,0,1-p45,1-p46;1-p15,1-p25,1-p35,...
    1-p45,0,0;1-p16,1-p26,1-p36,1-p46,0,0];

% Plot
figure(2); set(figure(2),'Units','inches','Position',[0.5 0.5 10 8])
imagesc(FMatrix); hold on;
%colormap(flipud(gray)); 
%colormap(flipud(pink));
colormap(pink);

% Add numbers to the boxes
textStrings = num2str(FMatrix(:),'%0.2E');
textStrings = strtrim(cellstr(textStrings));
textStrings = replace(textStrings,'0.00E+00','N/A');
[x, y] = meshgrid(1:6);
hStrings = text(x(:), y(:), textStrings(:), ...
                'HorizontalAlignment', 'center','FontName','Arial',...
                'FontSize',14);
midValue = mean(get(gca, 'CLim'));
% midValue = 0.05;
textColors = repmat(FMatrix(:) < midValue, 1, 3);
set(hStrings, {'Color'}, num2cell(textColors, 2));

% Formatting
colorbar;
xtl = ["Model 1";"Model 2";"Model 3";"Model 4";"Model 5";"Model 6"];
xt = [1:1:length(xtl)];
ytl = ["Model 1";"Model 2";"Model 3";"Model 4";"Model 5";"Model 6"];
yt = [1:1:length(xtl)];
set(gca,'XTick',xt,'XTickLabel',xtl,'YTick',yt,'YTickLabel',ytl);
set(gcf,'color','white');
title('P-value for model n vs. model m (F-Test)','FontName','Arial',...
    'FontSize',24,'FontWeight','Normal');
xlabel('Model m','FontName','Arial','FontSize',18);
ylabel('Model n','VerticalAlignment','bottom','FontName','Arial',...
    'FontSize',18);
set(gca,'FontName','Arial','FontSize',18);

%% Evidence Ratios
ER12 = EvRatio(AIC1,AIC2);
ER13 = EvRatio(AIC1,AIC3);
ER14 = EvRatio(AIC1,AIC4);
ER15 = EvRatio(AIC1,AIC5);
ER16 = EvRatio(AIC1,AIC6);

ER21 = EvRatio(AIC2,AIC1);
ER23 = EvRatio(AIC2,AIC3);
ER24 = EvRatio(AIC2,AIC4);
ER25 = EvRatio(AIC2,AIC5);
ER26 = EvRatio(AIC2,AIC6);

ER31 = EvRatio(AIC3,AIC1);
ER32 = EvRatio(AIC3,AIC2);
ER34 = EvRatio(AIC3,AIC4);
ER35 = EvRatio(AIC3,AIC5);
ER36 = EvRatio(AIC3,AIC6);

ER41 = EvRatio(AIC4,AIC1);
ER42 = EvRatio(AIC4,AIC2);
ER43 = EvRatio(AIC4,AIC3);
ER45 = EvRatio(AIC4,AIC5);
ER46 = EvRatio(AIC4,AIC6);

ER51 = EvRatio(AIC5,AIC1);
ER52 = EvRatio(AIC5,AIC2);
ER53 = EvRatio(AIC5,AIC3);
ER54 = EvRatio(AIC5,AIC4);
ER56 = EvRatio(AIC5,AIC6);

ER61 = EvRatio(AIC6,AIC1);
ER62 = EvRatio(AIC6,AIC2);
ER63 = EvRatio(AIC6,AIC3);
ER64 = EvRatio(AIC6,AIC4);
ER65 = EvRatio(AIC6,AIC5);

ER_Matrix = [0,ER12,ER13,ER14,ER15,ER16;ER21,0,ER23,ER24,ER25,...
    ER26;ER31,ER32,0,ER34,ER35,ER36;ER41,ER42,ER43,0,ER45,...
    ER46;ER51,ER52,ER53,ER54,0,ER56;ER61,ER62,ER63,ER64,ER65,0];

ER_Matrix_Log = [1,ER12,ER13,ER14,ER15,ER16;ER21,1,ER23,ER24,ER25,...
    ER26;ER31,ER32,1,ER34,ER35,ER36;ER41,ER42,ER43,1,ER45,...
    ER46;ER51,ER52,ER53,ER54,1,ER56;ER61,ER62,ER63,ER64,ER65,1];

% Plot
figure(3); set(figure(3),'Units','inches','Position',[0.5 0.5 10 8])
imagesc(log10(ER_Matrix_Log)); hold on;
%colormap(flipud(gray)); 
colormap(flipud(pink));

% Add numbers to the boxes
textStrings = num2str(ER_Matrix(:),'%0.2E');
textStrings = strtrim(cellstr(textStrings));
[x, y] = meshgrid(1:6);
hStrings = text(x(:), y(:), textStrings(:), ...
                'HorizontalAlignment', 'center','FontName','Arial',...
                'FontSize',14);
%midValue = mean(get(gca, 'CLim'));
midValue = log10(mean(mean(ER_Matrix)));
textColors = repmat(ER_Matrix(:) > midValue, 1, 3);
set(hStrings, {'Color'}, num2cell(textColors, 2));

% Formatting
c = colorbar;
lims = [-4:1:4];
c.Ticks = lims;
labels = cell(1,length(lims));
for i = 1:length(lims)
    %labels{1,i} = num2str(10^lims(i),'%0.0E');
    num = lims(i);
    string = sprintf('10^{%d}',num);
    labels{1,i} = string;
end
c.TickLabels = labels;
%c.TickLabels = 10.^lims;
xtl = ["Model 1";"Model 2";"Model 3";"Model 4";"Model 5";"Model 6"];
xt = [1:1:length(xtl)];
ytl = ["Model 1";"Model 2";"Model 3";"Model 4";"Model 5";"Model 6"];
yt = [1:1:length(xtl)];
set(gca,'XTick',xt,'XTickLabel',xtl,'YTick',yt,'YTickLabel',ytl);
set(gcf,'color','white');
title('Evidence ratio of model n vs. model m (AIC)',...
    'FontName','Arial','FontSize',18,'FontWeight','Normal');
ylabel('Model n','FontName','Arial','FontSize',18);
xlabel('Model m','FontName','Arial','FontSize',18);
set(gca,'FontName','Arial','FontSize',18);

rTest = runstest(SSE3);
disp(rTest);