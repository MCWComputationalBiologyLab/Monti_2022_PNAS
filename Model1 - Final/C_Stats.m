%% Stats Driver
warning('off','all');
%delete 'Stats.txt';

% Load in parameters 
load('Model_mpars_avg.mat'); %Change here
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

%% Find Jacobian matrix
Jacob = [];
for i = 1:length(vDNAin0)
   vDNAin0_curr = vDNAin0(i);
   [JTp5_1,JTp5_2,J_Virus] = Jacobian(mpars,Tp5_1tdata,Virustdata,...
       vDNAin0_curr);
   Jacob_curr = [JTp5_1;JTp5_2;J_Virus];
   Jacob = [Jacob;Jacob_curr];
end

%% Find Correlation Matrix
CM = Correlation_Matrix(Jacob,mpars);
CM = abs(CM);

% Plot
figure(3); set(figure(3),'Units','inches','Position',[0.5 0.5 10 8]);
imagesc(CM); hold on;
colormap(flipud(pink)); 

% Add numbers to the boxes
textStrings = num2str(CM(:),'%0.1E');
textStrings = strtrim(cellstr(textStrings));
[x, y] = meshgrid(1:length(mpars));
hStrings = text(x(:), y(:), textStrings(:), ...
                'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));
textColors = repmat(CM(:) > midValue, 1, 3);
set(hStrings, {'Color'}, num2cell(textColors, 2));

% Formatting
colorbar;
xtl = ["k_{s,1}";"K_{m,1}";"k_{d,1}";"k_{s,2}";"K_{m,2}";"k_{d,2}";...
    "k_{s,C}";"k_{s,P}";"k_{ex}"];
xt = [1:1:length(mpars)];
ytl = ["k_{s,1}";"K_{m,1}";"k_{d,1}";"k_{s,2}";"K_{m,2}";"k_{d,2}";...
    "k_{s,C}";"k_{s,P}";"k_{ex}"];
yt = [1:1:length(mpars)];
set(gca,'XTick',xt,'XTickLabel',xtl,'YTick',yt,'YTickLabel',ytl,...
    'FontName','Arial','FontSize',14,'TickDir','out');
box off
set(gcf,'color','white');
title('Correlation Matrix (Absolute Value)','FontName','Arial',...
    'FontWeight','normal','FontSize',18);

%% Sensitivity Analysis
Errs = SensitivityAnalysis(mpars,Data1,Data2,Data3,vDNAin0);

figure(7); set(figure(7),'Units','inches','Position',[0.5 0.5 12 8]);
names = ["k_{s,1}";"K_{m,1}";"k_{d,1}";"k_{s,2}";"K_{m,2}";"k_{d,2}";...
    "k_{s,C}";"k_{s,P}";"k_{ex}"];
cols = ["-k","-b","-k","-b","-k","-b","-k","-b","-r"];
for i = 1:length(Errs)
   xs_curr = Errs{1,i};
   idx = find(xs_curr == mpars(i));
   xs_curr = xs_curr/xs_curr(idx);
   ys_curr = Errs{2,i};
   ys_curr = ys_curr/ys_curr(idx);
   
   if (i == 1 || i == 4)
       subplot(2,2,1);
       plot(xs_curr,ys_curr,cols(i),'LineWidth',2); hold on;
       title('Synthesis Constants','FontName','Arial','FontWeight',...
           'Normal','FontSize',20);
       xlabel('p/p_{0}','FontName','Arial','FontSize',14);
       ylabel('SSE/SSE_{0}','FontName','Arial','FontSize',14);
   elseif(i == 2|| i == 5)
       subplot(2,2,2);
       plot(xs_curr,ys_curr,cols(i),'LineWidth',2); hold on;
       title('Michaelis-like Constants','FontName','Arial','FontWeight',...
           'Normal','FontSize',20);
       xlabel('p/p_{0}','FontName','Arial','FontSize',14);
       ylabel('SSE/SSE_{0}','FontName','Arial','FontSize',14);
   elseif(i == 3 || i == 6)
       subplot(2,2,3);
       plot(xs_curr,ys_curr,cols(i),'LineWidth',2); hold on;
       title('Degradation Constants','FontName','Arial','FontWeight',...
           'Normal','FontSize',20);
       xlabel('p/p_{0}','FontName','Arial','FontSize',14);
       ylabel('SSE/SSE_{0}','FontName','Arial','FontSize',14);
   else
       subplot(2,2,4);
       plot(xs_curr,ys_curr,cols(i),'LineWidth',2); hold on;
       title('Condensation Constants','FontName','Arial','FontWeight',...
           'Normal','FontSize',20);
       xlabel('p/p_{0}','FontName','Arial','FontSize',14);
       ylabel('SSE/SSE_{0}','FontName','Arial','FontSize',14);
   end
   
end

subplot(2,2,1);
legend(names(1),names(4),'FontName','Arial','FontSize',14,'Location',...
    'North','box','off');
axis([0.5 1.5 0.75 2.75]);
set(gca, 'FontName','Arial','FontSize',14,'box','off','XTick',...
    (0.5:0.25:1.5),'YTick',(0.75:0.5:2.75),'LineWidth',1.5);

subplot(2,2,2);
legend(names(2),names(5),...
    'FontName','Arial','FontSize',14,'Location','North','box','off');
axis([0.5 1.5 0.75 2.75]);
set(gca, 'FontName','Arial','FontSize',14,'box','off','XTick',...
    (0.5:0.25:1.5),'YTick',(0.75:0.5:2.75),'LineWidth',1.5);

subplot(2,2,3);
legend(names(3),names(6),'FontName','Arial','FontSize',14,'Location',...
    'North','box','off');
axis([0.5 1.5 0.75 2.75]);
set(gca, 'FontName','Arial','FontSize',14,'box','off','XTick',...
    (0.5:0.25:1.5),'YTick',(0.75:0.5:2.75),'LineWidth',1.5);

subplot(2,2,4);
legend(names(7),names(8),names(9),'FontName','Arial','FontSize',14,...
    'Location','North','box','off');
axis([0.5 1.5 0.75 2.75]);
set(gca, 'FontName','Arial','FontSize',14,'box','off','XTick',...
    (0.5:0.25:1.5),'YTick',(0.75:0.5:2.75),'LineWidth',1.5);

set(gcf,'color','white');