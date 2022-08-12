%% Driver
clear; clc; clf; close all;
format short e;
warning('off','all')
rng(1);

%% Load in data
vDNAData = load('vDNAData.txt');
vDNAin0 = vDNAData(1,2:length(vDNAData(1,:)));
Data1 = load('Tp5_1Data.txt');
Data2 = load('Tp5_2Data.txt');
Data3 = load('VirusData.txt');

disp("==================================================================");

%% Block of code for testing only
load('Pars2Set.mat');
load('Resids2Set.mat');
M = 200;
numPars = 13;
for i = 1:length(pars2)
    curr_pars = pars2{i};
    parEstSims(curr_pars,vDNAin0,vDNAData,Data1,Data2,Data3,1,'');
    % Will overwrite whatever is in figure 1!!
end


%% Minimum of Minimum SSE Strategy
% Find the index and value of SSE with lowest value
[mresids,idx] = min(resids2);

% Assign and save the parameters with the minimal SSE
mpars = pars2{idx};

% Plot final fittings based on minimum parameters in a separate figure
parEstSims(mpars,vDNAin0,vDNAData,Data1,Data2,Data3,2,''); hold on;

%% Average of Minimum SSE Strategy
% Extract each parameter set from pars2 and put it into a matrix
parsMat = [];
for i = 1:M
    currPars = pars2{i};
    parsMat = [parsMat; currPars];
end

% Find average of each parameter
mpars = [];
parsSDs = [];
parsSEs = [];
for i = 1:numPars
    currParsAvg = mean(parsMat(:,i));
    currParsSDs = std(parsMat(:,i));
    currParsSEs = currParsSDs/sqrt(M);
    mpars = [mpars currParsAvg];
    parsSDs = [parsSDs currParsSDs];
    parsSEs = [parsSEs currParsSEs];
end
residsAvg = mean(resids2);
residsSDs = std(resids2);
residsSEs = residsSDs/sqrt(M);

% Plot final fittings based on average parameters in a separate figure
parEstSims(mpars,vDNAin0,vDNAData,Data1,Data2,Data3,2,'FirstIter'); 

%% Save figures

% saveas(figure(1),'ModFits_Fig1.fig');
% saveas(figure(1),'ModFits_Fig1.tif');
% saveas(figure(2),'ModFits_Fig2.fig');
% saveas(figure(2),'ModFits_Fig2.tif');