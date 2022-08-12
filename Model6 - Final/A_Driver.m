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

%% fmincon optimization
options = optimset('fmincon');
%options = optimset(options,'algorithm','sqp');
options = optimset(options,'TolFun',1e-6,'TolX',1e-6);
options = optimset(options,'Display','off');
%options = optimset(options,'Maxiter',100,'MaxFunEvals',1000);
%options = optimset(options,'UseParallel','always');
%options = []; % We want to use only default optimizer options

numPars = 13;  % Number of parameters

% Initial guesses for parameters (order of magnitudes to be appropriate)
p0 = [6E2,9E2,1E-1,5E4,1E1,6E2,9E3,1E-1,2.5E4,1E1,1E-6,3E-5,1.5E0];

% Define parameters for initial guess randomization
delta = 0.2; %amount to change by
rNumLb = 1-delta; %lower bound
rNumUb = 1+delta; %upper boud

% Define parameters for loops
M = 200;
pars2 = cell(1,M);
resids2 = zeros(1,M);

tic 
for j = 1:M
    
    N = 10;
    pars = cell(1,N);
    resids = zeros(1,N);
    
    parfor i = 1:N
        rNum = rNumLb + (rNumUb-rNumLb)*rand(1,length(p0)); 
        %generates random number between 1-delta and 1+delta
        p = p0.*rNum; %p0 randomized between (1-delta)*p0 and (1+delta)*p0
        p(3) = p0(3); p(8) = p0(8); p(13) = p0(13);
    
        % Bounds. Parameters values can range from 0.01*p0 to 100*p0.
        lbp = p/10; ubp = p*10;
        lbp(3) = p(3)*1; ubp(3) = p(3)*1;
        lbp(8) = p(8)*1; ubp(8) = p(8)*1;
        lbp(13) = p(13)*1; ubp(13) = p(13)*1;
       
        % Parameter estimation
        [curr_pars,fun_val] = fmincon(@Error,p,[],[],[],[],...
            lbp,ubp,@mycon,options,Data1,Data2,Data3,vDNAin0);
    
        % Store parameters from this iteration
        pars{i} = curr_pars;
        resids(i) = fun_val;
    
%         % Display the current iteration number
%         fprintf("Inner Iteration %d \n",i);
%         fprintf('The current parameter set is: \n\n'); 
%         disp(curr_pars);
%         fprintf('With SSE: %.3f \n\n',fun_val);
    end

    %% Min SSE Strategy
    % Find the index and value of SSE with lowest value
    [mresids,idx] = min(resids);
    
    % Assign and save the parameters with the minimal SSE
    mpars = pars{idx};    
    
    % Display
    disp("..............................................................");
    fprintf('For outer iteration %d, the minimum parameter set is: \n\n', j); 
    disp(mpars);
    fprintf('With SSE: %.3f \n\n',mresids);

    pars2{j} = mpars;  % Store mpars for outer (j) loop
    resids2(j) = mresids;  % Store mresids for outer (j) loop
    p0 = mpars;  % initialize parameters for outer (j) loop
       
    %% Plotting model simulations and data after outer loop
    parEstSims(mpars,vDNAin0,vDNAData,Data1,Data2,Data3,1,'');
end

save('Pars2Set.mat','pars2');
save('Resids2Set.mat','resids2');
disp("==================================================================");
toc

%% Block of code for testing only
% load('Pars2Set.mat');
% load('Resids2Set.mat');
M = 200;
numPars = 13;
% for i = 1:length(pars2)
%     curr_pars = pars2{i};
%     parEstSims(curr_pars,vDNAin0,vDNAData,Data1,Data2,Data3,1,'');
%     % Will overwrite whatever is in figure 1!!
% end

if exist('Driver_Output.txt','file') ~= 0
    delete('Driver_Output.txt'); % get rid of old file if there is one
end
diary Driver_Output.txt;

% Initial guesses for parameters and SSE calculation
p0 = [6E2,9E2,1E-1,5E4,1E1,6E2,9E3,1E-1,2.5E4,1E1,1E-6,3E-5,1.5E0]; % Initial guess
Err_IG = Error(p0,Data1,Data2,Data3,vDNAin0); % SSE from initial guess
fprintf('The SSE when using the initial guess is: %.3f \n\n', Err_IG);

% Calculate corrected AIC for the model based on initial SSE
SSE = Err_IG;
Nd = 45; Np = numPars-3+1;
AIC = Nd*log(SSE/Nd) + 2*Np + (2*Np*(Np+1))/(Nd-Np-1);
fprintf('The AIC for the initial parameter set is: %.3f \n\n', AIC);

%% Minimum of Minimum SSE Strategy
% Find the index and value of SSE with lowest value
[mresids,idx] = min(resids2);

% Assign and save the parameters with the minimal SSE
mpars = pars2{idx};
save('Model_mpars_min.mat','mpars');

disp('The minimum parameter set is: '); 
disp(mpars);
fprintf('With SSE: %.3f \n\n', mresids);

% Calculate corrected AIC for the model based on minimum
SSE = mresids;
Nd = 45; Np = numPars-3+1;
AIC = Nd*log(SSE/Nd) + 2*Np + (2*Np*(Np+1))/(Nd-Np-1);
fprintf('The AIC for the minimum parameter set is: %.3f \n\n', AIC);

% Plot final fittings based on minimum parameters in a separate figure
parEstSims(mpars,vDNAin0,vDNAData,Data1,Data2,Data3,2,'');

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
save('Model_mpars_avg.mat','mpars');
save('Model_mpars_std.mat','parsSDs');
save('Model_mpars_sem.mat','parsSEs');

disp('The average parameter set is: ');
disp(mpars);
disp('The parameter set SD is: ');
disp(parsSDs);
disp('The parameter set SE is: ');
disp(parsSEs);

fprintf('The average SSE is: %.3f \n\n', residsAvg);
fprintf('The SSE SD is: %d \n\n', residsSDs);
fprintf('The SSE SE is: %d \n\n', residsSEs);

% Calculate corrected AIC for the model based on average SSE
SSE = residsAvg;
Nd = 45; Np = numPars-3+1;
AIC = Nd*log(SSE/Nd) + 2*Np + (2*Np*(Np+1))/(Nd-Np-1);
fprintf('The AIC for the average parameter set is: %.3f \n\n', AIC);

% Calculate actual SSE based on average parameters
Err_avgPars = Error(mpars,Data1,Data2,Data3,vDNAin0);
fprintf('The SSE using the average parameter set is: %.3f \n\n', Err_avgPars);

% Calculate corrected AIC for the model based on average SSE
SSE = Err_avgPars;
Nd = 45; Np = numPars-3+1;
AIC = Nd*log(SSE/Nd) + 2*Np + (2*Np*(Np+1))/(Nd-Np-1);
fprintf('The AIC for the average parameter set is: %.3f \n\n', AIC);

% Plot final fittings based on average parameters in a separate figure
parEstSims(mpars,vDNAin0,vDNAData,Data1,Data2,Data3,2,'FirstIter');

diary off
save

%% Save figures
figure(1)
sgtitle('Outer Loop Model Fits','FontName','Arial','FontWeight','normal');
figure(2);
sgtitle('Final Model Fits','FontName','Arial','FontWeight','normal');

saveas(figure(1),'ModFits_Fig1.fig');
saveas(figure(1),'ModFits_Fig1.tif');
saveas(figure(2),'ModFits_Fig2.fig');
saveas(figure(2),'ModFits_Fig2.tif');