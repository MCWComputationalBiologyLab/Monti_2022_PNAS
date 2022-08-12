function [JacobTp5_1,JacobTp5_2,JacobVirus] = Jacobian(mpar,Tdata1,Tdata2,vDNAin0)

% Initial conditions
Tp5_10 = 0; Tp5_20 = 0; Capsid0 = 0; Particle0 = 0; Virus0 = 1E-15;
y0 = [Tp5_10,Tp5_20,Capsid0,Particle0,Virus0]; %IC

options = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'InitialStep',1e-2,...
    'NonNegative',(1:5), 'MaxOrder',5, 'BDF','on', 'Stats','off');

% Calculate Jacobian matrix numerically by central difference method
JacobTp5_1 = zeros(length(Tdata1),length(mpar));
JacobTp5_2 = zeros(length(Tdata1),length(mpar));
JacobVirus = zeros(length(Tdata2),length(mpar));

dmpar = 1e-2; % Fractional change in each parameter values
for i = 1:length(mpar)
     % Plus end
     mpar1 = mpar;
     mpar1(i) = mpar(i) + dmpar*mpar(i);

    % Calculate model values for plus end mpars
    ODE_FH = @(t,y) Model(t,y,mpar1,vDNAin0);
    sol1 = ode15s(ODE_FH,Tdata1,y0,options); 
    %Sufficient for Tp5_1 and Tp5_2 because they have the same tspan
    sol2 = ode15s(ODE_FH,Tdata2,y0,options); %Has different tspan
    %Calculate model for time points a/w Tp5
    y1 = deval(Tdata1,sol1);
    %Calculate model for time points a/w Virus
    y2 = deval(Tdata2,sol2);

    % Assigning Models
    Tp5_1_ModelUp = y1(1,:)';
    Tp5_2_ModelUp = y1(2,:)';
    Virus_ModelUp = y2(5,:)';
    
    % Minus end
    mpar2 = mpar;
    mpar2(i) = mpar(i) - dmpar*mpar(i);

    % Calculate model values for plus end mpars
    ODE_FH = @(t,y) Model(t,y,mpar2,vDNAin0);
    sol1 = ode15s(ODE_FH,Tdata1,y0,options); 
    %Sufficient for Tp5_1 and Tp5_2 because they have the same tspan
    sol2 = ode15s(ODE_FH,Tdata2,y0,options); %Has different tspan
    %Calculate model for time points a/w Tp5
    y1 = deval(Tdata1,sol1);
    %Calculate model for time points a/w Virus
    y2 = deval(Tdata2,sol2);

    % Assigning Models
    Tp5_1_ModelDown = y1(1,:)';
    Tp5_2_ModelDown = y1(2,:)';
    Virus_ModelDown = y2(5,:)';

    % Calculate Jacobian
    JacobTp5_1(:,i) = (Tp5_1_ModelUp - Tp5_1_ModelDown)/(2*dmpar*mpar(i));
    JacobTp5_2(:,i) = (Tp5_2_ModelUp - Tp5_2_ModelDown)/(2*dmpar*mpar(i));
    JacobVirus(:,i) = (Virus_ModelUp - Virus_ModelDown)/(2*dmpar*mpar(i));
end