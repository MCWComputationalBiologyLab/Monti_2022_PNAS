t = [0:1:10];

a = t + 3;
b = 3*t - 9;
c = -2*t + 6;

plot(t,a,'-b','LineWidth',2); hold on;
plot(t,b,'-g','LineWidth',2); hold on;
plot(t,c,'-r','LineWidth',2); hold on;

legend("vDNA_{in,0} = 3","vDNA_{in,0} = 14","vDNA_{in,0} = 131",'FontName',...
    'Arial','FontSize',18,'Box','off'); 

