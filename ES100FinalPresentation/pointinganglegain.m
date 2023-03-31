close all; clear all; clc; clf;
theta = linspace(0,pi,1000);
LTheta = length(theta);
L_c=5;
k_a= (2*pi)/5;
h=3/2;

for i = 1:LTheta
    g(i) = abs(((cos(k_a*h*cos(theta(i)))-cos(k_a*h))/sin(theta(i))));
end

g = 20.*log(g);


figure()
hold on 
%S11EXP = (1-10.^(S11data/10));
%S11THE = (1-10.^(YData/10));


xlabel("Angle (Rad)")
ylabel("Gain dBi")
title("Cylindrical Dipole Gain Vs. Pointing Angle")
set(findall(gcf,'-property','FontSize'),'FontSize',15)
plot(theta,g,'LineWidth', 2.5)

xline(1.309,'--','LineWidth',1.75)
xline(1.8326,'--','LineWidth',1.75)
yline(2.12,'--','LineWidth',1.75 ,'Color',[0.6350 0.0780 0.1840]);

legend('Antenna Gain dBi','-15°','+15°','2.12dBi')
hold off

 set(findall(gcf,'-property','FontSize'),'FontSize',30)