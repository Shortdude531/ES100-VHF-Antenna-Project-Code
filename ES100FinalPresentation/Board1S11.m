close all; clear all; clc; clf;
%% define dipole (m)
d1 = dipoleCylindrical('Radius', 0.0006455,'Conductor',metal('Copper'),'Load',lumpedElement('Impedance',50),'ClosedEnd',1,"Length",2.39);
font=15;
figure()
show(d1)
set(findall(gcf,'-property','FontSize'),'FontSize',30)




%% % read S11 Experimental Data from CSv
S11 = readtable('S11BoardV1','NumHeaderLines',0);
datasize = size(S11);

%Delete Keysight Fieldfox Header
S11(1:17,:) = [];
S11(size(S11),:) = [];

S11 = S11{:,:};

fdata = S11(:,1);
explen = length(fdata);
S11data = S11(:,2);
S11EXP = (1-10.^(S11data/10));

f1 = 54e6;
f2 = 66e6;
fLower = 30e6;  %bandwidth lower frequency
fUpper = 90e6;  %bandwidth higher frequency
FCenter = 60e6; %Center frequency
f_interval = .1e6; %frequency interval
frequencies = fdata;

%% Calculate Simalted SParameters 
s1= sparameters(d1,frequencies,50);
%% plot S11
figure()
rfplot(s1);


D=get(gca,'Children'); %get the handle of the line object
XData=get(D,'XData'); %get the x data
YData=get(D,'YData'); %get the y data
Data=[XData' YData']; %join the x and y data on one array nx2


figure()
hold on 
%S11EXP = (1-10.^(S11data/10));
%S11THE = (1-10.^(YData/10));
plot(frequencies,S11data,frequencies,YData,'LineWidth', 2.5)
ylim([-50 0])
xlim([f1 f2])
xlabel("Frequency Hz")
ylabel("S11 Reflected Power DB")
title("Version 1 Board S11 Frequency Reponse")
set(findall(gcf,'-property','FontSize'),'FontSize',15)


xline(54e6,'--','LineWidth',1.75)
xline(66e6,'--','LineWidth',1.75)
yline(-3,'--','LineWidth',1.75 ,'Color',[0.6350 0.0780 0.1840]);
legend('Experimental S11','Simulated S11','54 MHz','66 MHz','-3 dB 50% Reflected Power')
hold off

set(findall(gcf,'-property','FontSize'),'FontSize',30)