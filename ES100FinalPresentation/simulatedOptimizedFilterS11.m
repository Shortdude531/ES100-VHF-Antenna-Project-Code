close all; clear all; clc; clf;
%% % read S11 Experimental Data from CSv
S11 = readtable('SimulatedOptimizedFilter','NumHeaderLines',0);
datasize = size(S11);

%Delete Keysight Fieldfox Header
S11(1,:) = [];
S11(size(S11),:) = [];

%1008 to 2009


S11 = S11{:,:};

len = length(S11);

fdata = S11(1:len,1);
explen = len;


S11data = S11(1:len,3);


S11EXP = (1-10.^(S11data/10));
S11EXP = smoothdata(S11EXP,"gaussian",20);

S12data = S11(1:len,2);

S12EXP = (1-10.^(S12data/10));
S12EXP = smoothdata(S12EXP,"gaussian",20);
f1 = 30e6;
f2 = 90e6;
fLower = 30e6;  %bandwidth lower frequency
fUpper = 90e6;  %bandwidth higher frequency
FCenter = 60e6; %Center frequency
f_interval = .1e6; %frequency interval
frequencies = fdata;


figure()
hold on 
%S11EXP = (1-10.^(S11data/10));
%S11THE = (1-10.^(YData/10));
% plot(frequencies,S11data,frequencies,S12data,'LineWidth', 2.5)
plot(frequencies,S11data,'LineWidth', 2.5)
colororder([0 0.4470 0.7410;0.9290 0.6940 0.1250])
%curveFitter(frequencies,S11data)

ylim([-40 0])
xlim([f1 f2])
xlabel("Frequency Hz")
ylabel("S11 Reflected Power DB")
title("Simulated Optimized Chebyshev Transformer S11 Frequency Response")
set(findall(gcf,'-property','FontSize'),'FontSize',15)


xline(54e6,'--','LineWidth',1.75)
xline(66e6,'--','LineWidth',1.75)
yline(-10,'--','LineWidth',1.75 ,'Color',[0.6350 0.0780 0.1840]);
legend('S11','54 MHz','66 MHz','-10 dB 90% Impedance Mismatch Efficiency')
hold off

set(findall(gcf,'-property','FontSize'),'FontSize',30)

%curveFitter(frequencies,S11data)

