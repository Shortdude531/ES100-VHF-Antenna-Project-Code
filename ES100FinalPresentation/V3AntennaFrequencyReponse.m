close all; clear all; clc; clf;
%% % read S11 Experimental Data from CSv
S11 = readtable('MyFilePrefix_S11_S12','NumHeaderLines',0);
datasize = size(S11);

%Delete Keysight Fieldfox Header
S11(1:17,:) = [];
S11(size(S11),:) = [];

%1008 to 2009


S11 = S11{:,:};

fdata = S11(1:1019,1);
fs12data = S11(1008:2009,1);
explen = length(fdata);
S11data = S11(1:1019,2);
S11EXP = (1-10.^(S11data/10));
S12data = S11(1008:2009,2);
S12EXP = (1-10.^(S12data/10));

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
% plot(frequencies,S11data,fs12data,S12data,'LineWidth', 2.5)
plot(fs12data,S12data,'LineWidth', 2.5)
colororder([0 0.4470 0.7410;0.9290 0.6940 0.1250])


ylim([-40 0])
xlim([f1 f2])
xlabel("Frequency Hz")
ylabel("S11 Reflected Power dB")
title("S12 Version 3 Board Frequency Response")
set(findall(gcf,'-property','FontSize'),'FontSize',15)


xline(54e6,'--','LineWidth',1.75)
xline(66e6,'--','LineWidth',1.75)

legend('S12','54 MHz','66 MHz')
hold off

set(findall(gcf,'-property','FontSize'),'FontSize',30)