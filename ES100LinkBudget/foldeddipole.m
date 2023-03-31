clc; clear all; close all;
%%Folded Dipole Radiation pattern and S11

f = 60*10^6;
c= 3e8;
c1 = .5;
c_offsetLegnth = -.255;
lambda =c/f;

L1 =.5*lambda+c_offsetLegnth



df = dipoleFolded('Length',2.2,'Width',0.1,'Load',lumpedElement('Impedance',75),'Spacing',2.23/51)


figure(1)
pattern(df,f), title('radiation a=pattern of dipole')
% 
% 
% impedance(d1,50e6:1e6:70e6)
% figure(2)
%range Mhz
f1 = 54e6;
f2 = 66e6;
f_interval = .1e6;
S = sparameters(df,f1:f_interval:f2);
% S1 = sparameters(d2,f1:1e6:f2);
% S2 = sparameters(d3,f1:1e6:f2);
figure(2)
title(sprintf('One dipole antennas with a radius of 10cm'))
hold on
rfplot(S)
% rfplot(S1)
% rfplot(S2)
hold off