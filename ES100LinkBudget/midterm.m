close all; clear all;
s= tf('s');

P = 1/(s^2+4*s+10);
C = ((1+.26*s)*(1+.26*s))/s


D = P/(1+P*C)

N = (P*C)/(1+P*C)



figure()

margin(N)
figure()
step(N)
figure

%sisotool(L/(L+1))


figure
%b= (s+100)/s;

nyquist(b)