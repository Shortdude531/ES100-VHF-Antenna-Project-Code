
close all; clear all; clc; clf;

%nitinol conductivity S/m 1219512

%% define dipole (m)
d1 = dipoleCylindrical('Radius', 0.0006455,'Load',lumpedElement('Impedance',50),'ClosedEnd',1,"Length",2.39);
font=15;
figure()
show(d1)
set(findall(gcf,'-property','FontSize'),'FontSize',30)

%% % read S11 Experimental Data from CSv


%% define frequnecy range
f1 = 54e6;
f2 = 66e6;
fLower = 54e6;  %bandwidth lower frequency
fUpper = 66e6;  %bandwidth higher frequency
FCenter = 60e6; %Center frequency
f_interval = .1e6; %frequency interval
frequencies = linspace(fLower,fUpper,1001);


%% Calculate Simalted SParameters 

%% simulate the S11 to get the impedence efficiency of the antenna





%% simlate teh radiation effiency of the antenna 


%% calculate the polarization loss

%% calculate the efficeincy of the antenna 
EffTotal = linspace(1,1,1001);

%% Matlab Plot Parameters
 font = 13;
 linewidth = 1.25;
 %

%% JPL Mission Parameters for REASON

%REASON transmit upper and lower bound (Watts)
    P_t_u = 10;

%REASON boresight gain upper and lower bound (dBi)

    G_B = 10;

%REASON Side lobe gain upper and lower bound (dBi)
    G_S = 10^(-15/20);
 

%REASON polarization loss 
    L_pol = .7;

%REASON Coherence loss 
    L_c = .5;

%REASON duty cycle
    d = .1;

%REASON allowable observation time for calibration (s)
    Obs = 20;
    
    Obs_a = linspace(0,600,length(fdata));
    Obs_a_s =Obs_a*d;
%Receive noise temperture upper and lower bound
    T_l = 2000;
    T_u = 5800;

%frequency badwidth upper and lower bound
    f_l = 54e6;
    f_c = 60e6;
    f_u = 66e6;
    B = f_u-f_l;

%min allowable seperation of CaliPer from Clipper
R = 75e3;

%% CaliPer Parameters
%Receive gain for isotropic antenna
G_r_max = 1.93;

%%S11 Receive atenna effiency at a certain frequency
%T = readtable('S11for2481000','NumHeaderLines',1);

%number of atenna elements    
N_elm = 1;

% T_obs
T_obs = d*Obs;
T_obs_a = Obs_a*d;

%% Universal constants

%Boltzman Constant
    k = .380649e-23;
    c = 3e8;
%% Wavelenth upper and lower bound
    lambda_l = c/f_l;
    lambda_u = c/f_u;
    lambda_c = c/f_c;

%k and h for 60e6 frequency

%REASON Coherence loss 
L_c = .5;
k_a= (2*pi)/5;
h=L_c/2;
%
theta = linspace(pi/6,((5*pi)/6),explen);
LTheta = length(theta);
for i = 1:LTheta
    g(i) = abs(((cos(k_a*h*cos(theta(i)))-cos(k_a*h))/sin(theta(i))));
end
%normalize g
g = g/max(g);

nadir = 75e3;
distance_r = nadir./abs(sin(theta));

freLen = length(frequencies);
lambda = c./frequencies;

%%SNR
P_iso_sidelobe_obs = zeros(freLen,freLen);
for i = 1:freLen
    for j = 1:freLen
        P_iso_sidelobe_obs(i,j) = G_r_max*EffTotal(i)*P_t_u*G_S*((lambda(i))^2/(4*pi))*N_elm*(1/(4*pi*distance_r(j)^2))*L_pol;
    end
end

%%Signal to Noise Power Ratio for tranmission (dB)
SNR_sidelobe_obs = zeros(freLen,freLen);
for i = 1:freLen
    for j = 1:freLen
        SNR_sidelobe_obs(i,j) = 10*log((2*P_iso_sidelobe_obs(i,j)*T_obs*L_c*(1/(k*T_u))));
    end
end


font=15;
figure()

meshnum=20;


theta_degrees = linspace(1,120,freLen);

x = frequencies(1:meshnum:num);
y = theta_degrees(1:meshnum:num);


[X, Y] = meshgrid(x, y);
surf(x, y, SNR_sidelobe_obs(1:meshnum:num,1:meshnum:num)','FaceAlpha',1)

font=30;


title("Antenna SNR of Reason Side Lobes in terms of Frequency & Distance",'FontSize',font)
xlabel("Frequency (Hz)",'FontSize',font)
ylabel("Angle (Degrees)",'FontSize',font)
zlabel("SNR (dB)",'FontSize',font)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

set(findall(gcf,'-property','FontSize'),'FontSize',font)
%assuming a max gain of 2.11, a distance of 1000km, and a obervation time of 15 minutes 


%% 
figure()


%S11EXP = (1-10.^(S11data/10));
%S11THE = (1-10.^(YData/10));
plot(frequencies,YData,frequencies,S11data,'LineWidth', 1.5)
ylim([-50 0])
xlim([54e6 66e6])
xlabel("Frequency Hz")
ylabel("S11 Reflected Power DB")
title("S11 Reflected Power Simulated vs. Experimental")
set(findall(gcf,'-property','FontSize'),'FontSize',13)
% %create a optimized 3 comonent Matchign network
% mnobj = matchingnetwork('LoadImpedance', S ,'Components',3,'CenterFrequency',60e6, 'Bandwidth',6e6);
% hline = rfplot(mnobj,frequencies,1);
% D=get(gca,'Children'); %get the handle of the line object
% XData=get(D(1,1),'XData'); %get the x data


