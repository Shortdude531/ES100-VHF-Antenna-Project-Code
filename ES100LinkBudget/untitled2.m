
close all; clear all; clc; clf;
%% define dipole (m)
d1 = dipoleCylindrical('Radius', 0.0006455,'Conductor',metal('Copper'),'Load',lumpedElement('Impedance',50),'ClosedEnd',1,"Length",2.39);

%% % read S11 Experimental Data from CSv
S11 = readtable('ES100 Data/MyFilePrefixJV2','NumHeaderLines',0);
datasize = size(S11);

%Delete Keysight Fieldfox Header
S11(1:17,:) = [];
S11(size(S11),:) = [];

S11 = S11{:,:};

fdata = S11(:,1);
S11data = S11(:,2);



%% define frequnecy range
f1 = 50e6;
f2 = 70e6;
fLower = 54e6;  %bandwidth lower frequency
fUpper = 66e6;  %bandwidth higher frequency
FCenter = 60e6; %Center frequency
f_interval = .1e6; %frequency interval
frequencies = fdata;

%calcuylate S parameters
S = sparameters(d1,frequencies);

%% simulate the S11 to get the impedence efficiency of the antenna
figure()
rfplot(S)

D=get(gca,'Children'); %get the handle of the line object
XData=get(D,'XData'); %get the x data
YData=get(D,'YData'); %get the y data
Data=[XData' YData']; %join the x and y data on one array nx2

%calculate impedence efficiency  
ImpEff = (1-10.^(YData/10));

%% 
figure()


%S11EXP = (1-10.^(S11data/10));
%S11THE = (1-10.^(YData/10));
plot(frequencies,YData,frequencies,S11data)
ylim([-50 0])
xlim([54e6 66e6])
%% simlate teh radiation effiency of the antenna 
figure()
efficiency(d1,frequencies)


D=get(gca,'Children'); %get the handle of the line object
XData1=get(D,'XData'); %get the x data
YData1=get(D,'YData'); %get the y data
Data1=[XData' YData']; %join the x and y data on one array nx2

%% calculate the polarization loss

%% calculate the efficeincy of the antenna 
EffTotal = ImpEff.*YData1;

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
    Obs = .25*60^2;
    
    Obs_a = linspace(0,8);
    Obs_a_s =Obs_a*60^2*d;
%Receive noise temperture upper and lower bound
    T_l = 2000;
    T_u = 5800;

%frequency badwidth upper and lower bound
    f_l = 54e6;
    f_c = 60e6;
    f_u = 66e6;
    B = f_u-f_l;

%min allowable seperation of CaliPer from Clipper
R = 1000e3;

%% CaliPer Parameters
%Receive gain for isotropic antenna
G_r_max = 2.11;

%%S11 Receive atenna effiency at a certain frequency
T = readtable('S11for2481000','NumHeaderLines',1);

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

    
%% calculate the SNR for the dipole antenna in terms of frequency at maximum dipole gain for boresight and side lobe of reason. 

%% SNR vs. distance

freLen = length(frequencies)
P_iso_sidelobe = zeros(freLen,1);
lambda = c./frequencies;



P_iso_sidelobe = zeros(freLen,1);
for i = 1:freLen
    P_iso_sidelobe_distance(i) = EffTotal(i)*P_t_u*G_S*((lambda(i))^2/(4*pi))*N_elm*(1/(4*pi*R^2))*L_pol;
    end

%%Signal to Noise Power Ratio for tranmission (dB)
for i = 1:freLen
    SNR_sidelobe_distance(i) = 10*log((2*P_iso_sidelobe_distance(i)*T_obs*L_c*(1/(k*T_u))));
end


P_iso_boresight = zeros(freLen,1);
for i = 1:freLen
   P_iso_boresight_distance(i) = EffTotal(i)*P_t_u*G_B*((lambda(i))^2/(4*pi))*N_elm*(1/(4*pi*R^2))*L_pol;
end

%%Signal to Noise Power Ratio for tranmission (dB)
for i = 1:freLen
    SNR_boresight_distance(i) = 10*log((2.11*P_iso_boresight_distance(i)*T_obs*L_c*(1/(k*T_u))));
end
    

figure()
plot(frequencies,SNR_sidelobe_distance,frequencies,SNR_boresight_distance,'LineWidth', linewidth)
title("Antenna SNR versus Frequency",'FontSize',font)
xlabel("Frequency (MHz)",'FontSize',font)
ylabel("SNR (dB)",'FontSize',font)
legend("Side Lobe -15dBi","Boresight 10dBi",'FontSize',font)
set(findall(gcf,'-property','FontSize'),'FontSize',font)
%assuming a max gain of 2.11, a distance of 1000km, and a obervation time of 15 minutes 

%% plot the SNR in terms of pointing angle and frequency

%% Pointing Dirrection and Recieve Gain


%k and h for 60e6 frequency
k_a= (2*pi)/5;
h=L_c/2;

%
theta = linspace(0,2*pi,100);
LTheta = length(theta);
for i = 1:LTheta
    g(i) = abs(((cos(k_a*h*cos(theta(i)))-cos(k_a*h))/sin(theta(i))));
end
%normalize g
g = g/max(g);




P_iso_sidelobe_point = zeros(freLen,1);
for i = 1:freLen
    for j = 1:LTheta
        P_iso_sidelobe_point(i,j) = g(j)*EffTotal(i)*P_t_u*G_S*((lambda(i))^2/(4*pi))*N_elm*(1/(4*pi*R^2))*L_pol;
    end
end

% Signal to Noise Power Ratio for tranmission (dB)
for i = 1:freLen
    for j = 1:LTheta
        SNR_sidelobe_point(i,j) = 10*log((2*P_iso_sidelobe_point(i,j)*T_obs*L_c*(1/(k*T_u))));
    end
end


    
figure()
x = frequencies;
y = theta;
[X, Y] = meshgrid(x, y);
surf(x, y, SNR_sidelobe_point','FaceAlpha',1)
colorbar
camlight('headlight')
camlight('left')


title("Antenna SNR of Reason Side Lobes in terms of Frequency & Pointing Angle",'FontSize',font)
xlabel("Frequncy (MHz)",'FontSize',font)
ylabel("\theta bounded from 0 to 2\pi (rad)",'FontSize',font)
zlabel("SNR (dB)",'FontSize',font)
set(findall(gcf,'-property','FontSize'),'FontSize',font)


%% 3D plot of SNR vs frequncy & distance
R_arr = linspace(1000e3,1000e6,100);
Rlen = length(R_arr)






P_iso_sidelobe_distance = zeros(freLen,Rlen);
for i = 1:freLen
    for j = 1:Rlen
        P_iso_sidelobe_distance(i,j) = EffTotal(i)*P_t_u*G_S*((lambda(i))^2/(4*pi))*N_elm*(1/(4*pi*R_arr(j)^2))*L_pol;
    end
end

%%Signal to Noise Power Ratio for tranmission (dB)
SNR_sidelobe_distance = zeros(freLen,Rlen);
for i = 1:freLen
    for j = 1:Rlen
        SNR_sidelobe_distance(i,j) = 10*log((2*P_iso_sidelobe_distance(i,j)*T_obs*L_c*(1/(k*T_u))));
    end
end


    
figure()
x = frequencies;
y = R_arr;
[X, Y] = meshgrid(x, y);
surf(x, y, SNR_sidelobe_distance','FaceAlpha',1)
colorbar



title("Antenna SNR of Reason Side Lobes in terms of Frequency & Distance",'FontSize',font)
xlabel("Frequncy (MHz)",'FontSize',font)
ylabel("Distance from REASON (km)",'FontSize',font)
zlabel("SNR (dB)",'FontSize',font)
set(findall(gcf,'-property','FontSize'),'FontSize',font)


%% 3D plot of SNR vs frequncy & Observation time
Obs_a = linspace(0,8);
Obs_a_s =Obs_a*60^2*d;



obsLen = length(Obs_a_s);






P_iso_sidelobe_obs = zeros(freLen,obsLen);
for i = 1:freLen

        P_iso_sidelobe_obs(i) = EffTotal(i)*P_t_u*G_S*((lambda(i))^2/(4*pi))*N_elm*(1/(4*pi*R^2))*L_pol;

end

%%Signal to Noise Power Ratio for tranmission (dB)
SNR_sidelobe_obs = zeros(freLen,obsLen);
for i = 1:freLen
    for j = 1:obsLen
        SNR_sidelobe_obs(i,j) = 10*log((2*P_iso_sidelobe_obs(i)*Obs_a_s(j)*L_c*(1/(k*T_u))));
    end
end


    
figure()
x = frequencies;
y = Obs_a;
[X, Y] = meshgrid(x, y);
surf(x, y, SNR_sidelobe_obs','FaceAlpha',1)
colorbar
camlight('headlight')
camlight('left')


title("Antenna SNR of Reason Side Lobes in terms of Frequency & Observation Time",'FontSize',font)
xlabel("Frequncy (MHz)",'FontSize',font)
ylabel("Observation Time (Hr)",'FontSize',font)
zlabel("SNR (dB)",'FontSize',font)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

%% Pointing Dirrection and Recieve Gain


%k and h for 60e6 frequency
k_a= (2*pi)/5;
h=L_c/2;

%
theta = linspace(0,2*pi,100);
LTheta = length(theta);
for i = 1:LTheta
    g(i) = abs(((cos(k_a*h*cos(theta(i)))-cos(k_a*h))/sin(theta(i))));
end
%normalize g
g = g/max(g);




P_iso_sidelobe_point = zeros(freLen,1);
for i = 1:freLen
    for j = 1:LTheta
        P_iso_sidelobe_point(i,j) = g(j)*EffTotal(i)*P_t_u*G_B*((lambda(i))^2/(4*pi))*N_elm*(1/(4*pi*R^2))*L_pol;
    end
end

% Signal to Noise Power Ratio for tranmission (dB)
for i = 1:freLen
    for j = 1:LTheta
        SNR_sidelobe_point(i,j) = 10*log((2*P_iso_sidelobe_point(i,j)*T_obs*L_c*(1/(k*T_u))));
    end
end


    
figure()
x = frequencies;
y = theta;
[X, Y] = meshgrid(x, y);
surf(x, y, SNR_sidelobe_point','FaceAlpha',1)
colormap turbo;
colorbar;
camlight('headlight')
camlight('left')



title("Antenna SNR of Reason Boresight in terms of Frequency & Pointing Angle",'FontSize',font)
xlabel("Frequncy (MHz)",'FontSize',font)
ylabel("\theta bounded from 0 to 2\pi (rad)",'FontSize',font)
zlabel("SNR (dB)",'FontSize',font)
set(findall(gcf,'-property','FontSize'),'FontSize',font)


%% 3D plot of SNR vs frequncy & distance
R_arr = linspace(1000e3,1000e6,100);
Rlen = length(R_arr)






P_iso_sidelobe_distance = zeros(freLen,Rlen);
for i = 1:freLen
    for j = 1:Rlen
        P_iso_sidelobe_distance(i,j) = EffTotal(i)*P_t_u*G_B*((lambda(i))^2/(4*pi))*N_elm*(1/(4*pi*R_arr(j)^2))*L_pol;
    end
end

%%Signal to Noise Power Ratio for tranmission (dB)
SNR_sidelobe_distance = zeros(freLen,Rlen);
for i = 1:freLen
    for j = 1:Rlen
        SNR_sidelobe_distance(i,j) = 10*log((2*P_iso_sidelobe_distance(i,j)*T_obs*L_c*(1/(k*T_u))));
    end
end


    
figure()
x = frequencies;
y = R_arr;
[X, Y] = meshgrid(x, y);
surf(x, y, SNR_sidelobe_distance','FaceAlpha',1)
colormap turbo;
colorbar;




title("Antenna SNR of Reason Boresight in terms of Frequency & Distance",'FontSize',font)
xlabel("Frequncy (MHz)",'FontSize',font)
ylabel("Distance from REASON (km)",'FontSize',font)
zlabel("SNR (dB)",'FontSize',font)
set(findall(gcf,'-property','FontSize'),'FontSize',font)


%% 3D plot of SNR vs frequncy & Observation time
Obs_a = linspace(0,8);
Obs_a_s =Obs_a*60^2*d;



obsLen = length(Obs_a_s);






P_iso_sidelobe_obs = zeros(freLen,obsLen);
for i = 1:freLen

        P_iso_sidelobe_obs(i) = EffTotal(i)*P_t_u*G_B*((lambda(i))^2/(4*pi))*N_elm*(1/(4*pi*R^2))*L_pol;

end

%%Signal to Noise Power Ratio for tranmission (dB)
SNR_sidelobe_obs = zeros(freLen,obsLen);
for i = 1:freLen
    for j = 1:obsLen
        SNR_sidelobe_obs(i,j) = 10*log((2*P_iso_sidelobe_obs(i)*Obs_a_s(j)*L_c*(1/(k*T_u))));
    end
end


    
figure()
x = frequencies;
y = Obs_a;
[X, Y] = meshgrid(x, y);
surf(x, y, SNR_sidelobe_obs','FaceAlpha',1)
colormap turbo;
colorbar;
camlight('headlight')
camlight('left')


title("Antenna SNR of Reason Boresight in terms of Frequency & Observation Time",'FontSize',font)
xlabel("Frequncy (MHz)",'FontSize',font)
ylabel("Observation Time (Hr)",'FontSize',font)
zlabel("SNR (dB)",'FontSize',font)
set(findall(gcf,'-property','FontSize'),'FontSize',font)


 for i =4:10
saveas(i,sprintf('SmallDipole%d.png',i));
 
 end

