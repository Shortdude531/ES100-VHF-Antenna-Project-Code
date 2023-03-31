
close all; clear all; clc;

%% Matlab Plot Parameters
 font = 14;
 linewidth = 1.25;
 %

%% JPL Mission Parameters for REASON

%REASON transmit upper and lower bound (Watts)
    P_t_l = 8;
    P_t_u = 10;

%REASON boresight gain upper and lower bound (dBi)
    G_0_l = 9;
    G_0_u = 10;

%REASON Side lobe gain upper and lower bound (dBi)
    G_ts_l = 10^(-15/20);
    G_ts_u = 10^(0/20);

%REASON polarization loss 
    L_pol = .7;

%REASON Coherence loss 
    L_c = .5;

%REASON duty cycle
    d = .1;

%REASON allowable observation time for calibration (s)
    Obs = .25*60^2;
    
    Obs_a = .01:.01:8;
    Obs_a_s =Obs_a*60^2;
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
%Receive gain
G_r_max = 2.15;

%%S11 Receive atenna effiency at a certain frequency
T = readtable('S11for2481000','NumHeaderLines',1);



%% Antenna radiation efficeincy (reasonable placeholder efficiency constant)
episilon_r = 1;

% L/50

%% Calulate Effective Area
% Get the Antenna width and height
    w_a = 2*table2array(T(:,2));
    h_a = table2array(T(:,3));
    length = height(w_a);

% convert S11 dB for 54,60,66MHz to decimal
f_a_l = zeros(length,1);
f_a_c = zeros(length,1);
f_a_h = zeros(length,1);

for i = 1:height(T(:,3))

    f_a_l(i) = 1-10^(table2array(T(i,4))/10);
    f_a_c(i) = 1-10^(table2array(T(i,6))/10);
    f_a_h(i) = 1-10^(table2array(T(i,5))/10);
end

%find A_eff
    A_eff = zeros(length);
    for i = 1:length
        A_eff(i) = G_r_max*f_a_h(i)*episilon_r;
    end



%simulated atenna efficenty from .01 to 100%
eff = .01:0.01:1;
leneff = width(eff);
%number of atenna elements    
N_elm = 1;

%% T_obs
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
%% Power received by one RF emission from REASON upper and lower bound side lobe
    P_r_l_e = zeros(length,leneff);
    
    for i = 1:length
        for j = 1:leneff
            P_r_l_e(i,j) = P_t_l*G_ts_l*(lambda_u)^2/(4*pi)*A_eff(i)*N_elm*(1/(4*pi*R^2))*L_pol*eff(j);
        end
    end

%%Signal to Noise Power Ratio for tranmission (dB)
    SNR = zeros(length, leneff);
    for i = 1:length
        for j = 1:leneff
            SNR(i,j) = 10*log((2*P_r_l_e(i,j)*T_obs*L_c*(1/(k*T_u))));
        end
    end
figure()
plot(eff,SNR,'LineWidth', linewidth)
title("SNR vs. Antenna efficiency -15dBi REASON signal 3000km and 1 hours",'FontSize',font)
xlabel("Antenna efficiency",'FontSize',font)
ylabel("SNR (dB)",'FontSize',font)
legend("4mm","48mm","100mm",'FontSize',font)

%%REASON Boresight power as a function of antenna diameter

for i = 1:length  
            P_r_l_bor(i) = P_t_l*G_0_u*(lambda_l)^2/(4*pi)*A_eff(i)*N_elm*(1/(4*pi*R^2))*L_pol*f_a_l(i);
   
end


for i = 1:length  
            P_r_u_bor(i) = P_t_l*G_0_u*(lambda_u)^2/(4*pi)*A_eff(i)*N_elm*(1/(4*pi*R^2))*L_pol*f_a_h(i);
   
end

for i = 1:length  
            P_r_c_bor(i) = P_t_l*G_0_u*(lambda_c)^2/(4*pi)*A_eff(i)*N_elm*(1/(4*pi*R^2))*L_pol*f_a_c(i);
   
end

figure()
plot(w_a,P_r_l_bor',w_a,P_r_c_bor',w_a,P_r_u_bor','LineWidth', linewidth)
title("Antenna Boresight Power at 3000km at 1 hours observation time",'FontSize',font)
xlabel("Antenna Width (m)",'FontSize',font)
ylabel("Power (W)",'FontSize',font)
legend("54e6 Mhz","60e6 Mhz","66e6 Mhz",'FontSize',font)

%%W SNR observation time, lower bound efficiency, width
   f_a_l_len = width(f_a_l);
   T_obs_a_len = width(Obs_a_s);
 
    
                for i = 1:length
                P_r_l_e4(i) = P_t_l*G_ts_l*(lambda_l)^2/(4*pi)*A_eff(i)*N_elm*(1/(4*pi*R^2))*L_pol;
                end

                for i = 1:length
                    for j = 1:T_obs_a_len
                        SNR1(i,j) = 20*log((2*P_r_l_e4(i)*Obs_a_s(j)*L_c*(1/(k*T_u))));
                    end
                end
figure()
 plot(Obs_a,SNR1,'LineWidth', linewidth)
 title("Antenna SNR vs observation time at 3000km -15dBi REASON signal",'FontSize',font)
 xlabel("Observation time (Hr)",'FontSize',font)
 ylabel("SNR (dB)",'FontSize',font)
legend("4mm","48mm","100mm",'FontSize',font)

%%SNR vs. distance
R_a = 100:100:100e3;
R_a_m = R_a*1e3;
R_a_len = width(R_a);
           
                for i = 1:length
                    for j = 1:R_a_len
                     P_r_l_e(i,j) = P_t_l*G_ts_l*(lambda_l)^2/(4*pi)*A_eff(i)*N_elm*(1/(4*pi*R_a_m(j)^2))*L_pol;
                    end
                end
                for i = 1:length
                    for j = 1:R_a_len
                        SNR5(i,j) = 20*log((2*P_r_l_e(i,j)*Obs*L_c*(1/(k*T_u))));
                    end
                end
    

figure()
plot(R_a,SNR5,'LineWidth', linewidth)
title("Antenna SNR vs distance -15dBi REASON signal 1 hour",'FontSize',font)
xlabel("Distance (km)",'FontSize',font)
ylabel("SNR (dB)",'FontSize',font)
legend("4mm","48mm","100mm",'FontSize',font)


%% Pointing Dirrection and Recieve Gain


%k and h for 60e6 frequency
k_a= (2*pi)/lambda_c;
h=L_c/2;

%
theta = linspace(0,2*pi,1000);
w_theta = width(theta);
for i = 1:w_theta
    g(i) = abs(((cos(k_a*h*cos(theta(i)))-cos(k_a*h))/sin(theta(i))));
end
%normalize g
g = g/max(g);

figure()
plot(theta,g)

%calculate power while changing pointing dirrection

for j = 1:length
    for i = 1:w_theta
        A_eff2(j,i) = (2.11*g(i)'*f_a_l(j)*episilon_r);
    end
end

 P_r_l_e2 = zeros(length,w_theta);
 for i = 1:length
    for j = 1:w_theta
        P_r_l_e2(i,j) = P_t_l*G_ts_l*(lambda_l)^2/(4*pi)*A_eff2(i,j)*N_elm*(1/(4*pi*R^2))*L_pol;
    end
 end

%% calculate SNR while changing pointing dirrection
    SNR2 = zeros(length,w_theta);
    for i =1:length
        for j = 1:w_theta
            SNR2(i,j) = 10*log((2*P_r_l_e2(i,j)*T_obs*L_c*(1/(k*T_u))));
        end
    end
figure()
plot(theta,SNR2,'LineWidth', linewidth)
title("SNR of REASON Side Lobe",'FontSize',font)
xlabel("\theta bounded from 0 to 2\pi (rad)",'FontSize',font)
ylabel("SNR (dB)",'FontSize',font)
legend("4mm","48mm","100mm",'FontSize',font)


%% Worst case SNR observation time, lower bound efficiency, width


A_eff3 = zeros(w_theta,1);

    for i = 1:w_theta
        A_eff3(i) = (2.11*g(i)'*.2*episilon_r);
    end

P_r_l_e3 = zeros(length,w_theta);

    for i = 1:w_theta
        P_r_l_e3(i) = P_t_l*G_ts_l*(lambda_l)^2/(4*pi)*A_eff3(i)*N_elm*(1/(4*pi*R^2))*L_pol;
    end

T_obs_worst = d*5*60;
%%calculate SNR while changing pointing dirrection
    SNR3 = zeros(w_theta,1);

        for i = 1:w_theta
            SNR3(i) = 10*log((2*P_r_l_e3(i)*T_obs_worst*L_c*(1/(k*T_u))));
        end

figure()



SNR_worst_max =  max(SNR3);
percet_cutoff = SNR_worst_max*.6;
SNR3_cutoff =SNR3;

indices = find(abs(SNR)>percet_cutoff);
    %SNR3_cutoff(indices) = NaN;
plot(theta,SNR3,theta, SNR3_cutoff,'LineWidth', linewidth)
title("Worth case Antenna SNR based on Pointing Direction",'FontSize',font)
xlabel("\theta bounded from 0 to 2\pi (rad)",'FontSize',font)
ylabel("SNR (dB)",'FontSize',font)

%code from https://www.mathworks.com/matlabcentral/answers/72396-x-value-on-y-max



acceptable_tilt_angle = (pi/2-.44)*(180/pi)

