
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

%% simlate teh radiation effiency of the antenna 
figure()
efficiency(d1,frequencies)


D=get(gca,'Children'); %get the handle of the line object
XData1=get(D,'XData'); %get the x data
YData1=get(D,'YData'); %get the y data
Data1=[XData' YData']; %join the x and y data on one array nx2



save('SimS11.mat','S','XData1','YData1','XData',"YData")

