close all; clear all; clc;


%define dipole (m)\
m = metal('lead')
d1 = dipoleCylindrical('Radius', 0.0006455,'Conductor',metal('copper'),'Load',lumpedElement('Impedance',50),'ClosedEnd',1,"Length",2.39)



%display antenna geometry



%define frequnecy range
f1 = 1e6;
f2 = 1000e6;
fLower = 54e6;  %bandwidth lower frequency
fUpper = 66e6;  %bandwidth higher frequency
FCenter = 60e6; %Center frequency
f_interval = .1e6; %frequency interval
frequencies = linspace(f1,f2);

%display the radiation pattern
%define dipole (100cm)
d2 = design(dipole,FCenter);

%calcuylate S parameters

figure()
efficiency(d1,frequencies)