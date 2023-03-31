close all; clear all; clc; clf;
%nitinol conductivity S/m 1219512
m = metal("Name",'Nitinol','Conductivity',1219512,'Thickness', 2e-07);
%% define dipole (m)
d1 = dipoleCylindrical('Radius', 1.2e-3,'Conductor',m,'Load',lumpedElement('Impedance',50),'ClosedEnd',1,"Length",2.39);

figure()
show(d1)
 set(findall(gcf,'-property','FontSize'),'FontSize',30)
figure()
 efficiency(d1,linspace(54e6,66e6,10))

 set(findall(gcf,'-property','FontSize'),'FontSize',30)