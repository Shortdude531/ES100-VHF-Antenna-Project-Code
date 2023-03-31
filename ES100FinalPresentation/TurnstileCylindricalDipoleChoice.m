close all; clear all; clc;
f1 = 50e6;
f2 = 70e6;
fLower = 54e6;  %bandwidth lower frequency
fUpper = 66e6;  %bandwidth higher frequency
FCenter = 60e6; %Center frequency
f_interval = .1e6; %frequency interval
frequencies = linspace(f1,f2,500);


% d1 = dipoleCylindrical('Radius', 0.0006455,'Conductor',metal('Copper'),'Load',lumpedElement('Impedance',50),'ClosedEnd',1,"Length",2.39);
p = PatternPlotOptions;
t1 = design(dipoleCrossed,FCenter);
d1 = design(dipoleCylindrical,FCenter);
%% 
figure()
p.Transparency = 0.5;
pattern(d1,FCenter,'patternOptions',p)
figure()
EHfields(t1,FCenter)
figure()
EHfields(d1,FCenter)

p.Transparency = 0.5;
figure()
pattern(t1,FCenter,'patternOptions',p)
set(findall(gcf,'-property','FontSize'),'FontSize',30)
figure()

set(findall(gcf,'-property','FontSize'),'FontSize',30)
 set(findall(gcf,'-property','FontSize'),'FontSize',30)


 figure()
 patternElevation(t1,FCenter)

 set(findall(gcf,'-property','FontSize'),'FontSize',30)

figure()
   patternElevation(d1,FCenter)

 set(findall(gcf,'-property','FontSize'),'FontSize',30)

  figure()
  patternAzimuth(t1,FCenter)
   set(findall(gcf,'-property','FontSize'),'FontSize',30)
   figure()
   patternAzimuth(d1,FCenter)


 set(findall(gcf,'-property','FontSize'),'FontSize',30)


 S_t1 = sparameters(t1,frequencies,50);
 figure()
rfplot(S_t1)
set(findall(gcf,'-property','FontSize'),'FontSize',30)

S_d1 = sparameters(d1,frequencies,50);
figure()
rfplot(S_d1)
set(findall(gcf,'-property','FontSize'),'FontSize',30)

 p = smithplot(S_d1,1,1, 'GridType','ZY');
