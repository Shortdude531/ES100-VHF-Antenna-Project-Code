close all; clear all; clc;
f1 = 50e6;
f2 = 70e6;
fLower = 54e6;  %bandwidth lower frequency
fUpper = 66e6;  %bandwidth higher frequency
FCenter = 60e6; %Center frequency
f_interval = .1e6; %frequency interval
frequencies = linspace(f1,f2);


% d1 = dipoleCylindrical('Radius', 0.0006455,'Conductor',metal('Copper'),'Load',lumpedElement('Impedance',50),'ClosedEnd',1,"Length",2.39);

t1 = design(dipoleCrossed,FCenter);
d1 = design(dipoleCylindrical,FCenter);
df1 = design(dipoleFolded,FCenter);

c1 = design(patchMicrostripCircular,FCenter);
dd1 = design(dipole,FCenter);

y1 = design(yagiUda,FCenter);

t1 = design(dipoleCrossed,FCenter);
he = design(helix,FCenter);
h1 = design(hornConical,FCenter);

m = design(monopoleCylindrical,FCenter);


font = 15;


figure()
show(d1)

set(findall(gcf,'-property','FontSize'),'FontSize',font)
figure()
show(df1)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
show(dd1)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
show(m)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
show(t1)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
show(c1)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
show(y1)
set(findall(gcf,'-property','FontSize'),'FontSize',font)


figure()
show(he)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
show(h1)
set(findall(gcf,'-property','FontSize'),'FontSize',font)





figure()
pattern(d1,FCenter)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
pattern(df1,FCenter)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
pattern(dd1,FCenter)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
pattern(m,FCenter)
set(findall(gcf,'-property','FontSize'),'FontSize',font)


figure()
pattern(t1,FCenter)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
pattern(c1,FCenter)

set(findall(gcf,'-property','FontSize'),'FontSize',font)


figure()
pattern(y1,FCenter)
set(findall(gcf,'-property','FontSize'),'FontSize',font)



figure()
pattern(he,FCenter)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
pattern(h1,FCenter)
set(findall(gcf,'-property','FontSize'),'FontSize',font)



figure()
S1= sparameters(d1,frequencies,50);
rfplot(S1)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
S2=sparameters(df1,frequencies,50);
rfplot(S2)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
S4=sparameters(dd1,frequencies,50);
rfplot(S4)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
S6=sparameters(t1,frequencies,50);
rfplot(S6,1,1)
set(findall(gcf,'-property','FontSize'),'FontSize',font)


figure()
S9=sparameters(m,frequencies,50);
rfplot(S9)
set(findall(gcf,'-property','FontSize'),'FontSize',font)



figure()
S3=sparameters(c1,frequencies,50);
rfplot(S3)
set(findall(gcf,'-property','FontSize'),'FontSize',font)


figure()
S5=sparameters(y1,frequencies,50);
rfplot(S5)
set(findall(gcf,'-property','FontSize'),'FontSize',font)


figure()
S7=sparameters(he,frequencies,50);
rfplot(S7)
set(findall(gcf,'-property','FontSize'),'FontSize',font)

figure()
S8=sparameters(h1,frequencies,50);
rfplot(S8)
set(findall(gcf,'-property','FontSize'),'FontSize',font)






 for i =1:27
saveas(i,sprintf('%d.png',i));
 
 end

