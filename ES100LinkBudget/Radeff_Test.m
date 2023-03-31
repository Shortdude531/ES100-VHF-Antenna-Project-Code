


d1 = dipoleCylindrical('Length',2.35,'Radius', 0.0006455,'Conductor',metal('Copper'))



%sobj = sparameters(d1,60e6,50);

f1 = 50e6;
f2 = 70e6;
fLower = 54e6;
fUpper = 66e6;
FCenter = 60e6;
f_interval = .1e6;
frequencies = linspace(f1,.1e6,f2);

efficiency(d1,frequencies)