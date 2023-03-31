%k and h for 60e6 frequency

%REASON Coherence loss 
L_c = .5;
k_a= (2*pi)/5;
h=L_c/2;
%
theta = linspace(pi/6,((5*pi)/6),100);
LTheta = length(theta);
for i = 1:LTheta
    g(i) = abs(((cos(k_a*h*cos(theta(i)))-cos(k_a*h))/sin(theta(i))));
end
%normalize g
g = g/max(g);

nadir = 75;
distance_r = nadir./abs(sin(theta));

figure()
plot(theta,g)