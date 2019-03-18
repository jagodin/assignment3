function [vx,vy] = newBoltDist()
%newBoltDist: Generates new x and y velocity component based on
%Maxwell-Boltzmann distribution
global k, global T, global mn;

xf = (normrnd(0, sqrt(k*T/mn)));
yf = (normrnd(0, sqrt(k*T/mn)));

dist = (xf^2 + yf^2)^0.5;
            
angle = rand*2*pi;
vx = dist.*cos(angle);
vy = dist.*sin(angle);
end

