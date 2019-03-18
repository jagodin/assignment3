function [vx_vec,vy_vec] = initBoltDist(num_e)
%initBoltDist: Initializes x and y velocity vectors in a Maxwell-Boltzman
%distribution taken from a Gaussian distribution.

global k, global T, global mn;

dist = ((sqrt(k*T/mn)*randn(1,num_e)).^2 + (sqrt(k*T/mn)*randn(1,num_e)).^2).^0.5;

angle = rand(1,num_e)*2*pi;
vx_vec = dist.*cos(angle);
vy_vec = dist.*sin(angle);
end

