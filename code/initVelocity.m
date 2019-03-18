function [vx_vec, vy_vec] = initVelocity(num_e, Vth)
%initPosition: Initializes x and y velocity vectors

angle = rand(1,num_e)*2*pi;

vx_vec = Vth*cos(angle);
vy_vec = Vth*sin(angle);

end
