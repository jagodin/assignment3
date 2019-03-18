function [x_vec, y_vec] = initPosition(num_e, dim_x, dim_y)
%initPosition: Initializes x and y position vectors

x_vec = rand(1,num_e)*dim_x;
y_vec = rand(1,num_e)*dim_y;

end

