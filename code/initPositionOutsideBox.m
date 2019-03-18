function [x_vec, y_vec] = initPositionOutsideBox(num_e, dim_x, dim_y, box_x, box_y)
%initPositionOutsideBox: Initializes x and y position vectors outside of
%box space

lim_x_low = (dim_x/2)-(box_x/2);
lim_x_high = (dim_x/2)+(box_x/2);

x_vec = zeros(1,num_e);
y_vec = zeros(1,num_e);

for i=1:num_e
    % Assign x and y position
    x_vec(i) = rand()*dim_x;
    y_vec(i) = rand()*dim_y;
    
    % Assign new x position while previous x is in the box limits
    % Only need to check x values for simplicity
    % TODO: Initialize position of electrons between the box limits
    while x_vec(i) > lim_x_low && x_vec(i) < lim_x_high
        x_vec(i) = rand().*dim_x;
    end
    
end

% figure(1);
% scatter(x_vec, y_vec)
% axis([0 dim_x 0 dim_y])
% rectangle('position', [lim_x_low 0 box_x box_y]);
% rectangle('position', [lim_x_low dim_y-box_y box_x box_y]);

end
