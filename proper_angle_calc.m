% Script to calculate max viewing angle in aperture

% measured or determined values [mm]
d = 19.36;
c = 3.8;
b = 24.18;
a = [57.33 56.75 56.09];        % r g b
p = .5;
y_p = [4 8 12];                   % true for some pinhole points (cross)

% bfl values
z = a - b - d - p

% angle calcs
theta_red = atand(y_p/z(1))
theta_grn = atand(y_p/z(2))
theta_blu = atand(y_p/z(3))

