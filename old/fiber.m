clc; clear;

a = 105e-6; % [m], fiber core diameter
bd = 0.8e-3; % beam diamter [m]
f = 1.85e-3; % lens focal length [m]
D = 2.8e-3; % lens diamter [m]
NA = bd / (2 * f); % numerical aperture [-]
theta = asin(NA);  % acceptance ange [rad]
da = a / f; % divergence angle [rad]

L = 15e-2; % [m]
r = L * da; % [m]