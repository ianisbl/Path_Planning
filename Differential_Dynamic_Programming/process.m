%% Process the continuous map into a discrete grid
clear all

%% Load the continuous elevation grid
load map.mat

%% Compute the discrete grid

Z_discrete = zeros(size(Z));
% change the holes into obstacles
Z_discrete = Z_discrete + (Z < 0);
% scale the map to get value between 0 and 1
p = 0.5;
Z_max = (1-p)*max(max(Z)); % lowest bound for an elevation being an obstacle
Z_discrete = Z_discrete + Z .* (Z>=0 & Z <= Z_max)/Z_max +  (Z>=0 & Z > Z_max);

surf(Z_discrete)
csvwrite('occ.csv', Z_discrete)