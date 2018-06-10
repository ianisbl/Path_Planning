clear all;

n=2;
m=2;
plot_init = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 25; % dimension of the grid
xx = 0:1:L;
yy = 0:1:L;
[XX, YY] = meshgrid(xx,yy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define terrain function %%%%%%%%%%%%%%%%%%%
% n_hills = 50;
n_holes = 2;

% add holes
holes = [5 7; 12 18];
radius = [3 3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define dynamics and cost-to-go
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [1;1];
x_goal = [L;L];


% Dynamics
f = @(x,y,u,v) [x+u; y+v];
fx = @(x,y,u,v) [1 0;0 1];
fu = @(x,y,u,v) [1 0;0 1];
fxx = @(x,y,u,v) [0 0;0 0];
%fxx = jacobian(fx,[x,y]);
%fuu = jacobian(fu,[u,v]);
%fux = jacobian(fu,[x,y]);

% Cost to go
R = eye(m);
l = @(x, y, u, v) [u;v]'*R*[u;v];
lx = @(x, y, u, v) [0; 0];
lxx = @(x, y, u, v) [0 0;0 0];
lu = @(x, y, u, v) [2*u; 2*v];
luu = @(x, y, u, v) [2 0; 0 2];
lux = @(x, y, u, v) [0 0;0 0];

% final cost
Qf = 50*eye(n);
lf = @(x,y) [x-x_goal(1);y-x_goal(2)]'*Qf*[x-x_goal(1);y-x_goal(2)];
lfx = @(x,y) [100*(x-x_goal(1)); 100*(y-x_goal(2))];
lfxx = @(x,y) [100 0;0 100];

% constraint functions
g = cell(n_holes,1);
gx = cell(n_holes,1);
gu = cell(n_holes,1);
for k=1:n_holes
    g{k} = @(x,y,u,v) radius(k)^2-norm([x,y]-holes(k,:),2)^2;
    gx{k} = @(x,y,u,v) [-2*(x-holes(k,1));-2*(y-holes(k,2))];
    gu{k} = @(x,y,u,v) [0;0];
end

save model.mat L x0 x_goal f fx fu l lx lxx lu luu ...
    lux lf lfx lfxx g gx gu holes radius %fxx fuu fux