clear all;

n=2;
m=2;
plot_init = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 25; % dimension of the grid
xx = linspace(0,L,100);
yy = linspace(0,L,100);
[XX, YY] = meshgrid(xx,yy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define terrain function %%%%%%%%%%%%%%%%%%%
n_hills = 1;
n_holes = 2;

hills = [13 11;23 21];

elev1 = @(x,y) 1e3/(3*L*sqrt(2*pi))*exp(-1/(3*L)*((x-hills(1,1))^2+(y-hills(1,2))^2));
elev1x = @(x,y) elev1(x,y)*[-2/(3*L)*(x-hills(1,1)); -2/(3*L)*(y-hills(1,2))];
elev1xx = @(x,y) elev1(x,y)*[-2/(3*L) 0; 0 -2/(3*L)];

elev2 = @(x,y) 10/(3*sqrt(2*pi))*exp(-1/(3)*((x-hills(2,1))^2+(y-hills(2,2))^2));
elev2x = @(x,y) elev2(x,y)*[-2/(3)*(x-hills(2,1)); -2/(3)*(y-hills(2,2))];
elev2xx = @(x,y) elev2(x,y)*[-2/(3) 0; 0 -2/(3)];

elev = @(x,y) elev1(x,y);% + elev2(x,y);
elevx = @(x,y) elev1x(x,y);% + elev2x(x,y);
elevxx = @(x,y) elev1xx(x,y);% + elev2xx(x,y);


Z = zeros(size(XX));
for i=1:length(XX)
    for j=1:length(XX)
        Z(i,j) = elev(XX(i,j),YY(i,j));
    end
end

% add holes
holes = [8 13; 13 17];
radius = [2 2];

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
l = @(x, y, u, v) [u;v]'*R*[u;v] + elev(x,y);
lx = @(x, y, u, v) [0; 0] + elevx(x,y);
lxx = @(x, y, u, v) [0 0;0 0] + elevxx(x,y);
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

save model.mat L x0 x_goal elev XX YY Z f fx fu l lx lxx lu luu ...
    lux lf lfx lfxx g gx gu holes radius %fxx fuu fux