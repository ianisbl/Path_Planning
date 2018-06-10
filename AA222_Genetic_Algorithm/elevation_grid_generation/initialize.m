clear all;

n=2;
m=2;
plot_init = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 50; % dimension of the grid
xx = 0:1:L;
yy = 0:1:L;
[XX, YY] = meshgrid(xx,yy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define terrain function %%%%%%%%%%%%%%%%%%%
n_hills = 50;
n_holes = 3;

% base terrain is a flat hill
elev = @(x1,x2) 750*mvnpdf([x1,x2],[L/2,L/2],L*eye(2));

% add some steeper hills
hills = L*rand(n_hills,2);
stds = exprnd(5*ones(n_hills));
for i=1:n_hills
    elev = @(x1,x2) elev(x1,x2) + 10*mvnpdf([x1,x2],hills(i,:),stds(i)*eye(2));
end

% add holes
holes = L*rand(n_holes,2);
radius = normrnd(3, 0.2, [n_holes,1]);

% display terrain
Z = zeros(size(XX));
for i=1:size(XX,1)
    for j=1:size(XX,2)
        Z(i,j) = elev(XX(i,j),YY(i,j));
        for k=1:n_holes
            if (XX(i,j)-holes(k,1))^2+(YY(i,j)-holes(k,2))^2 <= radius(k)^2
                Z(i,j)=-5;
            end
        end  
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define dynamics and cost-to-go
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [0;0];
x_goal = [L;L];


% Dynamics
syms x1 x2 u1 u2
M = [x1+u1; x2+u2];
f(x1, x2, u1, u2) = [M];
fx = jacobian(f,[x1,x2]);
%fxx = jacobian(fx,[x1,x2]);
fu = jacobian(f,[u1,u2]);
%fuu = jacobian(fu,[u1,u2]);
%fux = jacobian(fu,[x1,x2]);

% Cost-to-go
% cost-to-go is just the elevation + 1 (to 
% include distance/total time)
syms l(x1,x2,u1, u2)
l(x1, x2, u1, u2) = elev(x1,x2) + 1;
lx = gradient(l,[x1,x2]); lxx = jacobian(lx,[x1,x2]);
lu = gradient(l,[u1,u2]); luu = jacobian(lu,[u1,u2]);
lux = jacobian(lu,[x1,x2]);

% final cost
Qf = 10*eye(n);
syms lf(xl, x2)
lf(xl, x2) = [xl;x2]'*Qf*[xl;x2];
lfx = gradient(lf,[x1,x2]); lfxx = jacobian(lfx,[xl,x2]);

% constraint functions
syms x1 x2 u1 u2
g = cell(n_holes,1);
gx = cell(n_holes,1);
gu = cell(n_holes,1);
for k=1:n_holes
    g{k} = radius(k)^2-norm([x1,x2]-holes(k,:),2)^2;
    gx{k} = gradient(g{k}, [x1, x2]);
    gu{k} = gradient(g{k}, [u1, u2]);
end

if plot_init==true
    surf(Z);
    xlabel('x');
    ylabel('y');
    zlabel('elevation');
end

save map.mat elev Z

save model.mat L elev Z x0 x_goal f fx fu l lx lxx lu luu ...
    lux lf lfx lfxx g gx gu %fxx fuu fux