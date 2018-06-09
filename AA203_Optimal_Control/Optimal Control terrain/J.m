function [res] = J(X,U)
global l lf N;

res = 0;
for i=1:size(X,1)-1
    res = res + sqrt(l(X(i,1), X(i,2),U(i,1),U(i,2)));
end
res = res + lf(X(end,1),X(end,2));
end