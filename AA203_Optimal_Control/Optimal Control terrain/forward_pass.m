function [state_new, control_new, J_temp] = forward_pass(state, control)

global N  n f g gu mu1 mu2 Quu_all Qux_all Qu_all;

% Parameters to tune
BIG = 1;
alpha = 0.5;
alpha1 = 1.05;
alpha2 = 1.05;
beta1 = 0.95;
beta2 = 0.95;

% Forward pass
feasible = false;
J_int = J(state, control);
e = BIG;
while feasible==false
    feasible = true;
    x = state(1,:);
    X_temp = state;
    U_temp = control;
    for k=1:N-1
        delta_x = x - state(k,:);
        X_temp(k,:) = x;
        
        cvx_begin;
        cvx_quiet true;
        variable delta_u_opt(n);
        x1 = x(1); x2 = x(2); u1 = control(k,1); u2 = control(k,2);
        minimize (0.5*delta_u_opt'*Quu_all{k}*delta_u_opt + delta_u_opt'*Qux_all{k}*delta_x' + delta_u_opt'*Qu_all{k});
        subject to
            norm(delta_u_opt) <= e;
            g{1}(x1,x2,u1,u2) + gu{1}(x1,x2,u1,u2)'*delta_u_opt <= 0;
            g{2}(x1,x2,u1,u2) + gu{2}(x1,x2,u1,u2)'*delta_u_opt <= 0;
        cvx_end;

        if cvx_status == "Infeasible"
            feasible = false;
            disp("Infeasible... Reducing bounds on delta_u");
            e = alpha*e;
            break;
        end
        
        U_temp(k,:) = control(k,:) + delta_u_opt';
        x = f(x(1),x(2),U_temp(k,1),U_temp(k,2))';
        
    end
    X_temp(N,:) = x;
    J_temp = J(X_temp, U_temp);
    J_temp = double(J_temp);
end

if J_temp < J_int
    disp("cost improved");
    state_new = X_temp;
    control_new = U_temp;
    mu1 = beta1*mu1;
    mu2 = beta2*mu2;
else
    disp("cost is worse");
    disp(double(J_temp));
    disp(double(J_int));
    mu1 = alpha1*mu1;
    mu2 = alpha2*mu2;
    state_new = state;
    control_new = control;
end

end