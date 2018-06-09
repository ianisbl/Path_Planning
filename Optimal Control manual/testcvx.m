cvx_begin;
variable delta_u_opt(n);
x1 = 5; x2 = 6; u1 = 0; u2 = 0;
minimize (0.5*delta_u_opt'*Quu_all{1}*delta_u_opt + delta_u_opt'*Qux_all{1}*delta_x' + delta_u_opt'*Qu_all{1});
subject to
    norm(delta_u_opt) <= e;  
    g{1}(x1,x2,u1,u2) + gu{1}(x1,x2,u1,u2)'*delta_u_opt <= 0;
    g{2}(x1,x2,u1,u2) + gu{2}(x1,x2,u1,u2)'*delta_u_opt <= 0;
cvx_end;