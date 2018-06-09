function [] = backward_pass(state,control)
global n m g gx gu N lx lxx lu luu lux lfx lfxx fx fu ...
Quu_all Qux_all Qu_all mu1 mu2;

epsilon = 1e-3;
A = lfxx(state(end,1), state(end,2));
b = lfx(state(end,1), state(end,2));

for k = N-1:-1:1
    fxk = fx(state(k,1),state(k,2), control(k,1), control(k,2));
    fuk = fu(state(k,1),state(k,2), control(k,1), control(k,2));
    
    Qx = lx(state(k,1),state(k,2), control(k,1), control(k,2)) + fxk'*b;
    Qu = lu(state(k,1),state(k,2), control(k,1), control(k,2)) + fuk'*b;
    
    Qxx = lxx(state(k,1),state(k,2), control(k,1), control(k,2)) ...
        + fxk'*(A+mu1*eye(n))*fxk;% + b'*fxx(state(k,:),control(k,:));

    Quu = luu(state(k,1),state(k,2), control(k,1), control(k,2)) ...
        + fuk'*(A+mu1*eye(n))*fuk + mu2*eye(m); %+ ...
        %b'*fuu(state(k,:),control(k,:));

    Qux = lux(state(k,1),state(k,2), control(k,1), control(k,2)) ...
        + fuk'*(A+mu1*eye(n))*fxk; %+ ...
        %b'*fux(state(k,:),control(k,:));
    
    g_tilde = {}; % I think this is useless
    Ck = [];
    Dk = [];
    idx = 1;
    for i=1:length(g)
        gi = g{i};
        gui = gu{i};
        gxi = gx{i};
        x1 = state(k,1); x2 = state(k,2); u1 = control(k,1); u2 = control(k,2);
        if gi(x1,x2,u1,u2) >= -epsilon
%         if gi(state(k,1),state(k,2), control(k,1), control(k,2)) <= epsilon
            g_tilde{idx} = gi;  
            Ck = [Ck; gui(x1,x2,u1,u2)'];
            Dk = [Dk; -gxi(x1,x2,u1,u2)'];
            idx = idx+1;
        end
    end
    if length(Ck)==0
        Ck = zeros(1,n);
        Dk = zeros(1,n);
    end
    cvx_begin;
    cvx_quiet true
    variable delta_u(n);
    dual variable lambdas;
    minimize (0.5*delta_u'*Quu*delta_u+ delta_u'*Qu);
    subject to
        lambdas : Ck*delta_u == 0;
    cvx_end;

    Ck_tilde = [];
    Dk_tilde = [];
    for i=1:length(lambdas)
        if lambdas(i)>=0
            Ck_tilde = [Ck_tilde; Ck(i,:)];
            Dk_tilde = [Dk_tilde; Dk(i,:)];
        end
    end
        
    %W = ((Ck_tilde*(Quu\Ck_tilde'))\Ck_tilde)/(Quu);
    W = zeros(size(Ck_tilde,1),n);
    %H = Quu\(eye(n)-Ck_tilde'*W);
    H = Quu\(eye(n));
%     Kk = - H*Qux + W'*Dk_tilde;
    Kk = -H*Qux;
    jk = -H*Qu;
    A = Qxx + Kk'*Quu*Kk + Qux'*Kk + Kk'*Qux;
    b = Qx + Kk'*Quu*jk + Qux'*jk + Kk'*Qu;
    
    Quu_all{k} = Quu;
    Qux_all{k} = Qux;
    Qu_all{k} = Qu;
    
   
end

end