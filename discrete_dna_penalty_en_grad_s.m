function [energy,grad] = discrete_dna_penalty_en_grad_s(zvec)
% Given a z-vector (q's then r's, but not boundary values), compute energy
%  and its gradient

% Global vars (stiffnesses, bdy vals of r,q, contact params) set outside so known here
global avgs_for_mer
global stiffs_for_mer
global r0
global rn
global q0
global qn
global Q K rod_diam

% Length of z is 7(N-1); "nbp" is N+1, "njunc" is N 
s = size(zvec); zlen = s(1);
nbp = zlen/7+2; njunc = nbp-1;

penalty_weight=100; % Coefficient for energy term enforcing |q|=1

grad = zeros(zlen,1);

% Matrices q and r hold all q and r values, including bdy values
q = zeros(4,nbp); normq = zeros(1,nbp);
r = zeros(3,nbp);
id = eye(4);

% assign values to the vectors q and r for all i along the rod
for i=2:nbp-1
    q(:,i)=zvec(4*(i-2)+1:4*(i-1),1); normq(i)=norm(q(:,i));
    r(:,i)=zvec(4*(nbp-2)+3*(i-2)+1:4*(nbp-2)+3*(i-1),1);
end

% impose the boundary conditions
r(:,end) = rn'; r(:,1) = r0';
q(:,end) = qn'; normq(end) = norm(q(:,end));
q(:,1) = q0'; normq(1) = norm(q(:,1));

% define the B matrices used in theta and a definitions
b = zeros(3,4,4);
b(1,1,4) = 1; b(1,2,3) = 1; b(1,3,2)=-1; b(1,4,1) = -1;
b(2,1,3) = -1; b(2,2,4) = 1; b(2,3,1) = 1; b(2,4,2) = -1;
b(3,1,2) = 1; b(3,2,1) = -1; b(3,3,4) = 1; b(3,4,3) = -1;
bq = zeros(3,4,nbp);
for i=1:3
    bi = zeros(4,4);
    for j1 = 1:4
        for j2 = 1:4
            bi(j1,j2)=b(i,j1,j2);
        end
    end
    for k = 1:nbp
        dum = bi*q(:,k);
        for j1=1:4
            bq(i,j1,k)=dum(j1);
        end
    end
end

% cay = thetas, tr = a's, dfac = qa * qb (q dot its neighbor)
cay = zeros(3,njunc); tr = zeros(3,njunc); dfac = zeros(1,njunc);
for i=1:njunc
    dfac(i) = q(:,i+1)'*q(:,i);
    for k=1:3
        bk = zeros(4,4);
        for j1 = 1:4
            for j2 = 1:4
                bk(j1,j2)=b(k,j1,j2);
            end
        end
        cay(k,i) = 2/dfac(i)*q(:,i+1)'*(bk*q(:,i));        
    end
    dr = r(:,i+1)-r(:,i);
    dirs = compute_ds(q(:,i+1)/normq(i+1)+q(:,i)/normq(i));
    tr(1,i) = dr'*dirs(:,1);
    tr(2,i) = dr'*dirs(:,2);
    tr(3,i) = dr'*dirs(:,3);
end

% Elastic part of energy
energy = 0.0;
for i=1:njunc
    for j1 = 1:3
        energy = energy+stiffs_for_mer(i,j1+3)*(cay(j1,i)-avgs_for_mer(i,j1+3))^2/2;
        energy = energy+stiffs_for_mer(i,j1)*(tr(j1,i)-avgs_for_mer(i,j1))^2/2;
    end
    energy = energy+penalty_weight*(q(:,i+1)'*q(:,i+1)-1)^2;
end

% Self contact part of energy 
enq = 0; xnorm = 0;
for i=1:nbp
    for j=(i+1):nbp
        if j-i > 3.0*rod_diam*nbp && i+(nbp-j) > 3.0*rod_diam*nbp
            xnm = r(:,i) - r(:,j); xnorm = norm(xnm,2);
            enq = enq + exp(-K*(xnorm-rod_diam)) / xnorm;
        end
    end
end
energy = energy + Q*enq;

% Elastic portion of gradient
for i=1:njunc-1
    dirs = compute_ds(q(:,i)/normq(i)+q(:,i+1)/normq(i+1));
    dirsnext = compute_ds(q(:,i+1)/normq(i+1)+q(:,i+2)/normq(i+2));
    ddirs = compute_dds(q(:,i)/normq(i)+q(:,i+1)/normq(i+1));
    ddirsnext = compute_dds(q(:,i+1)/normq(i+1)+q(:,i+2)/normq(i+2));
    qthis = q(:,i+1);
    for j = 1:3
        bj = zeros(4,4);
        ddjdq = zeros(3,4);
        ddjdqnext = zeros(3,4);
        dj = dirs(:,j);
        djnext = dirsnext(:,j);
        for j1 = 1:4
            for j2 = 1:4
                bj(j1,j2)=b(j,j1,j2);
            end
        end
        for j1 = 1:3
            for j2 = 1:4
                ddjdq(j1,j2) = ddirs(j,j1,j2);
                ddjdqnext(j1,j2) = ddirsnext(j,j1,j2);
            end
        end
        grad(4*(i-1)+1:4*i,1) = grad(4*(i-1)+1:4*i,1) + stiffs_for_mer(i,j+3)*(cay(j,i)-avgs_for_mer(i,j+3))*(-cay(j,i)/dfac(i)*q(:,i)+2/dfac(i)*bj*q(:,i));
        grad(4*(i-1)+1:4*i,1) = grad(4*(i-1)+1:4*i,1) + stiffs_for_mer(i+1,j+3)*(cay(j,i+1)-avgs_for_mer(i+1,j+3))*(-cay(j,i+1)/dfac(i+1)*q(:,i+2)-2/dfac(i+1)*bj*q(:,i+2));
        mat = eye(4)/normq(i+1)-qthis*qthis'/(normq(i+1))^3;
        grad(4*(i-1)+1:4*i,1) = grad(4*(i-1)+1:4*i,1) + stiffs_for_mer(i,j)*(tr(j,i)-avgs_for_mer(i,j))*mat*ddjdq'*(r(:,i+1)-r(:,i));
        grad(4*(i-1)+1:4*i,1) = grad(4*(i-1)+1:4*i,1) + stiffs_for_mer(i+1,j)*(tr(j,i+1)-avgs_for_mer(i+1,j))*mat*ddjdqnext'*(r(:,i+2)-r(:,i+1));
        grad(4*(njunc-1)+3*(i-1)+1:4*(njunc-1)+3*i,1) = grad(4*(njunc-1)+3*(i-1)+1:4*(njunc-1)+3*i,1) + stiffs_for_mer(i,j)*(tr(j,i)-avgs_for_mer(i,j))*dj;
        grad(4*(njunc-1)+3*(i-1)+1:4*(njunc-1)+3*i,1) = grad(4*(njunc-1)+3*(i-1)+1:4*(njunc-1)+3*i,1) + stiffs_for_mer(i+1,j)*(tr(j,i+1)-avgs_for_mer(i+1,j))*(-djnext);
    end
    grad(4*(i-1)+1:4*i,1) = grad(4*(i-1)+1:4*i,1) + 2*penalty_weight*(q(:,i+1)'*q(:,i+1)-1)*2*q(:,i+1);
end

% Contact portion of gradient
xnorm = 0; cutoff = 3.0*rod_diam*nbp;
elecgrad = zeros(zlen,1); ofst = 4*(njunc-1);
for i=1:nbp
    for j=(i+1):nbp
        if j-i > cutoff && i+(nbp-j) > cutoff
          xdiff = r(:,i)-r(:,j); xnorm = norm(xdiff);
          fac = -exp(-K*(xnorm-rod_diam))*(1+K*xnorm)*xdiff/(xnorm^3);
          if i > 1 && i < nbp
             elecgrad(ofst+3*(i-2)+1:ofst+3*(i-1))=elecgrad(ofst+3*(i-2)+1:ofst+3*(i-1))+fac;
          end
          if j > 1 && j < nbp
             elecgrad(ofst+3*(j-2)+1:ofst+3*(j-1))=elecgrad(ofst+3*(j-2)+1:ofst+3*(j-1))-fac;
          end
        end
    end
end
elecgrad = Q*elecgrad;

grad = grad + elecgrad;

