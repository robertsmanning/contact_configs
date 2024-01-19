function [energy,grad,hess] = discrete_dna_penalty_en_grad_hess(zvec)
% Given a z-vector (q's then r's, but not boundary values), compute energy
%  and its gradient and Hessian

% Global vars (stiffnesses, bdy vals of r,q, contact params) set outside so
% known here
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

% Setup for building a sparse Hessian (ntriplets = number of times we'll put stuff in the Hessian)
ntriplets = 9*(nbp-2)*(nbp-2); % rr derivs
ntriplets = ntriplets + (nbp-2)*(3*12+3*12+3*16); % ii derivs for rq, qr, qq
ntriplets = ntriplets + (nbp-3)*(3*12+3*12+3*16); % i,i+1 derivs for rq, qr, qq
ntriplets = ntriplets + (nbp-3)*(3*12+3*12+3*16); % i,i-1 derivs for rq, qr, qq
ntriplets = ntriplets + (nbp-2)*(16);             % derivs of penalty function
I = zeros(ntriplets,1); J = zeros(ntriplets,1); X = zeros(ntriplets,1);

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
q(:,end) = qn'; normq(end)=norm(q(:,end));
q(:,1) = q0'; normq(1)=norm(q(:,1));

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

% Compute Hessian, starting with rr derivs
matrr = zeros(3*(nbp-2),3*(nbp-2)); % holds the rr derivs, later put into the Hessian

%   Compute contact portion of the rr derivs
cutoff = 3.0*rod_diam*nbp;
for i=1:nbp
    for j=(i+1):nbp
        if j-i > cutoff && i+(nbp-j) > cutoff
            xnm = r(:,i) - r(:,j); xnorm = norm(xnm,2);
            fpr = -Q*exp(-K*(xnorm-rod_diam))*(K*xnorm+1.0)/xnorm^2;
            fcombo = Q*exp(-K*(xnorm-rod_diam))*(K^2*xnorm^2+3*K*xnorm+3.0)/xnorm^2;
            contrib = fpr*eye(3)/xnorm+fcombo/xnorm^3*xnm*xnm';
            porti = 3*(i-2)+1:3*(i-1);
            portj = 3*(j-2)+1:3*(j-1);
	        if i > 1 && i < nbp  % means that ri is a variable
              matrr(porti,porti)=matrr(porti,porti)+contrib;
            end
	        if i > 1 && j > 1 && i < nbp && j < nbp  % means that ri and rj are variables
              matrr(porti,portj)=matrr(porti,portj)-contrib;
              matrr(portj,porti)=matrr(portj,porti)-contrib;
            end
 	        if j > 1 && j < nbp     % means that rj is a variable
              matrr(portj,portj)=matrr(portj,portj)+contrib;
	        end
        end
    end
end

%   Compute non-contact portion of the rr derivs
for i=1:nbp-2
    dirs = compute_ds(q(:,i)/normq(i)+q(:,i+1)/normq(i+1));
    dirsnext = compute_ds(q(:,i+1)/normq(i+1)+q(:,i+2)/normq(i+2));
    for j=1:3
        dj = dirs(:,j);
        djnext = dirsnext(:,j);
        range = 3*(i-1)+1:3*i;
        matrr(range,range) = matrr(range,range) + stiffs_for_mer(i,j)*dj*dj' + stiffs_for_mer(i+1,j)*djnext*djnext';
        if i < nbp-2
          matrr(range,range+3) = matrr(range,range+3) - stiffs_for_mer(i,j)*djnext*djnext';
          matrr(range+3,range) = matrr(range+3,range) - stiffs_for_mer(i,j)*djnext*djnext';
	    end
    end
end

% Assign the rr values to the (overall sparse) Hessian

ntriplets=0; % counter to number the Hessian entries as we enter them

for i=1:nbp-2
  for j=1:nbp-2
    for k1=1:3
      for k2=1:3
        ntriplets=ntriplets+1;
        if k1 <= k2
          X(ntriplets)=matrr(3*(i-1)+k1,3*(j-1)+k2);
        else
          X(ntriplets)=matrr(3*(i-1)+k2,3*(j-1)+k1);
        end
        I(ntriplets)=4*(nbp-2)+3*(i-1)+k1;
        J(ntriplets)=4*(nbp-2)+3*(j-1)+k2;
      end
    end
  end
end

% rq, qr, and qq derivs (NOTE: Hessian entries are sum over j=1,2,3; 
%   these are computed and recorded into I,J,X separately and then 
%   the "sparse" command will do the sum

for i=1:njunc-1
    dirs = compute_ds(q(:,i)/normq(i)+q(:,i+1)/normq(i+1));
    dirsnext = compute_ds(q(:,i+1)/normq(i+1)+q(:,i+2)/normq(i+2));
    dr = r(:,i+1)-r(:,i);
    drnext = r(:,i+2)-r(:,i+1);
    ddirs = compute_dds(q(:,i)/normq(i)+q(:,i+1)/normq(i+1));
    ddirsnext = compute_dds(q(:,i+1)/normq(i+1)+q(:,i+2)/normq(i+2));
    gmat = compute_drdotddir(q(:,i)/normq(i)+q(:,i+1)/normq(i+1),dr);
    gmatnext = compute_drdotddir(q(:,i+1)/normq(i+1)+q(:,i+2)/normq(i+2),drnext);
    d2dirs = compute_drdotd2ds(q(:,i)/normq(i)+q(:,i+1)/normq(i+1),dr);
    d2dirsnext = compute_drdotd2ds(q(:,i+1)/normq(i+1)+q(:,i+2)/normq(i+2),drnext);
    for j=1:3
        dj = dirs(:,j);
        djnext = dirsnext(:,j);
        ddjdq = zeros(3,4);
        ddjdqnext = zeros(3,4);
        d2djdq = zeros(4,4);
        d2djdqnext = zeros(4,4);
        gvec = gmat(j,:)';
        gvecnext = gmatnext(j,:)';
        bj = zeros(4,4);
        bjqbef = zeros(4,1);
        bjqmed = zeros(4,1);
        bjqaft = zeros(4,1);
        for j1 = 1:3
            for j2 = 1:4
                ddjdq(j1,j2) = ddirs(j,j1,j2);
                ddjdqnext(j1,j2) = ddirsnext(j,j1,j2);
            end
        end
        for j1 = 1:4
            for j2 = 1:4
                d2djdq(j1,j2) = d2dirs(j,j1,j2);
                d2djdqnext(j1,j2) = d2dirsnext(j,j1,j2);
            end
        end  
        for j1 = 1:4
            bjqbef(j1,1) = bq(j,j1,i);
            bjqmed(j1,1) = bq(j,j1,i+1);
            bjqaft(j1,1) = bq(j,j1,i+2);
        end
        for j1 = 1:4
            for j2 = 1:4
                bj(j1,j2)=b(j,j1,j2);
            end
        end
        qbef = q(:,i);
        qmed = q(:,i+1);
        qaft = q(:,i+2);
        mat = eye(4)/normq(i+1)-qmed*qmed'/(normq(i+1))^3;
        matnext = eye(4)/normq(i+2)-qaft*qaft'/(normq(i+2))^3;
        matprev = eye(4)/normq(i)-qbef*qbef'/(normq(i))^3;
        nmat = eye(4)/(normq(i+1))^3-3*qmed*qmed'/(normq(i+1))^5;
        range = 4*(njunc-1)+3*(i-1)+1:4*(njunc-1)+3*i;
        qrange = 4*(i-1)+1:4*i;
%                                 ii derivs (for rq/qr then qq)
        dum34 = stiffs_for_mer(i,j)*(tr(j,i)-avgs_for_mer(i,j))*ddjdq*mat;
        dum34=dum34-stiffs_for_mer(i+1,j)*(tr(j,i+1)-avgs_for_mer(i+1,j))*ddjdqnext*mat;
        dum34=dum34+stiffs_for_mer(i,j)*dj*dr'*ddjdq*mat-stiffs_for_mer(i+1,j)*djnext*drnext'*ddjdqnext*mat;
        for krow = 1:3
            for kcol = 1:4
                ntriplets=ntriplets+1;
                X(ntriplets) = dum34(krow,kcol);
                I(ntriplets) = range(krow);
                J(ntriplets) = qrange(kcol);
            end
        end        
        for krow = 1:4
            for kcol = 1:3
                ntriplets=ntriplets+1;
                X(ntriplets) = dum34(kcol,krow);
                I(ntriplets) = qrange(krow);
                J(ntriplets) = range(kcol);
            end
        end
        combo=mat*d2djdq*mat-(gvec*qmed'+qmed*gvec')/norm(qmed)^3-(qmed'*gvec)*nmat;
        dum44=stiffs_for_mer(i,j)*( (tr(j,i)-avgs_for_mer(i,j))*combo);
        dum44=dum44+stiffs_for_mer(i,j)*mat*ddjdq'*dr*dr'*ddjdq*mat;
        combonext=mat*d2djdqnext*mat-(gvecnext*qmed'+qmed*gvecnext')/norm(qmed)^3-(qmed'*gvecnext)*nmat;
        dum44=dum44+stiffs_for_mer(i+1,j)*( (tr(j,i+1)-avgs_for_mer(i+1,j))*combonext);
        dum44=dum44+stiffs_for_mer(i+1,j)*mat*ddjdqnext'*drnext*drnext'*ddjdqnext*mat;
        dum44=dum44+stiffs_for_mer(i,j+3)*(cay(j,i)-avgs_for_mer(i,j+3))/dfac(i)^2*(2*cay(j,i)*qbef*qbef'-2*bjqbef*qbef'-2*qbef*bjqbef');
        dum44=dum44+stiffs_for_mer(i+1,j+3)*(cay(j,i+1)-avgs_for_mer(i+1,j+3))/dfac(i+1)^2*(2*cay(j,i+1)*qaft*qaft'+2*bjqaft*qaft'+2*qaft*bjqaft');
        dum44=dum44+stiffs_for_mer(i,j+3)/dfac(i)^2*(2*bj-cay(j,i)*id)*qbef*qbef'*(-2*bj-cay(j,i)*id);
        dum44=dum44+stiffs_for_mer(i+1,j+3)/dfac(i+1)^2*(-2*bj-cay(j,i+1)*id)*qaft*qaft'*(2*bj-cay(j,i+1)*id);   
        for krow = 1:4
            for kcol = 1:4
                ntriplets=ntriplets+1;
                if krow <= kcol
                  X(ntriplets) = dum44(krow,kcol);
                else
                  X(ntriplets) = dum44(kcol,krow);
                end
                I(ntriplets) = qrange(krow);
                J(ntriplets) = qrange(kcol);
            end
        end 
%                   i,i+1 derivs (for rq/qr then qq)
        if i < njunc-1
            dum34=-stiffs_for_mer(i+1,j)*(tr(j,i+1)-avgs_for_mer(i+1,j))*ddjdqnext*matnext;
            dum34=dum34-stiffs_for_mer(i+1,j)*djnext*drnext'*ddjdqnext*matnext;
            for krow = 1:3
                for kcol = 1:4
                    ntriplets=ntriplets+1;
                    X(ntriplets) = dum34(krow,kcol);
                    I(ntriplets) = range(krow);
                    J(ntriplets) = qrange(kcol)+4;
                end
            end
            for krow = 1:4
                for kcol = 1:3
                    ntriplets=ntriplets+1;
                    X(ntriplets) = dum34(kcol,krow);
                    I(ntriplets) = qrange(krow)+4;
                    J(ntriplets) = range(kcol);
                end
            end
            dum44=stiffs_for_mer(i+1,j+3)*(cay(j,i+1)-avgs_for_mer(i+1,j+3))*( (2*cay(j,i+1)*qaft*qmed'-2*qaft*bjqmed'+2*bjqaft*qmed')/dfac(i+1)^2 + (-2*bj-cay(j,i+1)*id)/dfac(i+1) );
            dum44=dum44+stiffs_for_mer(i+1,j+3)/dfac(i+1)^2*(-2*bj-cay(j,i+1)*id)*qaft*qmed'*(-2*bj-cay(j,i+1)*id);
            dum44=dum44+stiffs_for_mer(i+1,j)*(tr(j,i+1)-avgs_for_mer(i+1,j))*mat*d2djdqnext*matnext;
            dum44=dum44+stiffs_for_mer(i+1,j)*mat*ddjdqnext'*drnext*drnext'*ddjdqnext*matnext;
            for krow = 1:4
                for kcol = 1:4
                    ntriplets=ntriplets+1;
                    X(ntriplets) = dum44(krow,kcol);
                    I(ntriplets) = qrange(krow);
                    J(ntriplets) = qrange(kcol)+4;
                end
            end
            for krow = 1:4
                for kcol = 1:4
                    ntriplets=ntriplets+1;
                    X(ntriplets) = dum44(kcol,krow);
                    I(ntriplets) = qrange(krow)+4;
                    J(ntriplets) = qrange(kcol);
                end
            end
        end
%               i,i-1 derivs (for rq/qr)
        if i>1
            dum34=stiffs_for_mer(i,j)*(tr(j,i)-avgs_for_mer(i,j))*ddjdq*matprev;
            dum34=dum34+stiffs_for_mer(i,j)*dj*dr'*ddjdq*matprev;
            for krow = 1:3
                for kcol = 1:4
                    ntriplets=ntriplets+1;
                    X(ntriplets) = dum34(krow,kcol);
                    I(ntriplets) = range(krow);
                    J(ntriplets) = qrange(kcol)-4;
                end
            end
            for krow = 1:4
                for kcol = 1:3
                    ntriplets=ntriplets+1;
                    X(ntriplets) = dum34(kcol,krow);
                    I(ntriplets) = qrange(krow)-4;
                    J(ntriplets) = range(kcol);
                end
            end            
        end
    end
%       ii derivs for penalty function
    dum44 = 4*penalty_weight*(q(:,i+1)'*q(:,i+1)-1)*id+8*penalty_weight*q(:,i+1)*q(:,i+1)';
    for krow = 1:4
        for kcol = 1:4
            ntriplets=ntriplets+1;
            if krow <= kcol
              X(ntriplets) = dum44(krow,kcol);
            else
              X(ntriplets) = dum44(kcol,krow);
            end
            I(ntriplets) = qrange(krow);
            J(ntriplets) = qrange(kcol);
        end
    end
end
hess = sparse(I,J,X,zlen,zlen);

