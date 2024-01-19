function [bend,twist,shear,stretch,contact,closest,where_closest,qnormcheck] = five_energies_plus(zvec)
% Given a z-vector (q's then r's, but not boundary values), compute five
%   components of energy, distance of closest approach, where it occurs,
%   and check on how close we are to having all |q| = 1

% Global vars (stiffnesses, bdy vals of r,q, contact params) set outside so known here
global avgs_for_mer
global stiffs_for_mer
global r0
global rn
global q0
global qn
global Q
global K
global rod_diam

% Length of z is 7(N-1); "nbp" is N+1, "njunc" is N 
s = size(zvec); zlen = s(1);
nbp = zlen/7+2; njunc = nbp-1;

penalty_weight=100; % Coefficient for energy term enforcing |q|=1

% Matrices q and r hold all q and r values, including bdy values
q = zeros(4,nbp);
r = zeros(3,nbp);
id = eye(4);

qnormcheck = 0.0;

% assign values to the vectors q and r for all i along the rod
for i=2:nbp-1
    q(:,i)=zvec(4*(i-2)+1:4*(i-1),1);
    qnorm = norm(q(:,i));
    if abs(qnorm-1) > qnormcheck
        qnormcheck = abs(qnorm-1);
    end
    r(:,i)=zvec(4*(nbp-2)+3*(i-2)+1:4*(nbp-2)+3*(i-1),1);
end

% impose the boundary conditions
r(:,end) = rn'; r(:,1) = r0';
q(:,end) = qn'; q(:,1) = q0';

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
    dirs = compute_ds(q(:,i+1)+q(:,i));
    tr(1,i) = dr'*dirs(:,1);
    tr(2,i) = dr'*dirs(:,2);
    tr(3,i) = dr'*dirs(:,3);
end
       
% Compute elastic energy components
bend = 0.0; twist = 0.0; shear = 0.0; stretch = 0.0;
for i=1:njunc
    bend = bend+stiffs_for_mer(i,4)*(cay(1,i)-avgs_for_mer(i,4))^2/2;
    bend = bend+stiffs_for_mer(i,5)*(cay(2,i)-avgs_for_mer(i,5))^2/2;
    twist = twist+stiffs_for_mer(i,6)*(cay(3,i)-avgs_for_mer(i,6))^2/2;
    shear = shear+stiffs_for_mer(i,1)*(tr(1,i)-avgs_for_mer(i,1))^2/2;
    shear = shear+stiffs_for_mer(i,2)*(tr(2,i)-avgs_for_mer(i,2))^2/2;
    stretch = stretch+stiffs_for_mer(i,3)*(tr(3,i)-avgs_for_mer(i,3))^2/2;
end

% Compute contact energy and dist of closest approach
contact = 0; xnorm = 0;
closest = 1; where_closest = [0 0];
cutoff = 3.0*rod_diam*nbp;
for i=1:nbp
    for j=(i+1):nbp
        if j-i > cutoff && i+(nbp-j) > cutoff
            xnm = r(:,i) - r(:,j); xnorm = norm(xnm,2);
            contact = contact + Q * exp(-K*(xnorm-rod_diam)) / xnorm;
            if xnorm < closest
               closest = xnorm;
               where_closest = [i j];
            end
        end
    end
end

