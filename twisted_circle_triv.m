warning('off','all')
global avgs_for_mer
global stiffs_for_mer
global Q K rod_diam
global r0
global rn
global q0
global qn
global nseg
eigtol = 1e-14;

diagflag = 1; % Flag for whether to output diagnostic info to screen
fig_flag = 0; % Flag for whether to create figures

% Number of segments in discrete rod
nseg = 200;

% Parameters for the energy
rod_diam = 0.03; contact_en_goal = 10;
contact_range_param = 10000;
% contact_range_param_final = 10000; kratio = contact_range_param_final/contact_range_param;
K = contact_range_param;
contact_strength = contact_en_goal/(4*pi*rod_diam);

Q = contact_strength/nseg^2;

avgs_for_mer = zeros(nseg+1,6);
% Column 3 of avgs is intrinsic stretch (usu. totalrodlength/nseg)
avgs_for_mer(:,3) = 1/nseg;   

stiffs_for_mer = zeros(nseg+1,6);
k1 = 1;
k2 = k1;
k3 = 1.0;
a1 = 4800*16;
a2 = a1;
a3 = 48000*16; 
lam = 0.1;  % shear force (which determines twist angle in the initial condition)
for i=1:nseg+1
    stiffs_for_mer(i,1:3) = [a1 a2 a3]*nseg;
    stiffs_for_mer(i,4:6) = [k1 k2 k3]*nseg;
end

timings = zeros(1,4); avg_steps = zeros(1,4);
tic
totsteps = 0; totruns = 0;

% Initial guess for energy-minimizing configuration
dt = 1/nseg;
rs = zeros(nseg+1,3);
qs = zeros(nseg+1,4);

sa = -2*(lam/a1)/(1+sqrt(1+4*(1-a1/a3)*(lam/a1)^2));
ca = sqrt(1-sa^2);
alp = asin(sa);
sa2 = sin(alp/2); 
ca2 = cos(alp/2);
R = (1+lam*sa/a3)/(2*pi*ca);
w = lam*R/(k3*ca)+(k1/k3-1)*2*pi*sa;
for i = 1:nseg+1
    rs(i,:) = [0 R*(cos(2*pi*(i-1)*dt)-1) R*sin(2*pi*(i-1)*dt)];
    s1 = sin(pi*(i-1)*dt); c1 = cos(pi*(i-1)*dt); 
    sw = sin(w*(i-1)*dt/2); cw = cos(w*(i-1)*dt/2);
    q1 = ca2*s1*cw+sa2*c1*sw;
    q2 =-ca2*s1*sw+sa2*c1*cw;
    q3 = ca2*c1*sw+sa2*s1*cw;
    q4 = ca2*c1*cw-sa2*s1*sw;
    qs(i,:) = [q1 q2 q3 q4];
end
r0 = rs(1,:);
rn = rs(nseg+1,:);
q0 = [0 0 0 1];
wapprox = w; % estimate of twist angle for initial condition
qn = [0 0 -sin(wapprox/2) -cos(wapprox/2)];

% Assemble all but the first and last r's and q's into zvecs (vector of
% variables)
zvecs = zeros(7*(nseg-1),1);
rTemp = rs(2:end-1,:);
qTemp = qs(2:end-1,:);
for i=1:nseg-1
    zvecs(4*(nseg-1) + 3*(i-1)+1:4*(nseg-1) + 3*i,1) = rTemp(i,:)';
    zvecs(4*(i-1)+1:4*i,1) = qTemp(i,:)';
end

% Minimize the discrete-rod energy
[zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,[],diagflag);
totsteps = totsteps + nsteps; totruns = totruns + 1;

fprintf('Solution with twist %4.1f degrees has energy components\n',wapprox*180/pi);
[e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f\n',e1,e2,e3,e4,e5);
[eneq,gradeq,hesseq] = discrete_dna_penalty_en_grad_hess(zeq);
Happrox = make_H_posdef(nseg,hesseq);

if fig_flag == 1
    plot_shape_to_screen(1,nseg,zeq,r0,rn,q0,qn);
end

% Save the solution to a file
filename = 'you_choose.txt';
fileID2 = fopen(filename,'wt');
z_with_bcs = [q0';zeq(1:4*(nseg-1));qn';r0';zeq(4*(nseg-1)+1:end);rn'];
fprintf(fileID2,'%f\r\n',z_with_bcs);
fclose(fileID2);

% Parameter continuation (increasing the twist)
%   (Want to continue until the twist energy jumps down, indicating we have
%     hopped from branch of circles to figure eights)
% Add twist relatively quickly up to one turn, then slowly to get 
%   close to jump (if k3/k1 > 1.375, do super-slowly to see buckled branch)
nk1 = 24; 
if k3/k1 <= 1.375
  nk2 = 144; dang = pi/36;
else
  nk2 = 360; dang = pi/360;
end
results_to_save = zeros(nk1+nk2+1,9);
results_to_save(1,:) = [180*wapprox/pi e1 e2 e3 e4 e5 dis qnormcheck norm(gradeq,"Inf")];

for k = 1:nk1
    wnew = k*pi/12;
    qn = [0 0 -sin(wnew/2) -cos(wnew/2)];
    zvecs = add_twist_to_z(nseg,zeq,pi/12);
    zvecs = zeq;
    [zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,[],diagflag);
    totsteps = totsteps + nsteps; totruns = totruns + 1;

    fprintf('Solution with twist %4.1f degrees has energy components\n',wnew*180/pi);
    [e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
    fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f Closest %5.3f (indices %d %d)\n',e1,e2,e3,e4,e5,dis,whe(1),whe(2));
    [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);

    results_to_save(k+1,:) = [wnew*180/pi e1 e2 e3 e4 e5 dis qnormcheck norm(gradeq,"Inf")];

        % Save the solution to a file
    filename = sprintf('iteration%i.txt',k);
    fileID2 = fopen(filename,'wt');
    z_with_bcs = [q0';zeq(1:4*(nseg-1));qn';r0';zeq(4*(nseg-1)+1:end);rn'];
    fprintf(fileID2,'%f\r\n',z_with_bcs);
    fclose(fileID2);

    if fig_flag == 1
      plot_shape_to_screen(1,nseg,zeq,r0,rn,q0,qn);
    end

end

k = 1;
while k <= nk2
    wnew = 2*pi+k*dang;
    qn = [0 0 -sin(wnew/2) -cos(wnew/2)];
    zvecs = add_twist_to_z(nseg,zeq,dang);
    zvecs = zeq;
    [zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,[],diagflag);
    totsteps = totsteps + nsteps; totruns = totruns + 1;

    fprintf('Solution with twist %4.1f degrees has energy components\n',wnew*180/pi);
    [e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
    fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f Closest %5.3f (indices %d %d)\n',e1,e2,e3,e4,e5,dis,whe(1),whe(2));
    [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);

    if e1 > 40
        fprintf('Jump from trivial branch detected\n')
        break
    end

    results_to_save(k+nk1+1,:) = [wnew*180/pi e1 e2 e3 e4 e5 dis qnormcheck norm(gradeq,"Inf")];

        % Save the solution to a file
    filename = sprintf('iteration%i.txt',k);
    fileID2 = fopen(filename,'wt');
    z_with_bcs = [q0';zeq(1:4*(nseg-1));qn';r0';zeq(4*(nseg-1)+1:end);rn'];
    fprintf(fileID2,'%f\r\n',z_with_bcs);
    fclose(fileID2);

    if fig_flag == 1
      plot_shape_to_screen(1,nseg,zeq,r0,rn,q0,qn);
    end

    k=k+1;         
end

elapsedtime = toc;
timings(1,1) = elapsedtime; avg_steps(1,1) = totsteps/totruns;

writematrix(results_to_save)

format long
timings_are = timings
avg_steps_are = avg_steps

[eneq,grade,hesseq] = discrete_dna_penalty_en_grad_hess(zeq);
w = eigs(hesseq,7*(nseg-1),'largestabs','Tolerance',eigtol);
w2 = mink(w,4);
fprintf('   Smallest eigs %15.10f %15.10f %15.10f %15.10f\n',w2(1),w2(2),w2(3),w2(4));