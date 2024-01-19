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
mintol = 1e-5;
opttol = 1e-4;

diagflag = 1; % Flag for whether to output diagnostic info to screen
              % 1 shows minimal, 2 shows more
fig_flag = 0; % Flag for whether to create figures

% Number of segments in discrete rod
nseg = 50; nseg_final = 200;

% Parameters for the energy
rod_diam = 0.03; contact_en_goal = 10;
contact_range_param = 10000;
K = contact_range_param;
contact_strength = contact_en_goal/(4*pi*rod_diam);
Q = contact_strength/nseg^2;

avgs_for_mer = zeros(nseg+1,6);
% Column 3 of avgs is intrinsic stretch (usu. totalrodlength/nseg)
avgs_for_mer(:,3) = 1/nseg;  
intrins_bend = pi/3; % Gets turned on later

stiffs_for_mer = zeros(nseg+1,6);
k1 = 1;
k2 = k1;
k3 = 1.0;
a1 = 4800*16;
a2 = a1;
a3 = 48000*16; 
for i=1:nseg+1
    stiffs_for_mer(i,1:3) = [a1 a2 a3]*nseg;
    stiffs_for_mer(i,4:6) = [k1 k2 k3]*nseg;
end

timings = zeros(1,4); avg_steps = zeros(1,4);
tic
totsteps = 0; totruns = 0;

% Initial guess: slightly buckled rod
dt = 1/nseg;
rs = zeros(nseg+1,3);
qs = zeros(nseg+1,4);

eps = 0.1;
for i = 1:nseg+1
    s = (i-1)*dt;
    rs(i,:) = [0 eps*(1-cos(2*pi*s))/(2*pi) s+eps^2*(sin(4*pi*s)-4*pi*s)/(16*pi)];
    qs(i,:) = [0 0 sin(eps*sin(2*pi*s)) cos(eps*sin(2*pi*s))];
end
r0 = rs(1,:);
rn = rs(nseg+1,:);

q0 = [0 0 0 1];
qn = [0 0 0 1];

% Assemble all but the first and last r's and q's into zvecs (vector of
% variables)
zvecs = zeros(7*(nseg-1),1);
rTemp = rs(2:end-1,:);
qTemp = qs(2:end-1,:);
for i=1:nseg-1
    zvecs(4*(nseg-1) + 3*(i-1)+1:4*(nseg-1) + 3*i,1) = rTemp(i,:)';
    zvecs(4*(i-1)+1:4*i,1) = qTemp(i,:)';
end

% Minimize the discrete-rod energy for original boundary conditions
[zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,[],diagflag);
totsteps = totsteps + nsteps; totruns = totruns + 1;

fprintf('Initial Solution has energy components\n');
[e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f\n',e1,e2,e3,e4,e5);
[eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);

if fig_flag == 1
    plot_shape_to_screen(1,nseg,zeq,r0,rn,q0,qn);
end

% Save the solution to a file
filename = 'you_choose.txt';
fileID2 = fopen(filename,'wt');
z_with_bcs = [q0';zeq(1:4*(nseg-1));qn';r0';zeq(4*(nseg-1)+1:end);rn'];
fprintf(fileID2,'%f\r\n',z_with_bcs);
fclose(fileID2);

% Parameter continuation (lowering the top of the strut)
nk = 6;
for k = 1:nk
    zold = rn(3);
    rn = [0 0 (1-eps^2/4)*(nk-k)/nk+(0.25)*(k/nk)];
    znew = rn(3);
    qn = [0 0 0 1];
    zvecs = zeq;
    zvecs(4*(nseg-1)+3:3:7*(nseg-1))=zvecs(4*(nseg-1)+3:3:7*(nseg-1))*znew/zold;
    zvecs(4*(nseg-1)+1:3:7*(nseg-1))=zvecs(4*(nseg-1)+1:3:7*(nseg-1))*zold/znew;
    zvecs(4*(nseg-1)+2:3:7*(nseg-1))=zvecs(4*(nseg-1)+2:3:7*(nseg-1))*zold/znew;
    [e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zvecs);
    [zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,Happrox,diagflag);
    totsteps = totsteps + nsteps; totruns = totruns + 1;
    fprintf('%d of %d steps of pushing down the top; energy components\n',k,nk);
    [e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
    fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f Closest %5.3f (indices %d %d)\n',e1,e2,e3,e4,e5,dis,whe(1),whe(2));
    [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);
    Happrox = make_H_posdef(nseg,Happrox);

    if fig_flag == 1
        plot_shape_to_screen(k+1,nseg,zeq,r0,rn,q0,qn);
    end
end

% Use interpolation to get from orig nseg to final, and turn on intrinsic bend
nkint = 6;
nsegvals = round(10.^(linspace(log10(nseg),log10(nseg_final),nkint+1)));
for k = 1:nkint

    nsegold = nsegvals(k); nsegnew = nsegvals(k+1);
    zeq = interpolate_z(zeq,nsegold,nsegnew);
    nseg = nsegnew;

    fprintf("Interpolating from %d segs to %d\n",nsegold,nsegnew);

    Q = contact_strength/nseg^2;
    avgs_for_mer = zeros(nseg+1,6);
    avgs_for_mer(:,3) = 1/nseg;
    avgs_for_mer(:,4) = (k/nkint)*intrins_bend/nseg;
    for i=1:nseg+1
      stiffs_for_mer(i,1:3) = [a1 a2 a3]*nseg;
      stiffs_for_mer(i,4:6) = [k1 k2 k3]*nseg;
    end

    [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);
    fprintf("Energy of interpolant is %f",eneq)
    fprintf(" Gradient is %f\n",norm(gradeq,Inf))
    [e1,e2,e3,e4,e5,dis,whe]=five_energies_plus(zeq);
    fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f Closest %5.3f (indices %d %d)\n',e1,e2,e3,e4,e5,dis,whe(1),whe(2));

    zvecs = zeq;
    [zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,[],diagflag);
    totsteps = totsteps + nsteps; totruns = totruns + 1;
    fprintf('Solution has energy components\n');
    [e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
    fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f Closest %5.3f (indices %d %d)\n',e1,e2,e3,e4,e5,dis,whe(1),whe(2));
    [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);
end

elapsedtime = toc;
timings(1,1) = elapsedtime; avg_steps(1,1) = totsteps/totruns;
totsteps = 0; totruns = 0;

tic

% Increase the twist to 1.5 turns in steps of pi/6
nk = 18; nk2 = 72;
results_to_save = zeros(nk+nk2+1,9);
results_to_save(1,:) = [0.0 e1 e2 e3 e4 e5 dis qnormcheck norm(gradeq,Inf)];
qsep = qn;
for k = 1:nk
    wnew = k*pi/6;
    qtw = [0 0 sin(-wnew/2) cos(-wnew/2)];
    qn = quatmult(qsep,qtw);
    zvecs = add_twist_to_z(nseg,zeq,-pi/6);

    [zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,Happrox,diagflag);
    totsteps = totsteps + nsteps; totruns = totruns + 1;

    fprintf('Solution with twist %4.1f degrees has energy components\n',wnew*180/pi);
    [e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
    results_to_save(k+1,:) = [wnew*180/pi e1 e2 e3 e4 e5 dis qnormcheck norm(gradeq,Inf)];
    fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f Closest %5.3f (indices %d %d)\n',e1,e2,e3,e4,e5,dis,whe(1),whe(2));
    [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);
    Happrox = make_H_posdef(nseg,Happrox);

    % Save the solution to a file
    filename = sprintf('iteration%i.txt',k);
    fileID2 = fopen(filename,'wt');
    z_with_bcs = [q0';zeq(1:4*(nseg-1));qn';r0';zeq(4*(nseg-1)+1:end);rn'];
    fprintf(fileID2,'%f\r\n',z_with_bcs);
    fclose(fileID2);
end

% Increase the twist to three turns in steps of pi/24
Happrox = [];
for k = 1:nk2
    wnew = 3*pi+k*pi/24;
    qtw = [0 0 sin(-wnew/2) cos(-wnew/2)];
    qn = quatmult(qsep,qtw);

    zvecs = add_twist_to_z(nseg,zeq,-pi/24);
    [zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,Happrox,diagflag);
    totsteps = totsteps + nsteps; totruns = totruns + 1;

    [eneq,gradeq]=discrete_dna_penalty_en_grad_s(zeq);
    [e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
    fprintf('Solution with twist %4.1f degrees has energy components\n',wnew*180/pi);
    fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f Closest %5.3f (indices %d %d)\n',e1,e2,e3,e4,e5,dis,whe(1),whe(2));
    results_to_save(k+nk+1,:) = [wnew*180/pi e1 e2 e3 e4 e5 dis qnormcheck norm(gradeq,Inf)];

% Save the solution to a file
    filename = sprintf('iteration%i.txt',k+nk);
    fileID2 = fopen(filename,'wt');
    z_with_bcs = [q0';zeq(1:4*(nseg-1));qn';r0';zeq(4*(nseg-1)+1:end);rn'];
    fprintf(fileID2,'%f\r\n',z_with_bcs);
    fclose(fileID2);

    if fig_flag == 1
        plot_shape_to_screen(k,nseg,zeq,r0,rn,q0,qn);
    end
end
writematrix(results_to_save)
elapsedtime = toc;
timings(1,2) = elapsedtime; avg_steps(1,2) = totsteps/totruns;

format long
timings_are = timings
avg_steps_are = avg_steps

[eneq,gradeq,hesseq]=discrete_dna_penalty_en_grad_hess(zeq);
w = eigs(hesseq,7*(nseg-1),'largestabs','Tolerance',eigtol);
w2 = mink(w,4);
fprintf('   Smallest eigs %15.10f %15.10f %15.10f %15.10f\n',w2(1),w2(2),w2(3),w2(4));
