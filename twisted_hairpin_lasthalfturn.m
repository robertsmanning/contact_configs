% Follows the twisted compressed semicircle as we increase twist
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
nseg_final = 200; nseg = 50;

% Parameters for the energy
rod_diam = 0.03; contact_en_goal = 10;
cont_par_init = 100; cont_par_final = 10000;
K = cont_par_init;
contact_strength = contact_en_goal/(4*pi*rod_diam);
Q = contact_strength/nseg^2;

avgs_for_mer = zeros(nseg+1,6);
avgs_for_mer(:,3) = 1/nseg; % Intrinsic stretch (usu. rodlength/nseg) 

stiffs_for_mer = zeros(nseg+1,6);
k1 = 1;
k2 = k1;
k3 = 1.0;
a1 = 4800*16;
a2 = a1;
a3 = a1*10;
for i=1:nseg+1
    stiffs_for_mer(i,1:3) = [a1 a2 a3]*nseg;
    stiffs_for_mer(i,4:6) = [k1 k2 k3]*nseg;
end

timings = zeros(1,4); avg_steps = zeros(1,4);
tic
totsteps = 0; totruns = 0;

% Initial guess for energy-minimizing configuration: semicircle no twist
dt = 1/nseg;
R = (1)/(pi);
for i = 1:nseg+1
    rs(i,:) = [0 R*(cos(pi*(i-1)*dt)-1) R*sin(pi*(i-1)*dt)];
    s1 = sin(pi*(i-1)*dt/2); c1 = cos(pi*(i-1)*dt/2); 
    qs(i,:) = [s1 0 0 c1];
end
r0 = rs(1,:);  rn = rs(nseg+1,:);
q0 = [0 0 0 1]; qn = [1 0 0 0];

% Assemble all but first and last r's & q's into zvecs (vec of vars)
zvecs = zeros(7*(nseg-1),1);
rTemp = rs(2:end-1,:); qTemp = qs(2:end-1,:);
for i=1:nseg-1
    zvecs(4*(nseg-1) + 3*(i-1)+1:4*(nseg-1) + 3*i,1) = rTemp(i,:)';
    zvecs(4*(i-1)+1:4*i,1) = qTemp(i,:)';
end

% Minimize the discrete-rod energy
[zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,[],diagflag);
totsteps = totsteps + nsteps; totruns = totruns + 1;

[e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
[eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);

if diagflag > 0
  fprintf('Initial Solution has energy components\n');
  fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f\n',e1,e2,e3,e4,e5);
end

if fig_flag == 1
    plot_shape_to_screen(1,nseg,zeq,r0,rn,q0,qn);
end

% Save the solution to a file
filename = 'you_choose.txt';
fileID2 = fopen(filename,'wt');
z_with_bcs = [q0';zeq(1:4*(nseg-1));qn';r0';zeq(4*(nseg-1)+1:end);rn'];
fprintf(fileID2,'%f\r\n',z_with_bcs);
fclose(fileID2);

% Push the two ends toward each other
for k = 1:6
    rn = [0 -2/pi+(k/6)*(2/pi-0.1) 0];
    zvecs = zeq;
    [zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,[],diagflag);
    totsteps = totsteps + nsteps; totruns = totruns + 1;
    [e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
    [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);

    if diagflag > 0
      fprintf('%d of 6 steps of pushing ends together; energy components\n',k);
      fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f Closest %5.3f (indices %d %d)\n',e1,e2,e3,e4,e5,dis,whe(1),whe(2));
    end

    if fig_flag == 1
        plot_shape_to_screen(k,nseg,zeq,r0,rn,q0,qn);
    end

end

% Shifting from soft-contact to hard-contact
nk = 4;
for k = 1:nk
    par_ratio = cont_par_final/cont_par_init;
    K = cont_par_init*(par_ratio)^(1.0*k/4);

    zvecs = zeq;
    [zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,[],diagflag);
    totsteps = totsteps + nsteps; totruns = totruns + 1;
    [e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
    [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);

    if diagflag > 0
      fprintf('%d of 4 steps shifting to hard contact; energy components\n',k);
      fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f Closest %5.3f (indices %d %d)\n',e1,e2,e3,e4,e5,dis,whe(1),whe(2));
    end

end

elapsedtime = toc;
timings(1,1) = elapsedtime; avg_steps(1,1) = totsteps/totruns;
totsteps = 0; totruns = 0;

tic

% Increasing the twist to 5 pi
nk = 60;
qsep = qn;
for k = 1:nk
    wnew = k*pi/12;
    qtw = [0 0 sin(-wnew/2) cos(-wnew/2)];
    qn = quatmult(qsep,qtw);
    zvecs = add_twist_to_z(nseg,zeq,-pi/12);
    [zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,Happrox,diagflag);
    totsteps = totsteps + nsteps; totruns = totruns + 1;
    [e1,e2,e3,e4,e5,dis,whe]=five_energies_plus(zeq);
    [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);

    if diagflag > 0
      fprintf('Solution with twist %4.1f degrees has energy components\n',wnew*180/pi); 
      fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f Closest %5.3f (indices %d %d)\n',e1,e2,e3,e4,e5,dis,whe(1),whe(2));
    end

    if fig_flag == 1
        plot_shape_to_screen(k,nseg,zeq,r0,rn,q0,qn);
    end
end

elapsedtime = toc;
timings(1,2) = elapsedtime; avg_steps(1,2) = totsteps/totruns;
totsteps = 0; totruns = 0;

tic

% Interpolating
nk = 6;
nsegvals = round(10.^(linspace(log10(nseg),log10(nseg_final),nk+1)));
for k = 1:nk
    nsegold = nsegvals(k); nsegnew = nsegvals(k+1);
    fprintf("Interpolating from %d segs to %d\n",nsegold,nsegnew);
    zeq = interpolate_z(zeq,nsegold,nsegnew);
    nseg = nsegnew;
    avgs_for_mer = zeros(nseg+1,6);
    avgs_for_mer(:,3) = 1/nseg; 
    for i=1:nseg+1
      stiffs_for_mer(i,1:3) = [a1 a2 a3]*nseg;
      stiffs_for_mer(i,4:6) = [k1 k2 k3]*nseg;
    end
    Q = contact_strength/nseg^2;

    [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);
    fprintf("Energy of interpolant is %15.10f, Gradient is %15.10f\n",eneq,norm(gradeq,Inf))
    [e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
    fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f Closest %5.3f (indices %d %d)\n',e1,e2,e3,e4,e5,dis,whe(1),whe(2));

    zvecs = zeq;
    [zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,[],diagflag);
    totsteps = totsteps + nsteps; totruns = totruns + 1;
    [e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
    [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);

    if diagflag > 0
      fprintf('Interpolation step %d, energy components\n',k);
      fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f Closest %5.3f (indices %d %d)\n\n',e1,e2,e3,e4,e5,dis,whe(1),whe(2));
    end

end

elapsedtime = toc;
timings(1,3) = elapsedtime; avg_steps(1,3) = totsteps/totruns;
totsteps = 0; totruns = 0;

tic

% Increasing the twist to 6 pi
nk = 12;
results_to_save = zeros(nk,9);
% Happrox = [];
for k = 1:nk
    wnew = 5*pi+k*pi/12;
    qtw = [0 0 sin(-wnew/2) cos(-wnew/2)];
    qn = quatmult(qsep,qtw);
    zvecs = add_twist_to_z(nseg,zeq,-pi/12);
    [zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zvecs,Happrox,diagflag);
    totsteps = totsteps + nsteps; totruns = totruns + 1;
    [e1,e2,e3,e4,e5,dis,whe,qnormcheck]=five_energies_plus(zeq);
    [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);

    if diagflag > 0
      fprintf('Solution with twist %4.1f degrees has energy components\n',wnew*180/pi); 
      fprintf('Be %7.3f Tw %7.3f Shear %8.3f Stret %10.6f Contact %7.3f Closest %5.3f (indices %d %d)\n',e1,e2,e3,e4,e5,dis,whe(1),whe(2));
    end

    results_to_save(k,:) = [wnew*180/pi e1 e2 e3 e4 e5 dis qnormcheck norm(gradeq,"Inf")];

% Save the solution to a file
    filename = sprintf('iteration%i.txt',k);
    fileID2 = fopen(filename,'wt');
    z_with_bcs = [q0';zeq(1:4*(nseg-1));qn';r0';zeq(4*(nseg-1)+1:end);rn'];
    fprintf(fileID2,'%20.15f\r\n',z_with_bcs);
    fclose(fileID2);

    if fig_flag == 1
        plot_shape_to_screen(k,nseg,zeq,r0,rn,q0,qn);
    end
end

writematrix(results_to_save)

elapsedtime = toc;
timings(1,4) = elapsedtime; avg_steps(1,4) = totsteps/totruns;

format long
timings_are = timings
avg_steps_are = avg_steps