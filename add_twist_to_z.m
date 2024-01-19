function znew = add_twist_to_z(nseg,zold,tw_angle)
% Adds twist angle "tw_angle" uniformly to a given z-vector "zold"
znew = zold; 
for i=1:nseg-1
    qold = zold(4*(i-1)+1:4*i);
    qnew = quatmult(qold,[0;0;sin(tw_angle*i/(2*nseg));cos(tw_angle*i/(2*nseg))]);
    znew(4*(i-1)+1:4*i) = qnew;
end
return