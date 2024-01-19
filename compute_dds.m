function ddirs = compute_dds(qvec)
% Given quaternion "qvec", compute derivatives of dj with respect to q
% Stored in 3x3x4 array: ddirs(i,j,k) is deriv of (di)_j wrt qk 
ddirs = zeros(3,3,4);
q1 = qvec(1); q2 = qvec(2); q3 = qvec(3); q4 = qvec(4);
qnormsq = q1^2+q2^2+q3^2+q4^2;
b1 = [0 0 0 1;0 0 1 0;0 -1 0 0;-1 0 0 0];
b2 = [0 0 -1 0;0 0 0 1;1 0 0 0;0 -1 0 0];
b3 = [0 1 0 0;-1 0 0 0;0 0 0 1;0 0 -1 0];
d1 = (1/qnormsq)*[q1^2-q2^2-q3^2+q4^2;2*q1*q2+2*q3*q4;2*q1*q3-2*q2*q4];
d2 = (1/qnormsq)*[2*q1*q2-2*q3*q4;-q1^2+q2^2-q3^2+q4^2;2*q1*q4+2*q2*q3];
d3 = (1/qnormsq)*[2*q1*q3+2*q2*q4;-2*q1*q4+2*q2*q3;-q1^2-q2^2+q3^2+q4^2];
dd1 = (2/qnormsq)*(d2*(b3*qvec)'-d3*(b2*qvec)');
dd2 = (2/qnormsq)*(d3*(b1*qvec)'-d1*(b3*qvec)');
dd3 = (2/qnormsq)*(d1*(b2*qvec)'-d2*(b1*qvec)');
for j1 = 1:3
    for j2 = 1:4
       ddirs(1,j1,j2)=dd1(j1,j2);
       ddirs(2,j1,j2)=dd2(j1,j2);
       ddirs(3,j1,j2)=dd3(j1,j2);
    end
end
return