function gmat = compute_drdotddir(qvec,dr)
% Given quaternion "qvec" and 3-vector "dr", compute q-derivs of dr'*di
% Stored in 3x4 matrix: gmat(i,j) is d(dr'*di)/dqj
gmat = zeros(3,4);
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
for j1 = 1:4
    gmat(1,j1)=dd1(1,j1)*dr(1)+dd1(2,j1)*dr(2)+dd1(3,j1)*dr(3);
    gmat(2,j1)=dd2(1,j1)*dr(1)+dd2(2,j1)*dr(2)+dd2(3,j1)*dr(3);
    gmat(3,j1)=dd3(1,j1)*dr(1)+dd3(2,j1)*dr(2)+dd3(3,j1)*dr(3);
end