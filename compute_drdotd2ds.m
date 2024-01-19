function d2ds = compute_drdotd2ds(qvec,drvec)
% Given quaternion "qvec" and 3-vector "drvec"
%   compute second derivs of (drvec'*di) wrt q
% Stored in 3x4x4 array: d2ds(i,j,k) is (qj,qk) deriv of (drvec'*di)
d2ds = zeros(3,4,4);
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
d2d11 = (2/qnormsq)*(d2(1)*b3-d3(1)*b2+b3*qvec*dd2(1,:)-b2*qvec*dd3(1,:)-dd1(1,:)'*qvec');
d2d12 = (2/qnormsq)*(d2(2)*b3-d3(2)*b2+b3*qvec*dd2(2,:)-b2*qvec*dd3(2,:)-dd1(2,:)'*qvec');
d2d13 = (2/qnormsq)*(d2(3)*b3-d3(3)*b2+b3*qvec*dd2(3,:)-b2*qvec*dd3(3,:)-dd1(3,:)'*qvec');
d2d1 = drvec(1)*d2d11+drvec(2)*d2d12+drvec(3)*d2d13;
for j1=1:4
    for j2=1:4
        d2ds(1,j1,j2)=d2d1(j1,j2);
    end
end
d2d21 = (2/qnormsq)*(d3(1)*b1-d1(1)*b3+b1*qvec*dd3(1,:)-b3*qvec*dd1(1,:)-dd2(1,:)'*qvec');
d2d22 = (2/qnormsq)*(d3(2)*b1-d1(2)*b3+b1*qvec*dd3(2,:)-b3*qvec*dd1(2,:)-dd2(2,:)'*qvec');
d2d23 = (2/qnormsq)*(d3(3)*b1-d1(3)*b3+b1*qvec*dd3(3,:)-b3*qvec*dd1(3,:)-dd2(3,:)'*qvec');
d2d2 = drvec(1)*d2d21+drvec(2)*d2d22+drvec(3)*d2d23;
for j1=1:4
    for j2=1:4
        d2ds(2,j1,j2)=d2d2(j1,j2);
    end
end
d2d31 = (2/qnormsq)*(d1(1)*b2-d2(1)*b1+b2*qvec*dd1(1,:)-b1*qvec*dd2(1,:)-dd3(1,:)'*qvec');
d2d32 = (2/qnormsq)*(d1(2)*b2-d2(2)*b1+b2*qvec*dd1(2,:)-b1*qvec*dd2(2,:)-dd3(2,:)'*qvec');
d2d33 = (2/qnormsq)*(d1(3)*b2-d2(3)*b1+b2*qvec*dd1(3,:)-b1*qvec*dd2(3,:)-dd3(3,:)'*qvec');
d2d3 = drvec(1)*d2d31+drvec(2)*d2d32+drvec(3)*d2d33;
for j1=1:4
    for j2=1:4
        d2ds(3,j1,j2)=d2d3(j1,j2);
    end
end
return