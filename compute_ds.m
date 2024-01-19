function dirs = compute_ds(qvec)
% Given quaternion "qvec", compute dj (returned as 3x3 matrix with dj as columns) 
dirs = zeros(3,3);
q1 = qvec(1); q2 = qvec(2); q3 = qvec(3); q4 = qvec(4);
qnormsq = q1^2+q2^2+q3^2+q4^2;
dirs(:,1) = (1/qnormsq)*[q1^2-q2^2-q3^2+q4^2;2*q1*q2+2*q3*q4;2*q1*q3-2*q2*q4];
dirs(:,2) = (1/qnormsq)*[2*q1*q2-2*q3*q4;-q1^2+q2^2-q3^2+q4^2;2*q1*q4+2*q2*q3];
dirs(:,3) = (1/qnormsq)*[2*q1*q3+2*q2*q4;-2*q1*q4+2*q2*q3;-q1^2-q2^2+q3^2+q4^2];
return