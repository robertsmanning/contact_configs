function product = quatmult(qa,qb)
% Given quaternions "qa" and "qb", compute their product
  product = zeros(1,4);
  product(1) = qa(1)*qb(4)+qa(2)*qb(3)-qa(3)*qb(2)+qa(4)*qb(1);
  product(2) = -qa(1)*qb(3)+qa(2)*qb(4)+qa(3)*qb(1)+qa(4)*qb(2);
  product(3) = qa(1)*qb(2)-qa(2)*qb(1)+qa(3)*qb(4)+qa(4)*qb(3);
  product(4) = -qa(1)*qb(1)-qa(2)*qb(2)-qa(3)*qb(3)+qa(4)*qb(4);
end