function Hnew = make_H_posdef(nseg,Hold)
  [V,D]=eigs(Hold,7*(nseg-1));
  for i = 1:7*(nseg-1)
      if D(i,i) < 0.001
          D(i,i)=0.1;
      end
  end
  Hnew = real(V*D*V');
return