function znew = interpolate_z(zold,nsegold,nsegnew)
% Given a z-vector (q then r, without bdy values) for "nsegold" number
%    of segments, use interpolation to generate a new z-vector "znew"
%    with "nsegnew" number of segments

  % Bdy values passed in as global variables
  global r0
  global rn
  global q0
  global qn

  znew = zeros(7*(nsegnew-1),1);

  % Build padded vectors that incorporate bdy values
  roldpadded = [r0' ; zold(4*(nsegold-1)+1:end) ; rn'];
  qoldpadded = [q0' ; zold(1:4*(nsegold-1)) ; qn'];

  dsold = 1/nsegold; % Arclength width of old-discretization
  iold = 1; % Counter for where we are in old-discretization
  sold = dsold; % Arclength value in old-discretization
  % Loop over counter for new-discretization
  for inew = 1:(nsegnew-1)
      snew = inew/nsegnew; % Arclength value for new-discretization
      while snew > sold % Advance old-discretization arclength until it
                        %   equals or exceeds new-discretization arclength
                        %   (also update old-discretization counter)
          iold=iold+1; sold=sold+dsold;
      end
      part = (snew-(sold-dsold))/dsold; % How far snew is within [sold-dsold,sold]
      qprev = qoldpadded(4*(iold-1)+1:4*iold,1); % Prev q in old
      qnext = qoldpadded(4*iold+1:4*iold+4,1); % Next q in old
      znew(4*(inew-1)+1:4*inew,1) = qprev+part*(qnext-qprev); % Interpolate q (then make it norm 1)
      fac = norm(znew(4*(inew-1)+1:4*inew,1));
      znew(4*(inew-1)+1:4*inew,1) = znew(4*(inew-1)+1:4*inew,1)/fac; 
      rprev = roldpadded(3*(iold-1)+1:3*iold,1); % Prev r in old
      rnext = roldpadded(3*iold+1:3*iold+3,1); % Next r in old
      znew(4*(nsegnew-1)+3*(inew-1)+1:4*(nsegnew-1)+3*inew,1)=rprev+part*(rnext-rprev); % Interpolate r
  end
end