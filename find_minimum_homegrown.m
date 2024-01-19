function [zeq,nsteps,Happrox] = find_minimum_homegrown(nseg,zguess,initH,diagflag)
% Given initial guess "zguess" with "nseg" rod segments
%   and initial guess "initH" of inverse-Hessian (can be empty)
% find local minimum "zeq".
%
% "diagflag" is additional input argument controlling what gets printed to
%   screen (0 means nothing, 1 means some, 2 means max)
%
% Function also outputs number of steps and final approx of inv-Hessian

gradtol = 10^(-6); % Tolerance (stop when norm(grad) is below this)
zeq = zguess;

% Compute E,gradient for initial guess (needed to start while loop below)
[eneq,gradeq] = discrete_dna_penalty_en_grad_s(zguess);
Happrox = eye(7*(nseg-1));
if diagflag > 0
  fprintf("\n");
end

% If inverse-Hessian approx that was input was empty, call function
%   to use true inverse-Hessian to generate reasonable positive-definite
%   initial guess
if size(initH) == 0
  initH = set_H_to_inv_of_adjusted_true_hessian(nseg,zguess);
end

% Call minimization algorithm (see below) to find local minimum "zeq"
%   Also return:
%     final approx of inverse-Hessian, 
%     counters for number of steps of three types,
%     flag "fail" saying if algorithm succeeded
[zeq,Happrox,count1,count2,count3,fail] = min_alg(nseg,zguess,initH,gradtol,diagflag);
nsteps = count1+count2+count3;
if diagflag > 1
  fprintf("  Done with min_alg after %d steps, fail flag is %d\n",nsteps,fail);
end
if fail == -1 % This flag value indicates we did not reach gradtol but it was
              %   determined that we cannot improve our estimate of local min
  [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);
  if diagflag > 0
    fprintf("  %d steps (%d BFGS, %d bisec, %d stpst desc), normgrad %12.10f (cannot be improved)\n",count1+count2+count3,count1,count2,count3,norm(gradeq));
  end
  return
end
if fail == 1 % This flag value indicates we did not succeed
  [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);
  fprintf("  Alg failed, %d steps (%d BFGS, %d bisec, %d stpst desc), normgrad %12.10f\n",count1+count2+count3,count1,count2,count3,norm(gradeq));
  return
end

% If we reach this point, fail flag = 0, so we hit grad tolerance 
[eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);
if diagflag > 0
  fprintf("  %d steps (%d BFGS, %d bisec, %d stpst desc), normgrad %12.10f\n",nsteps,count1,count2,count3,norm(gradeq));
end
return

function [zeq,H,count1,count2,count3,fail] = min_alg(nseg,zguess,initH,gradtol,diagflag)
% Given initial guess "zguess" with "nseg" rod segments
%   initial guess "Hinit" of inverse-Hessian, and
%   desired tolerance "gradtol" for gradient,
% find local minimum "zeq".
%
% "diagflag" is additional input argument controlling what gets printed to
%   screen (0 means nothing, 1 means some, 2 means max)
%
% Function also outputs final approx of inv-Hessian, number of steps
%   of three different types that were taken, and a flag indicating status

% Global variables used for line-search function
  global current_guess
  global direc

  warning('off','all')
  fail=0;
  zeq=zguess; H = initH; 

% Timers for how long algorithm spends doing various types of work
  fminbnd_time = 0.0; bfgsupdates_time = 0.0; grad_time = 0.0; 
  hessupdate_time = 0.0; wolfebisec_time = 0.0; steep_time = 0.0;

  allclock = tic; tmpclock = tic; 
  [eneq,gradeq] = discrete_dna_penalty_en_grad_s(zeq);
  grad_time = grad_time + toc(tmpclock);
  if diagflag > 1
    fprintf("  Initial Energy %15.10f Gradient %20.10f\n",eneq,norm(gradeq))
  end

  fminbnd_tol = 10^(-2); count1 = 0; count2 = 0; count3 = 0;

  % Continue steps below until norm(grad) below its tolerance, or
  %   we determine no further progress is possible, or we declare failure
  while norm(gradeq) > gradtol
    % Compute search direction; update "current_guess"
    tmpclock = tic;
    direc = real(-H*gradeq); 
    current_guess = zeq;
    bfgsupdates_time = bfgsupdates_time + toc(tmpclock);
    
    % Call fminbnd to try to find alpha so that 
    %   zeqnew = current_guess + alpha*direc is an acceptable new guess
    tmpclock = tic;
    opts=optimset('TolX',fminbnd_tol,'Display','off');
    [alpha,fval,exitflag,output] = fminbnd(@line_search_g,0,2,opts); 
    fminbnd_time = fminbnd_time + toc(tmpclock);

    % Compute E and grad for zeqnew
    tmpclock = tic;
    zeqnew=zeq+alpha*direc;
    [ennew,gradnew] = discrete_dna_penalty_en_grad_s(zeqnew);
    grad_time = grad_time + toc(tmpclock);
    if diagflag > 1
        fprintf("    fminbnd alpha = %f, %d iterations, exitflag = %d, En %15.10f, Gr %20.10f\n",alpha,output.iterations,exitflag,ennew,norm(gradnew));
    end

    % Assess Wolfe conditions by computing wolfe1 and wolfe2
    tmpclock = tic;
    dotprod = direc'*gradeq; dotprodnew = direc'*gradnew;
    wolfe1=ennew-eneq+0.0001*alpha*dotprod;
    wolfe2 = -dotprodnew+0.9*dotprod;
    bfgsupdates_time = bfgsupdates_time + toc(tmpclock);

    status = 0; steepest_status = 0;
    if wolfe1>0 || wolfe2>0 % If Wolfe test fails, 
      if norm(gradeq) < 0.1 % For small gradient, try a steepest descent step
        tmpclock = tic;
        [alp,zeqnew,steepest_status] = steepest_descent_step(zeq,eneq,-gradeq,diagflag);
        steep_time = steep_time + toc(tmpclock);
        if steepest_status == 0 % If steepest fails, we are as close as poss
                                %   so exit function
          fail = -1;
          all_time = toc(allclock);
          if diagflag > 0
            fprintf("Total time %15.10f\n",all_time);
            fprintf("  Gradient    %9.3f  Hessian %9.3f\n",grad_time,hessupdate_time);
            fprintf("  BFGSupdates %9.3f  fminbnd %9.3f  Wolfebisec %9.3f  Steepest %9.3f\n", ...
              bfgsupdates_time,fminbnd_time,wolfebisec_time,steep_time);
          end
          return
        else % If steepest succeeds, update ennew and gradnew and continue on
          [ennew,gradnew] = discrete_dna_penalty_en_grad_s(zeqnew);
        end
      else % For larger gradient, try Wolfe-bisection first, then steepest as backup
        if diagflag > 1
          fprintf("    Switching to bisection, fminbnd alpha was %10.8f, dotprod is %20.15f\n",alpha,dotprod);
          fprintf("     wolfe1 = %20.15f, wolfe2 = %20.15f\n",wolfe1,wolfe2);
          fprintf("     ennew %20.15f eneq %20.15f\n",ennew,eneq);
        end
        tmpclock = tic;
        [alpha,status] = wolfe_bisection(zeq,direc,diagflag);
        wolfebisec_time = wolfebisec_time + toc(tmpclock);
        if status > 0 % If bisection worked, update zeqnew, ennew, gradnew and move on
          tmpclock = tic;
          zeqnew = zeq + alpha*direc;
          [ennew,gradnew] = discrete_dna_penalty_en_grad_s(zeqnew);
          bfgsupdates_time = bfgsupdates_time + toc(tmpclock);
        else % If bisection didn't work, try steepest descent step
          tmpclock = tic;
          [alp,zeqnew,steepest_status] = steepest_descent_step(zeq,eneq,-gradeq,diagflag);
          steep_time = steep_time + toc(tmpclock);
          if steepest_status == 0 % If steepest fails, we are as close as poss
                                  %   so exit function
            fail = -1;
            all_time = toc(allclock);
            if diagflag > 0
              fprintf("Total time %15.10f\n",all_time);
              fprintf("  Gradient    %9.3f  Hessian %9.3f\n",grad_time,hessupdate_time);
              fprintf("  BFGSupdates %9.3f  fminbnd %9.3f  Wolfebisec %9.3f  Steepest %9.3f\n", ...
                bfgsupdates_time,fminbnd_time,wolfebisec_time,steep_time);
            end
            return
          else % If steepest succeeds, update ennew and gradnew and move on
            tmpclock = tic;
            [ennew,gradnew] = discrete_dna_penalty_en_grad_s(zeqnew);
            grad_time = grad_time + toc(tmpclock);
          end
        end
      end
    end

    % If we reach here, we were able to take a step, so assess what kind of
    %   step it was, and continue to next step
    if status > 0 
      if diagflag > 1
        fprintf("  Bisec step, alpha is %10.8f, En %15.10f, Gr %20.10f\n",alpha,ennew,norm(gradnew))
      end
      count2 = count2+1;
    else
      if steepest_status > 0 % Tried steepest descent and it worked
        if diagflag > 1
          fprintf("  Steepest descent step, alpha is %15.14f, En %15.10f, Gr %20.10f\n",alp, ennew,norm(gradnew))
        end
        count3 = count3+1;
      else % otherwise, BFGS step worked
        if diagflag > 1
          fprintf("  BFGS step, alpha is %10.8f, En %15.10f, Gr %20.10f\n",alpha,ennew,norm(gradnew));
        end
        count1 = count1 + 1;
      end
    end
    if count1+count2+count3 >= 2000 % If we reach 2000 steps of any type, give up
      fail=1;
      return
    end

    % Having taken a step, update approximate inv-Hessian via BFGS scheme
    tmpclock = tic;
    y = gradnew-gradeq; w = zeqnew-zeq; term = y'*w;
    if term > 0
      vq = H*y; outprod = (vq*w')/term;
      H = H + (term+y'*vq)/term^2*(w*w') - outprod - outprod';
    end
    bfgsupdates_time = bfgsupdates_time + toc(tmpclock);

    zeq=zeqnew; gradeq=gradnew; eneq = ennew;

    % Every 20 steps, update approx inv-Hess to true inv-Hess if Hess is pos-definite
    if mod(count1+count2+count3,20) == 0
        if diagflag > 0
          fprintf("  After %d steps, grad is %12.8f   ",count1+count2+count3,norm(gradeq))
        end
        tmpclock = tic;
        H = set_H_to_true_hessian_inv_if_posdef(nseg,zeq,H,diagflag);
        hessupdate_time = hessupdate_time + toc(tmpclock);
        if diagflag > 0
          fprintf("\n");
        end
    end 
    % If did a Bisection or steepest descent step and grad is large, 
    %   use true inv-Hess to update H as reasonable pos-def approx
    if status > 0 || steepest_status > 0 
        if norm(gradeq) > 0.1
          tmpclock = tic;
          H = set_H_to_inv_of_adjusted_true_hessian(nseg,zeq);
          hessupdate_time = hessupdate_time + toc(tmpclock);
        end
    end
  end
  all_time = toc(allclock);
  if diagflag > 0
    fprintf("Total time %15.10f\n",all_time);
    fprintf("  Gradient    %9.3f  Hessian %9.3f\n",grad_time,hessupdate_time);
    fprintf("  BFGSupdates %9.3f  fminbnd %9.3f  Wolfebisec %9.3f  Steepest %9.3f\n", ...
      bfgsupdates_time,fminbnd_time,wolfebisec_time,steep_time);
  end
return


function [a,znew,status] = steepest_descent_step(zstart,enstart,dir,diagflag)
  % Given current guess "zstart", with energy "enstart"
  %   and search direction "dir", use backtracking strategy to take
  %   steepest descent step that decreases energy by more than 10^(-14)
  %
  % If no such step can be found, declare that we are done by setting 
  %   output variable "status" to 0; otherwise, "status" is returned as 1
  %
  % Additional input variable "diagflag" controls what we print to screen
  a = 1; status = 0;
  while a > 10^(-13)
    ztest = zstart + a*dir;
    entest = discrete_dna_penalty_en(ztest);
    if diagflag > 1
      fprintf("a %20.15f enstart %20.15f entest %20.15f\n",a,enstart,entest)
    end
    if entest<enstart-10^(-14) % productive step, accept
      znew = ztest; status = 1;
      return
    end
    a=a/10;
  end
  znew = zstart; % We found no better z; status will return as 0
return

function [alpha,status] = wolfe_bisection(zstart,direc,diagflag)
% Use a bisection approach to try to satisfy Wolfe criteria 
%  starting at "zstart" with search direction "direc", i.e., considering
%  znew = zstart + alpha*direc

% Output alpha value, plus a flag "status" = -1 if failed, 1 if succeeded

% Additional input argument "diagflag" controls what we print to screen

  maxsteps = 15; % If we take more steps than this, declare failure
  steps = 0;

  [enstart,gradstart] = discrete_dna_penalty_en_grad_s(zstart);

  % [gamma,beta] is bracket for possible solution; alpha is current guess
  gamma = 0; alpha = 1; beta = 10^(4);

  % Compute two Wolfe quantities
  znew = zstart + alpha*direc;
  [ennew,gradnew] = discrete_dna_penalty_en_grad_s(znew);
  dotprod = direc'*gradstart; dotprodnew = direc'*gradnew;
  wolfe1=ennew-enstart+0.0001*alpha*dotprod;
  wolfe2 = -dotprodnew+0.9*dotprod;

  while (wolfe1 > 0 || wolfe2 > 0) && steps < maxsteps
    if wolfe1 >= 0 % If E increases or does not decrease enough
                   %  shift bracket to smaller values
      beta = alpha; alpha = (gamma+beta)/2;
    elseif wolfe2 >= 0 % If E decreases enough, but dot product decreases
                       %  or does not increase enough, 
                       %  shift bracket to larger values
      gamma = alpha;
      if beta > 5*10^3 
        alpha = 2*gamma;
      else
        alpha = (gamma+beta)/2;
      end
    else % If both Wolfe conditions hold, declare victory and exit function
      status = 1; 
      return
    end
    % Update everything and continue 
    znew = zstart + alpha*direc; steps=steps+1;
    [ennew,gradnew] = discrete_dna_penalty_en_grad_s(znew);
    dotprodnew = direc'*gradnew;
    wolfe1=ennew-enstart+0.0001*alpha*dotprod;
    wolfe2 = -dotprodnew+0.9*dotprod;
  end
  if steps >= maxsteps || abs(beta-gamma) < 0.001
      if diagflag > 1
          fprintf("Failed bisection after %d steps, alpha is %15.10f\n",steps,alpha);
          fprintf(" En is %15.10f, Wolfes are %15.10f %15.10f\n",ennew,wolfe1,wolfe2)
      end
      status = -1; % Reached end without achieving wolfe condition
  else
      status = 1;  % Satisfied both wolfe conditions
      if diagflag > 1
          fprintf("Finished bisection in %d steps, alpha is %15.10f\n",steps,alpha);
          fprintf(" En is %15.10f, Wolfes are %15.10f %15.10f\n",ennew,wolfe1,wolfe2)
      end
  end
return

function output = line_search_g(t)
% Function that we seek to minimize in line search
  global current_guess
  global direc
  output= discrete_dna_penalty_en(current_guess+t*direc);
return

function H = set_H_to_inv_of_adjusted_true_hessian(nseg,zeq)    
  % Given a z-vector "zeq" with "nseg" rod segments
  %   compute eigenvalue decomp of true Hessian at zeq, 
  %   then use that to find a positive-definite matrix
  %   "close" to (true Hessian)^(-1)
  [eneq,gradeq,hesseq] = discrete_dna_penalty_en_grad_hess(zeq);
  [V,D]=eigs(hesseq,7*(nseg-1));
  for i = 1:7*(nseg-1)
      if D(i,i) < 0.001
          D(i,i)=1;
      else
          D(i,i)=1/D(i,i);
      end
  end
  H = V*D*V';
return

function H = set_H_to_true_hessian_inv_if_posdef(nseg,zeq,currentH,diagflag)
  % Given a z-vector "zeq" with "nseg" rod segments
  %   compute eigenvalue decomp of true Hessian at zvec
  %   and if that is positive definite, return its inverse
  %   Otherwise, return current "approx of H^(-1)"
  [eneq,gradeq,hesseq] = discrete_dna_penalty_en_grad_hess(zeq);
  [V,D]=eigs(hesseq,7*(nseg-1));
  if min(diag(D)) > 0.001
    for i = 1:7*(nseg-1)
        D(i,i)=1/D(i,i);
    end
    if diagflag > 0
        fprintf("Setting H to true inv-Hess")
    end
    H = V*D*V';
  else
    H = currentH;
  end
return
