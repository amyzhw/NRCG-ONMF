function [x, histout, itc, fail] = LRGeomCG(prob, opts, x0)
% LRGEOMCG    Low-rank NMF by Geometric CG.
%   Uses Armijo rule, polynomial linesearch on embedded submanifold of
%   fixed rank matrices.
%
% Input: prob     = problem instance, see MAKE_PROB.
%        opts     = options, see DEFAULT_OPTS.
%        x0       = starting guess.
%
% Output: x       = solution.
%         histout = iteration history. Each row of histout is
%                   [rel_norm(grad), rel_err_on_omega, relative_change, ...
%                        number of step length reductions, restarts]
%         fail    = flag for failure
%
% See also SMALL_EXAMPLE, MAKE_PROB, DEFAULT_OPTS.

% (C) Bart Vandereycken, 2011-2012.
% Adapted from steep by C. T. Kelley, Dec 20, 1996.
%
% Paper: Low-rank matrix completion by Riemannian optimization.
%
% Info and bugs: bart.vandereycken@epfl.ch

addpath('.\manopt\tools');
if nargin<5
    rel_inner_tol = 1e-10;
end

if nargin<4
    fc0 = F_cost(prob,x0);
end

%beta_type = 'F-R';
% beta_type = 'P-R';

fail = true;


opts.orth_value = 0.1; % the search directions should be almost orthogonal
opts.beta_type = 'H-S';
itc = 1;
xc = x0;
fc = fc0;
gc = grad(prob,xc);
ip_gc = inner_product(xc,gc,gc);
fold = fc; reschg = -1;
% first search-dir is steepest gradient
dir = scaleTxM_stiefel(gc,-1);
% rel_grad = sqrt(ip_gc)/max(1,norm(gc,'fro'));

ithist=zeros(opts.maxit,5);
beta = 0;
rel_grad = 1;
ithist(itc,1) = rel_grad;
ithist(itc,2) = (fc);%/norm_M_Omega
ithist(itc,3) = reschg;
ithist(itc,4) = 0;
prob.t0 = clock;
ithist(itc,5) = etime(clock,prob.t0);
tinit = 1;
congu  = 1;
for itc=2:opts.maxit
    %   tinit = exact_search_onlyTxM(prob, xc,dir);
    fc_b = fc;
    [xc_new,fc,succ,~,iarm,tc] = armijo_search_stiefel(prob, xc, fc, gc, dir, tinit);
    tinit = 5*tc;
    if ~succ && beta ~= 0
        if opts.verbosity > 0; disp('Line search failed on CG. Resetting to gradient.'); end
        beta = 0;
        % if we did cg, reset to steepest descent and try again
        dir = scaleTxM_stiefel(gc,-1);
        %     tinit(itc) = exact_search_onlyTxM(prob, xc,dir);
        [xc_new,fc,succ,~,iarm, tc] = armijo_search_stiefel(prob, xc,fc_b,gc,dir, tinit);
        tinit = 5*tc;
    end
    
    
    
    % if it still fails (beta is always 0 here -> steepest descent)
    % then we give up
    if ~succ
        x = xc_new;
        ithist(itc,1) = rel_grad;
        ithist(itc,2) = (fc);%/norm_M_Omega
        ithist(itc,3) = reschg;
        ithist(itc,4) = iarm;
        ithist(itc,5) = etime(clock,prob.t0);
        histout=ithist(1:itc,:);
        if opts.verbosity > 0; disp('Line search failed on steepest descent. Exiting...'); end
        return
    end
    
    
    %update V
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update V
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if mod(itc,2)==0
        otx=prob.D'*xc_new.U;
        xc_new.V = updateV(otx,xc_new.V,prob.lamda + 1,0.01, 10); % use cpp update V
%     end % end if rem
    
    % grad(new x)
    gc_new = grad(prob,xc_new);
    ip_gc_new = inner_product(xc_new,gc_new,gc_new);
    
    
    % Test for stopping
    if (itc>3)&&((ithist(itc-2,2)-ithist(itc-1,2))/ithist(itc-2,2) < rel_inner_tol)
        break;
    end
    
    
    rel_grad = 1;
    %   if rel_grad < opts.rel_grad_tol
    %     if opts.verbosity > 0; disp('Relative gradient tol reached.'); end
    %     fail = false;
    %     break;
    %   end
    
    % for 'stagnation stopping criterion' after 10 iters
    reschg = abs(1-sqrt(2*fc)/sqrt(2*fold) );  % LMARank's detection
    %reschg = abs(sqrt(2*fc) - sqrt(2*fold)) / max(1,norm_M_Omega);
    if itc > 10 && reschg < opts.rel_tol_change_res
        if opts.verbosity > 0; disp('Iteration stagnated rel_tol_change_res.'); end
        fail = true;
        break;
    end
    if congu
        
        dir = construct_search_direction(prob,xc,xc_new,gc,gc_new,gc_new,dir,ip_gc_new,ip_gc,opts);
%         % new dir = -grad(new x)
%         %           + beta * vectTransp(old x, old dir, tmin*old dir)
%         gc_old = transpVect(prob,xc,gc,xc_new,1);
%         dir = transpVect(prob,xc,dir,xc_new,1);
%         
%         % we test how orthogonal the previous gradient is with
%         % the current gradient, and possibly reset the to the gradient
%         ip_gc_old_new = ip(xc_new,gc_old,gc_new);
%         orth_grads = ip_gc_old_new / ip_gc_new;
%         
%         if abs(orth_grads) >= ORTH_VALUE
%             if opts.verbosity > 1; disp('New gradient is almost orthogonal to current gradient. This is good, so we can reset to steepest descent.'); end
%             beta = 0;
%             dir = plusTxM(gc_new, dir, -1, beta);
%             
%         else % Compute the CG modification
%             % beta
%             if strcmp(beta_type, 'F-R')  % Fletcher-Reeves
%                 beta = ip_gc_new / ip_gc;
%                 % new dir
%                 dir = plusTxM(gc_new, dir, -1, beta);
%             elseif strcmp(beta_type, 'P-R')  % Polak-Ribiere+
%                 % vector grad(new) - transported grad(current)
%                 beta = (ip_gc_new - ip_gc_old_new) / ip_gc;
%                 beta = max(0,beta);
%                 dir = plusTxM(gc_new, dir, -1, beta);
%             end
%         end
%         
%         % check if dir is descent, if not take -gradient (i.e. beta=0)
%         g0 = ip(xc_new,gc_new,dir);
%         if g0>=0
%             if opts.verbosity > 1;
%                 disp('New search direction not a descent direction. Resetting to -grad.');
%             end
%             dir = scaleTxM(gc_new, -1);
%             beta = 0;
%         end
        
    else
        dir = scaleTxM_stiefel(gc_new, -1);
    end
    
    % update _new to current
    gc = gc_new;
    ip_gc = ip_gc_new;
    xc = xc_new;
    fold = fc;
    
    ithist(itc,1) = rel_grad;
    ithist(itc,2) = (fc);
    ithist(itc,3) = reschg;
    ithist(itc,4) = iarm;
    ithist(itc,5) = etime(clock,prob.t0);
end
%/norm_M_Omega
x = xc_new;
ithist(itc,1) = rel_grad;
ithist(itc,2) = (fc);
ithist(itc,3) = reschg;
ithist(itc,4) = iarm;
ithist(itc,5) = etime(clock,prob.t0);
histout=ithist(1:itc,:);

