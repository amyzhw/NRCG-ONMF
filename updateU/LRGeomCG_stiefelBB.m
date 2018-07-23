function [x, histout, itc, fail] = LRGeomCG_stiefelBB(prob, opts, x0)

addpath('.\manopt\tools');
if nargin<5
    rel_inner_tol = 1e-10;
end

if nargin<4
    fc0 = F_cost(prob,x0);
end

fail = true;

opts.orth_value = 0.20; % the search directions should be almost orthogonal
itc = 1;
xt = x0;
fc = fc0;
gc = grad(prob,xt);
ip_gc = inner_product(xt,gc,gc);
fold = fc; reschg = -1;
dir = scaleTxM_stiefel(gc,-1);

ithist=zeros(opts.maxit,5);
beta = 0;
rel_grad = 1;
ithist(itc,1) = rel_grad;
ithist(itc,2) = (fc);
ithist(itc,3) = reschg;
ithist(itc,4) = 0;
prob.t0 = clock;
ithist(itc,5) = etime(clock,prob.t0);
tinit = 10;


[n, k] = size(xt.U);
Q = 1; 
Cval = fc0;  
tau = 0.001;
opts.eta = 0.15;
eta = opts.eta;
gamma = 0.85;
beta = 1e-4;
congu  = 1;
for itc=2:opts.maxit
    
     fc_old = fc;
    
     nls = 1;
     qp0 = inner_product(xt,gc,dir);
     xc_new = xt;
     while 1 % wei: find stepsize
        % calculate G, F,
        xc_new.U = retraction_stiefel(xt.U,dir,tau);
        fc = F_cost(prob,xc_new);
%         out.nfe = out.nfe + 1;
        
        if (fc <= Cval + tau*beta*qp0) || nls >= 10
            break;
        end
        tau = eta*tau;          nls = nls+1;
    end  % end of while 1
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update V
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     if mod(itc,2)==0
    
    otx=prob.D'*xc_new.U;
    xc_new.V = updateV(otx,xc_new.V,prob.lamda + 1,0.01, 10); % use cpp update V

    [gc_new, fc] = grad(prob,xc_new, 1);

    
    rel_grad = 1;
    
    % for 'stagnation stopping criterion' after 10 iters
    reschg = abs(1-sqrt(2*fc)/sqrt(2*fold) );  % LMARank's detection
     if itc > 40 && reschg < opts.rel_tol_change_res
        if opts.verbosity > 0; disp('Iteration stagnated rel_tol_change_res.'); end
        fail = true;
        break;
    end
    if congu
         ip_gc_new = inner_product(xc_new,gc_new,gc_new);
        [dir, Y]= construct_search_direction(prob,xt,xc_new,gc,gc_new,gc_new,dir,ip_gc_new,ip_gc,opts);
     
    else
         ip_gc_new = 0;
         opts.beta_type = 'S-D';
        [dir, Y]= construct_search_direction(prob,xt,xc_new,gc,gc_new,gc_new,dir,ip_gc_new,ip_gc,opts);
    end
    
    S = xc_new.U - xt.U;    
    SY = abs(sum(sum(S.*Y)));
    if mod(itc-1,2)==0;
        tau = sum(sum(S.*S))/SY;
    else
        tau  = SY/sum(sum(Y.*Y));
    end
    tau = max(min(tau, 1e20), 1e-20); % wei:step size for next iter
    Qp = Q; Q = gamma*Qp + 1; Cval = (gamma*Qp*Cval + fc)/Q;
    
    
    % update _new to current
    gc = gc_new;
    ip_gc = ip_gc_new;
    xt = xc_new;
    fold = fc;
    
    ithist(itc,1) = rel_grad;
    ithist(itc,2) = (fc);
    ithist(itc,3) = reschg;
%     ithist(itc,4) = iarm;
    ithist(itc,5) = etime(clock,prob.t0);
    ithist(itc,6) = Re_Fnorm(prob, xc_new);

end

%/norm_M_Omega
x = xc_new;
ithist(itc,1) = rel_grad;
ithist(itc,2) = (fc);
ithist(itc,3) = reschg;
% ithist(itc,4) = iarm;
ithist(itc,5) = etime(clock,prob.t0);
ithist(itc,6) = Re_Fnorm(prob, x);
histout=ithist(1:itc,:);

iterno=itc;
