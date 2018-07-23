function [xt,ft, succ,numf,iarm,t] = armijo_search_stiefel(prob, xc,fc,gc,dir,t)

% Armijo line search with polynomial interpolation
% Adapted from steep by C. T. Kelley, Dec 20, 1996

% xc, current point
% fc = F(sys,xc)
% gc is gradient
% dir is search direction
%
% returns: xt,ft: Armijo point and F(sys,xt)
%          succ = true if success
%          numf: nb of F evals
%          iarm: nb of Armijo backtrackings




%%% constants
alp=1e-4; blow=0.1; bhigh=0.5;
%alp = 1e-8; bhigh=.5; blow=.1;
%%%
MAX_IARM = 20;

numf = 0;

% trial evaluation at full step
xt.U = retraction_stiefel(xc.U,dir,t);
xt.V = xc.V;
ft = F_cost(prob,xt);
numf = numf+1; iarm = 0;

%fgoal = fc-alp*lambda*ip(xc,gc,gc); %norm(gc.Y,2)^2;
%
%       polynomial line search
%
q0=fc; qp0=inner_product(xc,gc,dir); lamc=t; qc=ft;
fgoal=fc+alp*t*qp0;

while(ft > fgoal)    
    iarm=iarm+1;
    if iarm==1
        t=polymod(q0, qp0, lamc, qc, blow, bhigh);
    else
        t=polymod(q0, qp0, lamc, qc, blow, bhigh, lamm, qm);
    end
    qm=qc; lamm=lamc; lamc=t;
    xt.U = retraction_stiefel(xc.U,dir,t);
    xt.V = xc.V;
    ft = F_cost(prob,xt);
    numf = numf+1; qc=ft;
    if(iarm > MAX_IARM)
        succ = false;        
        return
    end    
    fgoal=fc+alp*t*qp0;
end
succ = true;


end
