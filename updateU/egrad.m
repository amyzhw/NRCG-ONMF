function [egrad, f]= egrad(prob,x,cost_flag)
if nargin<3
    cost_flag = 0;
end
%GRAD   Computes the gradient on the manifold

egrad = x.U*(x.V'*x.V)-prob.D*x.V;
if cost_flag
    DTU = prob.D'*x.U;
    f = 0.5*prob.Dsqure - sum(sum(DTU.*x.V)) + 0.5*(1+prob.lamda)*norm(x.V,'fro')^2;
end


