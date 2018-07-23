function [grad, f]= grad(prob,x,cost_flag)
if nargin<3
    cost_flag = 0;
end

g = x.U*(x.V'*x.V)-prob.D*x.V; %gradient in vector space %+prob.lamba*x.U
grad = projection_tangent(x.U, g);

if cost_flag
    DTU = prob.D'*x.U;
    f = 0.5*prob.Dsqure - sum(sum(DTU.*x.V)) + 0.5*(1+prob.lamda)*norm(x.V,'fro')^2;
end
% toc
