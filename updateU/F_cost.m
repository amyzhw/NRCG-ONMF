function [f]= F_cost(prob,x)

DTU = prob.D'*x.U;
% tic
f = 0.5*prob.Dsqure - sum(sum(DTU.*x.V)) + 0.5*(1+prob.lamda)*norm(x.V,'fro')^2;
