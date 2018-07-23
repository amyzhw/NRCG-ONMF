function [ren]= Re_Fnorm(prob,x)

DTU = prob.D'*x.U;
ren = 0.5*prob.Dsqure - sum(sum(DTU.*x.V)) + 0.5*norm(x.V,'fro')^2;
ren= ren/(0.5*norm(prob.D,'fro')^2);
