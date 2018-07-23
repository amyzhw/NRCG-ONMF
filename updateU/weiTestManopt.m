% function weiTestManopt()

clear all;
%-----initialize------
%%% set the size/sparsity of the data
% small
m = 100
n = 100;
k=5;

lamda=0.0001;
maxiteration = 50;

    
    %------prepare data---
    %%% generating X from random U and V
    U_org = orth(rand(m,k));
    %I = U_org'* U_org;
    V_org = abs(rand(n,k));
    X = abs(U_org * V_org');

    [U, S, V] = svds(X,k);
    x0.U = U;
    x0.V  = max(S*V',0)'; %

    %-----parameters and pre computed values------
    problem.lamda=lamda;
    problem.lstype =1; % NMF without closed form alpha 0: line_search_adaptive; 2: NMF with closed form for alpha
    problem.D = X;
    problem.Dsqure = norm(problem.D(:))^2;

    opts.maxit = maxiteration;
    opts.rel_inner_tol = 1e-5;
    opts.verbosity = 1;
    opts.rel_tol_change_res = 1e-5;
    
    
    %----Stiefel HS-------
    opts.beta_type='H-S';
    tic
    [x, histout, itc, fail] = LRGeomCG_stiefelBB(problem, opts, x0);
    toc

% ---- vs time       
  figure;

   semilogy(histout(:,5),histout(:,2), '-r','LineWidth',2); % object value vs time -dr
   set(gca,'FontSize',16);
   xlabel('Time (seconds)','FontSize',16)
   ylabel('Object Value','FontSize',16);
   
    
  %--------vs iteration  
   figure;  
   hold on
   semilogy(histout(:,2), '-r','LineWidth',2); %object value vs iteration
   
    set(gca,'FontSize',16);
   xlabel('Iteration','FontSize',16)
   ylabel('Object Value','FontSize',16);

          
% end
