function [desc_dir_new, diff]= construct_search_direction(prob,xc,xc_new,gc,gc_new,Pgc_new,desc_dir_old,gc_newPgc_new,gradPgrad,options)

if strcmpi(options.beta_type, 'steep') || ...
        strcmpi(options.beta_type, 'S-D')              % Gradient Descent
    
    beta = 0;
    desc_dir_new = scaleTxM_stiefel(gc_new, -1);
    desc_dir = projection_tangent(xc_new.U, desc_dir_old);
    diff = matrixlincomb(xc_new.U, 1, desc_dir, -1, desc_dir_new);
    %diff = matrixlincomb(xc_new.U, 1, gc_new, -1, oldgrad);%
else
    
    oldgrad = projection_tangent(xc_new.U, gc);
    orth_grads = inner_product(xc_new.U, oldgrad, Pgc_new) / gc_newPgc_new;
    
    % Powell's restart strategy (see page 12 of Hager and Zhang's
    % survey on conjugate gradient methods, for example)
    if abs(orth_grads) >= options.orth_value,
        beta = 0;
        desc_dir_new = scaleTxM_stiefel(gc_new, -1);
        diff = matrixlincomb(xc_new.U, 1, gc_new, -1, oldgrad);
    else % Compute the CG modification
        
        desc_dir = projection_tangent(xc_new.U, desc_dir_old); % wei ?? old desc_dir!
        
        switch upper(options.beta_type)
            
            case 'F-R'  % Fletcher-Reeves
                beta = gc_newPgc_new / gradPgrad;
                diff = matrixlincomb(xc_new.U, 1, gc_new, -1, oldgrad);
            case 'P-R'  % Polak-Ribiere+
                % vector grad(new) - transported grad(current)
                diff = matrixlincomb(xc_new.U, 1, gc_new, -1, oldgrad);
                ip_diff = inner_product(xc_new.U, Pgc_new, diff);
                beta = ip_diff / gradPgrad;
                beta = max(0, beta);
                
            case 'H-S'  % Hestenes-Stiefel+
                diff = matrixlincomb(xc_new.U, 1, gc_new, -1, oldgrad);
                ip_diff = inner_product(xc_new.U, Pgc_new, diff);
                beta = ip_diff / inner_product(xc_new.U, diff, desc_dir);
                beta = max(0, beta);
                
            case 'H-Z' % Hager-Zhang+
                diff = matrixlincomb(xc_new.U, 1, gc_new, -1, oldgrad);
                Poldgrad = projection_tangent(xc_new.U, Pgrad);
                Pdiff = matrixlincomb(xc_new.U, 1, Pgc_new, -1, Poldgrad);
                deno = inner_product(xc_new.U, diff, desc_dir);
                numo = inner_product(xc_new.U, diff, Pgc_new);
                numo = numo - 2*inner_product(xc_new.U, diff, Pdiff)*...
                    inner_product(xc_new.U, desc_dir, gc_new) / deno;
                beta = numo / deno;
                
                % Robustness (see Hager-Zhang paper mentioned above)
                desc_dir_norm = norm(desc_dir(:));
                eta_HZ = -1 / ( desc_dir_norm * min(0.01, gradnorm) );
                beta = max(beta, eta_HZ);
                
            otherwise
                error(['Unknown options.beta_type. ' ...
                    'Should be steep, S-D, F-R, P-R, H-S or H-Z.']);
        end
        
        desc_dir_new = matrixlincomb(xc_new.U, -1, Pgc_new, beta, desc_dir);
%         diff = matrixlincomb(xc_new.U, 1, desc_dir, -1, desc_dir_new);
    end
    
end