% Function finds the best analysis operator Omega that lies on the oblique
% manifold via a riemannian conjugate gradient method. The cost function to
% be minimized is f(Omega) = ||Omega*Y||_p - kappa*log(det(Omega'*Omega))
% (c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
% Muenchen, 2012. Contact: simon.hawe@tum.de
function para = goal(Y, para)
% Normalize the initial Operator
Omega       = bsxfun(@times,para.Omega,1./sqrt(sum(para.Omega.^2,2)));
nu 			= para.nu;
q 			= [para.q,para.p];
Sp_type     = para.Sp_type;
kappa       = 1/(log(size(para.Omega,2))*size(para.Omega,2))*para.kappa;
mu          = para.mu*2/(size(para.Omega,1)^2-size(para.Omega,1));
M           = numel(Y);%size(Y,2);
% Linesearch parameters
alpha 		= 1e-1;
beta  		=   .9;
Qk          =    0;
Ck          =    0;
nuk         =   .0;

SUB_LOG_DET = size(Omega,2)*log(size(Omega,1));

I = eye(length(Omega)) + ones(length(Omega));
I_ext_grad = I == 2;
OY = Omega*Y;
[f0,q_w] = Sparsifying_functions(Sp_type, 'Evaluate', OY, q, nu);

if para.verbose
    Sp_input = f0;
end
f0 = 1/M*f0;

if mu ~= 0
    k_term = - mu*sum(sum(log(I-(Omega*Omega').^2)));
    f0           =  f0 + k_term;
end

if kappa ~= 0
    [~,Si,Vi]       = svd(Omega);
    if para.zmean
        Vi = Vi(:,1:end-1);           % rank n-1
        Si = Si(1:end-1,1:end-1);     % rank n-1
    end
    LogDetTerm      =  log(prod(diag(Si).^2))-SUB_LOG_DET;
    g_term          = -kappa*(LogDetTerm);
    f0              =  f0 + g_term;
end


for k = 1:para.max_iter
    if para.verbose && (k==1 || mod(k,10) == 0)
        draw_atoms(Omega,[para.p_sz,para.p_sz],para.p_sz*3);
    end
    %% Computing the gradient
    Omega_grad = Sparsifying_functions(Sp_type, 'Derivative', OY, q, nu, [], q_w);
    
    Omega_grad  = 1/M*Omega_grad*Y';
    if mu ~= 0
        Mat = (Omega*Omega');
        Mat(I_ext_grad) = 0;
        Omega_grad  = Omega_grad  + 4*mu*((Mat./(I-Mat.^2))*Omega);
    end
    if kappa ~= 0
        Omega_grad  = Omega_grad - 2*kappa*Omega*Vi*diag(1./diag(Si).^2)*Vi';
    end
            
    %% Transpose for update step on the Oblique manifold
    Omega_grad = Omega_grad';
    Omega      = Omega';
    
    if para.zmean   % project operator atoms onto 1-normal plane
        Omega_grad = bsxfun(@minus, Omega_grad, mean(Omega_grad, 1));
    end
    
    % Projection of the gradient onto the tangent space via
    % dO = dO - O*ddiag(O'*dO);
    Omega_grad = Omega_grad - bsxfun(@times,Omega,sum(Omega.*Omega_grad));
    % Omega_grad(:,1)=0;
    if numel(Omega_grad(isnan(Omega_grad)))
        Omega_grad(isnan(Omega_grad)) = 1e10;
        %return
    end
    
    if para.zmean   % project operator atoms onto 1-normal plane (again for numerical reasons)
        Omega_grad = bsxfun(@minus, Omega_grad, mean(Omega_grad, 1));
    end
    
    %% Computation of Conjugate Direction
    if k == 1 || ~mod(k,numel(Omega))
        dx          = -Omega_grad;
        t = 2*pi/max(sqrt(sum(dx.^2)));
        t_initial   = t;
    else
        tau_g  =  parallel_transport( g0, Omega_old, dx, t_prev, Norm_dx );
        tau_dx =  parallel_transport([], Omega_old, dx, t_prev, Norm_dx);
        
        yk     =  Omega_grad - tau_g;
        denom  = tau_dx(:)'*yk(:);
        
        cg_beta     = (Omega_grad(:)'*yk(:))/denom;
        cg_beta_dy  = norm(Omega_grad(:))^2/denom;
        cg_beta     = max(0,min(cg_beta_dy, cg_beta));
        
        dx       = -Omega_grad + cg_beta*tau_dx;
    end
    
    val      = Omega_grad(:)'*dx(:);
    Norm_dx  = sqrt(sum(dx.^2));
    sel      = Norm_dx > 0;
    ls_iter  = 0;
    
    Ck = (nuk*Qk*Ck+f0)/(nuk*Qk+1);
    Qk =  nuk*Qk+1;
    

    %% Backtracking Linesearch
    while (ls_iter == 0) || ... % Quasi Do while loop
            (f0 > Ck + alpha * t * val) && ... % Check Wolfe Condition
            (ls_iter < 100) % Check Maximum number of Iterations
        
        % Update the step size
        t  = t*beta;
        % Store the step length as it is required for the parallel
        % transport and the retraction
        t_prev    = t;
        
        % Exponential mapping and transpose back to get standing operator
        Omega_c = exp_mapping(Omega, dx, t, Norm_dx,sel)';
        % Evalute the costfunction
        OY = Omega_c*Y;
        [f0_sp,q_w] = Sparsifying_functions(Sp_type, 'Evaluate', OY, q, nu);
        f0         = 1/M*f0_sp;
        if mu ~= 0
            k_term = - mu*sum(sum(log(I-(Omega_c*Omega_c').^2)));
            f0         = f0 + k_term;
        end
        
        if kappa ~= 0
            [~,Si,Vi]      = svd(Omega_c);
            if para.zmean
                Vi = Vi(:,1:end-1);           % rank n-1
                Si = Si(1:end-1,1:end-1);     % rank n-1
            end
            g_term      =  -kappa*(log(prod(diag(Si).^2))-SUB_LOG_DET);
            f0         = f0 + g_term;
        end

        % Increase the number of linesearch iterates as we only
        % allow a certain amount of steps
        ls_iter = ls_iter + 1;
    end
    para.logger = [para.logger,f0];
    
    Omega = Omega';
    
    if ~mod(k,para.verbose)
        fprintf('Change of the operator: %f\n',norm(Omega(:)-Omega_c(:)))
        fprintf('Current Sparsity %e ~ Input Sparsity %e\n', f0_sp, Sp_input)
        fprintf('Number of Linesearchsteps %d ~ Stepsize: %e ~estimate %e\n',ls_iter,t_prev,t_initial)
        OO=Omega_c*(eye(size(Omega_c,2))-1/size(Omega_c,2)*ones(size(Omega_c,2)));
        Sings = svd(OO);
        Cohe = sort(abs(OO*OO'));
        fprintf('Condition: %f ~ Mutual Coherence: %f\n',Sings(1)/Sings(end-1),max(Cohe(end-1,:)))
        fprintf('Mean Correlation: %f ~ Variance Correlation: %f\n',mean(mean(Cohe(1:end-1,:))),var(var(Cohe(1:end-1,:))))
    end
    
    if  norm(Omega-Omega_c,'fro') < 1e-6
        fprintf('************ NO PROGRESS **********\n');
        break;
    end
    
    para.Omega = Omega_c;
    
    if ls_iter == 100
        fprintf('************ FINISHED **********\n');
        break;
    else
        Omega_old = Omega';
        g0      = Omega_grad;
        Omega   = Omega_c;
        t       = t/beta^2;
    end
    if ~mod(k,para.verbose)
        fprintf('************ NEXT ITERATION %d **********\n',k);
    end
end


end





