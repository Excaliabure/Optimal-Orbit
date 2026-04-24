%{
@file    quasai_newton.m
@author  Armando Mendoze Cruz/UCMerced (email)
@date    2026-04-23
@brief   Simple orbital propegation solver with quasai-newton method
 
@section This program was made to attempt to optimize a satellite (Low
earth orbit) orbit trajectory given the solar zenith constraint of 30 < SZA
< 50. 
%}



function [x_opt, f_opt, iter] = quasai_newton(obj_fun, x0, t0, tf)
    %{
    Params: 
        obj_fun - function handle (use @obj_fun when calling method)
        x0 - [Inclination_0
             R.A.A.N_0
             Argument of Perigee (w)_0
             Initial true anomaly_0
            ]
        
        t0 - [year, month, day, hour, minute, second]
        tf - # of days elapsed till end
    Return Val:
        x_opt - optimal x values per iteration
        f_opt - optimal total function values
        iter  - iteratin number on 
    %}

    tol = 1e-4;
    maxiter = 30; 
    h = 1; % Step size for numerical gradient (larger for satellite scenarios)
    c1 = 1e-1;
    
    n_vars = length(x0);
    x_k = x0(:);
    H_inv = eye(n_vars);
    
    % Internal wrapper to handle the 3-argument objective_fun and negate it
    % This ensures t0 and tf are passed correctly every time
    get_f = @(x_vec) -obj_fun(x_vec, t0, tf);
    
    % Numerical Gradient Function (Central Difference)
    function g = compute_grad(x)
        g = zeros(size(x));
        for j = 1:length(x)
            x_plus = x; x_plus(j) = x_plus(j) + h;
            x_minus = x; x_minus(j) = x_minus(j) - h;
            g(j) = (get_f(x_plus) - get_f(x_minus)) / (2 * h);
        end
    end

    f_k = get_f(x_k);
    grad_k = compute_grad(x_k);
    
    iter = 0;
    while norm(grad_k) > tol && iter < maxiter
        p_k = -H_inv * grad_k;
        
        % Simple Backtracking Line Search
        alpha = 1;
        while get_f(x_k + alpha * p_k) > f_k + c1 * alpha * (grad_k' * p_k)
            alpha = alpha * 0.5;
            if alpha < 1e-5, break; end
        end
        
        x_next = x_k + alpha * p_k;
        f_next = get_f(x_next);
        grad_next = compute_grad(x_next);
        
        % BFGS Update
        s_k = x_next - x_k;
        y_k = grad_next - grad_k;
        if y_k' * s_k > 1e-10
            rho = 1 / (y_k' * s_k);
            H_inv = (eye(n_vars) - rho*s_k*y_k') * H_inv * (eye(n_vars) - rho*y_k*s_k') + rho*(s_k*s_k');
        end
        
        x_k = x_next;
        f_k = f_next;
        grad_k = grad_next;
        iter = iter + 1;
        fprintf('Iteration %d: Time = %.2f seconds\n', iter, -f_k);
       
    end
    
    x_opt = x_k;
    f_opt = -f_k;
    fprintf('Inclination %.2f\nRAAN : %.2f\nArgument of Perigee (w) : %.2f\nInital true anomaly (nu) : %.2f', x_opt(1),x_opt(2),x_opt(3),x_opt(4));
end


