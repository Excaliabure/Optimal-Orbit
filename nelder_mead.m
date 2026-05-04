% Nelder-Mead

% Nelder-Mead works by using n+1 points which are the vertices of a
% simplex,the simplest possible geometric shape in any dimension.
% All of the points are ranked until the worst , second worse, and best 
% candidates are identified. Next the centroid, or the average position of 
% all points in the current simplex except for the worst point, is also 
% identified. Then there are a series of transformations that can be 
% applied to the worst point such as reflection, expansion, contraction,
% and shrink. Reflection reflects the worst point across the centroid, if
% the reflected candidate is better than the original, we expand in that
% direction. Contraction is split into two cases, outside contraction occurs
% when the reflected point is better than the original worst candidate but 
% not the second worse we take a smaller step in the reflected direction.
% If the reflected candidate is worse than the original then we compute an 
% inside contraction in the opposite direction. Finally, shrink takes
% place when none of the prior steps have bettered the conditions, we
% shrink by moving in the direction of the best candidate. 



% reference : https://brandewinder.com/2022/03/31/breaking-down-Nelder-Mead/
%http://www.scholarpedia.org/article/Nelder-Mead_algorithm
function[x_opt, f_opt, iter] = nelder_mead(obj_fun, x0, t0, tf)

alpha = 1;
gamma = 2;
rho = 0.5;
sigma = 0.5;

tol = 1e-4;
maxiter = 60;
n = length(x0);

get_f = @(x_vec) - obj_fun(x_vec, t0, tf); % Convert to negative to convert from maximize f(x) to minimze -f(x)
step_size = [30; 60; 60; 60];

simplex = zeros(n, n+1);
simplex(:,1) = x0(:);

for i = 2:n+1
    y = x0(:);
    y(i-1) = y(i-1) + step_size(i-1);
    simplex(:,i) = y;
end 

iter = 0;

%Initializing fvals

fvals = zeros(1, n+1);
for i = 1:n+1
    fvals(i) = get_f(simplex(:,i));
end

% Sorting

while iter < maxiter
    [fvals, idx] = sort(fvals);
    simplex = simplex(:, idx);

    best = simplex(:,1);
    worst = simplex(:,end);

    if max(abs(fvals - fvals(1))) < tol
        break;
    end 

% Finding Centroid 

    centroid = mean(simplex(:,1:end-1), 2);

% Reflecting 

    xr = centroid + alpha*(centroid - worst);
    %xr = wrap_angles(xr);
    fr = get_f(xr);

   
% Expansion 

    if fr < fvals(1)
        xe = centroid + gamma*(xr - centroid);
        %xe = wrap_angles(xe);
        fe = get_f(xe);

        if fe < fr
            simplex(:,end) = xe;
            fvals(end) = fe;
        else 
            simplex(:,end) = xr;
            fvals(end) = fr;
        
        end 

    elseif fr < fvals(end-1)
        simplex(:,end) = xr;
        fvals(end) = fr;

% Contraction 

    else
        %xc = centroid + rho*(worst - centroid);
        %xc = wrap_angles(xc);
        %fc = get_f(xc);

        %if fc < fvals(end)
            %simplex(:,end) = xc;
            %fvals(end) = fc;
        %else 

        if fr < fvals(end)
            xc = centroid + rho*(xr- centroid);
        else
            xc = centroid + rho*(worst - centroid);
        end 
        fc = get_f(xc);

        if fc < fvals(end)
            simplex(:,end) = xc;
            fvals(end) = fc;
        else
% Shrink 

            for i = 2:n+1
                simplex(:,i) = best + sigma*(simplex(:,i) - best);
                %simplex(:,i) = wrap_angles(simplex(:,i));
                fvals(i) = get_f(simplex(:,i));
            end 
        end 
    end 
    iter = iter + 1;
    fprintf('Iterations %d: f(x) = %.6f\n', iter, fvals(1));
end 

x_opt = simplex(:,1);
f_opt = -fvals(1);

fprintf('\nFinal Solution:\n');
fprintf('Inclination: %.2f\n', x_opt(1));
fprintf('RAAN: %.2f\n', x_opt(2));
fprintf('Argument of Perigee (w): %.2f\n', x_opt(3));
fprintf('Initial True Anomaly (nu): %.2f\n', x_opt(4));
fprintf('Total SZA Time: %.2f seconds\n', f_opt);

end

function x_wrapped = wrap_angles(x)
x_wrapped = mod(x,360);
end 
