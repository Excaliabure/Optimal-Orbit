%{
@file    BOB1_Optimizer_SZA_constraints_simplex_newton_method_20260423.m
@author  Zander Abid / UCMerced (zabid@ucmerced.edu)
@date    2026-04-23
@brief   Main file for executing Bobsat-1 orbital optimization
 
@section This file is to apply both the quasai-newton and simplex (altered to Nedler-mead) 
algorithm
in attempting to find a feasable region of solutions for the orbital
parameters of bobsat-1
%}



clc; clear all; close all

% PARAMETERS
t0 = [2024,10,15,19,0,0]; % year, month, day, hour, minute, second
tf = 7; % Total amount of days elapsed till de-orbit, (for testing purposes keep low)
x0 = [64; % Inclination
      40; % RAAN
      0; % argument of perigee (w)
      30 % inital true anomoly (nu)
      ]';



% IC = Initial Conditions
% --- Optimizing (Engineering Method) ---

% optins for tuning minimizer 
% main_optimizer_options = optimoptions("fmincon", ... % 
%                                         "Algorithm","active-set");

% Solver
% [best_IC] = main_optimizer(@objective_fun, x0, t0, tf, main_optimizer_options);



% --- Optimizing (Class Implimentation)---

% nelder_mead(@objective_fun,x0,t0,tf)

% Simplex (needs to update to nedler-mead)
% simplexmethod([1],[1],[1])

% Quasai-Newton
% Set maxiter = 30, if want to change go into quasi_newton.m

% quasai_newton(@objective_fun,x0,t0,tf); % Displays optimal result after completion

% Quasi Newton Testing Log Reading
data = readmatrix('EXP_quasi_newton_timestart_20241015_7_days_20260501.txt');

% Nedler Mead Testing Log Reading
% data = readmatrix('EXP_nedler_mead_timestart_20241015_7_days_20260501.txt');

i = data(:,1);
RAAN = data(:,2);
w = data(:,3);
nu = data(:,4);

figure;
subplot(2,2,1)
plot(i)
title("Quasi: i")
xlabel("Objective Function Call")

subplot(2,2,2)
plot(RAAN)
title("Quasi: RAAN")
xlabel("Objective Function Call")

subplot(2,2,3)
plot(w)
title("Quasi: w")
xlabel("Objective Function Call")

subplot(2,2,4)
plot(nu)
title("Quasi: nu")
xlabel("Objective Function Call")

x = data;
fprintf('%.6f %.6f %.6f %.6f\n', ...
        x(1), x(2), x(3), x(4));

