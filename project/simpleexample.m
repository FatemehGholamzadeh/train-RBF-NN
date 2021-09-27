close all
clc
clear

%%  Rosenbrock's function
    %   Minimum: f(1,1) = 0
%     f      = @(x,u) (1-x(1,:)).^2 + 100*(x(2,:)-x(1,:).^2).^2;
     f      = @(x,u) (atan(x(1,:))./cosh(x(2,:)));
%      f      = @(x,u) (1-x(1,:)).^2 + 100*(x(2,:)-x(1,:).^2).^2;
     f      = @(x,u) (cosh(x(1,:))./sinh(x(2,:)));
%      f      = @(x,u) (cosh(x(1,:))./sinh(x(2,:)));
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)
%%
%% Setting initial parameters
nf      = 1;                 % length of the output vector 'f(x,y)'
mu      = 100;               % parent population size
lambda  = 100;               % offspring population size
gen     = 20;               % number of generations
sel     = '+';               % Selection scheme (Pag. 78 in (BACK))
rec_obj = 2;                 % Type of recombination to use on object
                             % variables (Pag. 74 in (BACK))
                             % See 'recombination.m'
rec_str = 4;                 % Type of recombination to use on strategy
                             % parameters (Pag. 74 in (BACK))
u       = 0;                 % external excitation
%%
[min_x, min_f, off, EPS,idx] = evolution_strategy(f, mu, lambda, gen, sel, rec_obj, rec_str, u, obj, nf, n_x, limits);

