function [min_x, min_f, off, EPS,j] = evolution_strategy(mu, lambda, gen, sel, rec_obj, rec_str, u, obj, nf, n,X,Ystar,limits, x_0, sigma)
%%  [min_x, min_f, off, EPS,j] = evolution_strategy(f, mu, lambda, gen, sel, rec_obj, rec_str, u, obj, nf, n, limits, x_0, sigma)
%
%   This function tries to optimize the function 'f' using "evolution
%   strategies". The arguments "x_0" and "sigma" are optional.
%   Taken from "Algorithm 4", Section 2.1.6, pag. 81 in (BACK)
%
%   INPUT DATA:
%
%   - f:       Objective function (handle function: f(x,u))
%   - mu:      Parent population size (positive integer number)
%   - lambda:  Offspring population size (positive integer number)
%   - gen:     Number of generations (positive integer number)
%   - sel:     Selection scheme (Pag. 78 in (BACK)):
%                  * ',' = (mu, lambda)-selection scheme
%                  * '+' = (mu + lambda)-selection scheme
%   - rec_obj: Type of recombination to use on objective variables
%              (Pag. 74 in (BACK)):
%                  * 1   = No recombination
%                  * 2   = Discrete recombination
%                  * 3   = Panmictic discrete recombination
%                  * 4   = Intermediate recombination
%                  * 5   = Panmictic intermediate recombination
%                  * 6   = Generalized intermediate recombination
%                  * 7   = Panmictic generalized intermediate recombination
%   - rec_str: Type of recombination to use on strategy parameters
%              (Pag. 74 in (BACK)).
%   - u:       External excitation (if it does not exist, type 0 (zero))
%   - obj:     Vector with the desired results
%   - nf:      Length of the handle function vector (length(f) x 1 vector)
%   - n:       Length of the vector x_0 (positive integer number)
%   - limits:  Matrix with the limits of the variables (nx2 matrix). The
%              first column is the lower boundary, the second column is
%              the upper boundary.
%   - x_0:     Starting point (optional) (nxmu matrix)
%   - sigma:   Cell with covariance matrices (Optional) (1 x mu cell; each
%              cell has to have an nxn symmetric matrix)
%
%   OUTPUT DATA:
%
%   - min_x:   Cell with the parent population, and whose last component
%              minimizes the objective function 'f'
%              vector)
%   - min_f:   Cell with the values of the Objective Function 'f'
%              (length(f) x 1 vector)
%   - off:     Cell with the offspring population in each generation
%   - EPS:     Vector with the minimum error of each generation (gen x 1
%              vector)
%   - j:       Number of iterations the algorithm ran (Final number of
%              generations)
%
%   Bibliography:
%
%   - BACK, Thomas. "Evolutionary algorithms in theory and practice".
%     Oxford University Press. New York. 1996.
%
% -------------------------------------------------------
% | Developed by:   Gilberto Alejandro Ortiz Garcia     |
% |                 gialorga@gmail.com                  |
% |                 Universidad Nacional de Colombia    |
% |                 Manizales, Colombia.                |
% -------------------------------------------------------
%
%   Date: 20-Sep-2011

%% Beggining

%% Initialization:
if ((sel ~= ',') || (sel ~= '+'))
  error('not supported selection scheme')
end

if exist('x_0','var')
  if exist('sigma','var')
    [x_0, sigma, alpha] = validate_sizes(mu, n, limits, x_0, sigma);
  else
    [x_0, sigma, alpha] = validate_sizes(mu, n, limits, x_0);
  end
else
  [x_0, sigma, alpha] = validate_sizes(mu, n, limits);
end

e     = 1e-8;                         % maximum permissible error
min_x = cell(1,gen);                  % allocate space in memory for min_x
min_f = cell(1,gen);                  % allocate space in memory for min_f
off   = cell(1,gen);                  % allocate space to store the offspring population

min_x{1} = x_0;                       % first point
value = zeros(nf,mu);                 % allocate space for function evaluation
for i = 1:mu
  value(:,i) = fun(X,Ystar,x_0(:,i));
end
min_f{1} = value;                     % first approximation
off{1}   = zeros(n,1);

j      = 1;                           % generations counter
jj     = 0;                           % stagnation counter
eps    = abs(obj - value(1,:));       % initial error
EPS    = zeros(gen,1);                % allocate space in memory for minimum error every generation
EPS(1) = min(eps);

%if n == 1
  %% Plot function
  %figure
  %X = linspace(limits(1,1),limits(1,2),100);
  %Y = f(X,[]);
  %plot(X,Y);
  %hold on;
%elseif n == 2
  %% Plot function
  %figure
  %[X,Y] = meshgrid(linspace(limits(1,1),limits(1,2),100),linspace(limits(2,1),limits(2,2),100));
  %Z = reshape(f([X(:) Y(:)]',[]), 100, 100);
  %contour(X,Y,Z,30,'k');      % Contour plot
  %grid on
  %xlabel('x','FontSize',16);
  %ylabel('f(x)','FontSize',16);
  %title('Offspring evolution','FontSize',18);
  %hold on;
  %pcolor(X,Y,Z);              % is really a SURF with its view set to directly above
  %shading interp
%end

%% Begin ES
while ((j < gen) && (min(eps) > e))
  %% Print report:
  if mod(j,5) == 0
    fprintf('\tGeneration j = %4d,  fitness = %g\n',j,min(eps));
  end;
  
  %% Recombine:
  [xr,sigmar,alphar] = recombination(rec_obj,rec_str,n,mu,lambda,min_x{j},sigma,alpha);
  off{j+1}           = xr;            % offspring population
  
  %% Mutation:
  [xm,sigmam,alpham] = mutation(n,lambda,xr,sigmar,alphar,limits);
  
  %% Evaluation:
  phie = zeros(nf,lambda);
  for i = 1:lambda
    phie(:,i) = fun(X,Ystar,xm(:,i));
  end
  epse = abs(obj - phie(1,:));
  
  %% Selection:
  [min_x{j+1}, sigma, alpha] = selection(sel, mu, lambda, epse, eps, xm, sigmam, min_x{j}, sigma, alpham, alpha);

  %% Store better results:
  value = zeros(nf,mu);               % allocate space for function evaluation
  for i = 1:mu
    value(:,i) = fun(X,Ystar,min_x{j+1}(:,i));
  end
  min_f{j+1} = value;                 % next approximation
  eps = abs(obj - value(1,:));        % error
  
  EPS(j+1) = min(eps);
  
  %% Stagnation criterion:
  if (EPS(j) == EPS(j+1))
    jj = jj+1;
  else
    jj = 0;
  end
  
  %% Increase generation counter:
  j = j+1;
  
  %% Plot preliminary results
  %if n == 1
    %% Plot offspring population at each generation
   % h = plot(off{j-1}(1,:),min_f{j-1}(1,:),'*r');
   % pause(0.2)
  %  if (((EPS(j) > e) && (j < gen)) && (jj<30))
  %    delete(h);
   % end
 % elseif n == 2
    %% Plot offspring population at each generation
 %   h = plot(off{j-1}(1,:),off{j-1}(2,:),'*r');
 %   axis([limits(1,1) limits(1,2) limits(2,1) limits(2,2)]);
    %pause(0.2)
    %if (((EPS(j) > e) && (j < gen)) && (jj<30))
      %delete(h);
    %end
 % end

  if (jj == 30)
    fprintf('\n\n\tError remains constant for 30 consecutive generations\n\n');
    break
  end

end


end
%% END