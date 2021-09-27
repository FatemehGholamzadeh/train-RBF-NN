function [xm,sigmam,alpham] = mutation(n,lambda,xr,sigmar,alphar,limits)
%%  [xm,sigmam,alpham] = mutation(n,lambda,xr,sigmar,alphar,limits)
%   This function performs the mutation required in evolutionary strategies
%   to obtain the next offspring population.
%
%   INPUT DATA:
%
%   - n: Dimension of the state vector 'x_0' (positive integer number)
%   - lambda: Offspring population size (positive integer number)
%   - xr:     Recombined state (n x lambda vector)
%   - sigmap: Recombined covariances (1 x lambda cell; each cell has an
%             n x n matrix)
%   - alphap: Recombined rotation angles (1 x lambda cell; each cell has an
%             n x n matrix)
%   - limits: Matrix with the limits of the variables (nx2 matrix). The
%             first column is the lower boundary, the second column is
%             the upper boundary.
%
%   OUTPUT DATA:
%
%   - xm:     Mutated state
%   - sigmam: Mutated covariances
%   - alpham: Mutated rotation angles
%
%   Bibliography:
%
%   - BACK, Thomas. "Evolutionary algorithms in theory and practice". Oxford
%     University Press. New York. 1996.
%
% -------------------------------------------------------
% | Developed by:   Gilberto Alejandro Ortiz Garcia     |
% |                 gialorga@gmail.com                  |
% |                 Universidad Nacional de Colombia    |
% |                 Manizales, Colombia.                |
% -------------------------------------------------------
%
%   Date: 20 - Sep - 2011

%% Beginning:

%% Mutation factors (Eq. (2.18) in BACK):
tau   = 1/(sqrt(2*sqrt(n)));          % learning rate
taup  = 1/(sqrt(2*n));                % learning rate
beta  = 5*pi/180;                     % 5 degrees (in radians)

%% Mutate:
xm     = zeros(n,lambda);
sigmam = cell(1,lambda);
alpham = cell(1,lambda);

for i = 1:lambda
  tmp       = randn(n,n);
  sigmam{i} = sigmar{i}.*exp(taup*randn + tau*(tmp + tmp'));
  tmp       = rand(n,n);
  alpham{i} = alphar{i} + beta*triu((tmp + tmp'),1);
  
  %% Coordinate transformation with respect to axes 'i' and 'j' and angle
  %  'alpha_ij' (Eq. 2.14)
  R = eye(n);
  for m = 1:n-1
    for q = m+1:n
      T               =  eye(n);
      T([m q], [m q]) =  [  cos(alpham{i}(m,m))     -sin(alpham{i}(m,q))
                            sin(alpham{i}(q,m))      cos(alpham{i}(q,q)) ];
      R               =  R*T;
    end
  end

  xm(:,i) = xr(:,i) + R*sqrt(diag(diag(sigmam{i})))*randn(n,1);
  
  %% take in account boundaries (limits)
  for ii = 1:n
    % Lower boundary
    if xm(ii,i) < limits(ii,1)
      xm(ii,i) = limits(ii,1);
    end
    % Upper boundary
    if xm(ii,i) > limits(ii,2)
      xm(ii,i) = limits(ii,2);
    end
  end

end

end
%% END