function [x_0, sigma, alpha] = validate_sizes(mu, n, limits, x_0, sigma)
%%  [x_0, sigma, alpha] = validate_sizes(mu, n, limits, x_0, sigma)
%
%   This function just validates the sizes of the matries x_0 and sigma. Also,
%   computes the rotation angles matrix 'alpha'.
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

%% Beginning

%%  Validate 'x_0' size (if it exists), if 'x_0' doesn't exist, create the
%   vector 'x_0':
if ((nargin < 3) || (nargin > 5))
  error('Number of input arguments is incomplete (< 3) or there are more (> 5)');

elseif nargin == 3
  %% Generate random parent population:
  x_0 = zeros(n,mu);
  for i = 1:n
    x_0(i,:) = unifrnd(limits(i,1), limits(i,2), 1, mu);      % initialization
  end

  %% Generate covariance matrix
  sigma  = cell(1,mu);
  for i = 1:mu
    % generate random symmetric covariance matrices
    tmp      = rand(n,n);
    sigma{i} = tmp + tmp';
  end

elseif nargin == 4
  %% Validate size of 'x_0'
  if ~isequal(size(x_0), [n mu])
    error('The dimension of x_0 must be n x mu')
  end
  %% Generate covariance matrix
  sigma  = cell(1,mu);
  for i = 1:mu
    % generate random symmetric covariance matrices
    tmp      = rand(n,n);
    sigma{i} = tmp + tmp';
  end
else
  %% Validate size of 'x_0'
  if ~isequal(size(x_0), [n mu])
    error('The dimension of x_0 must be n x mu')
  end
  
  %% Validate size of 'sigma'
  if ~isequal(size(sigma), [1 mu])
    error('Sigma must be a "cell" of size 1xmu. Each cell must contain an nxn symmetric covariance matrix')
  end

end

%% Compute 'alpha' matrices from the covariance matrices (Eq. 2.13 from BACK)
alpha = cell(1,mu);
for i = 1:mu
  alpha{i} = zeros(n);
  for j = 1:n-1
    for k = j+1:n
      alpha{i}(j,k) = 0.5*atan2(2*sigma{i}(j,k),(sigma{i}(j,j)^2 - sigma{i}(k,k)^2));
    end
  end
end

end
%% END