function [xip,sigmap,alphap] = recombination(type_obj,type_str,n,mu,lambda,xi,sigma,alpha)
%%  [xip,sigmap,alphap] = recombination(type_obj,type_str,n,mu,lambda,xi,sigma,alpha)
%   This function performs the recombination required in evolutionary strategies
%   to obtain the next offspring population.
%
%   INPUT DATA:
%
%   - type_obj: Type of recombination to use on object variables
%               (Pag. 74 in (BACK)):
%                 * 1 = No recombination
%                 * 2 = Discrete recombination
%                 * 3 = Panmictic discrete recombination
%                 * 4 = Intermediate recombination
%                 * 5 = Panmictic intermediate recombination
%                 * 6 = Generalized intermediate recombination
%                 * 7 = Panmictic generalized intermediate recombination
%   - type_str: Type of recombination to use on strategy parameters
%               (Pag. 74 in (BACK))
%   - n:        Length of 'xi'
%   - mu:       Parents population size
%   - lambda:   Offspring population size
%   - xi:       State at time 'i'
%   - sigma:    Covariances at time 'i'
%   - alpha:    Rotation angles at time 'i'
%
%   OUTPUT DATA:
%
%   - xip:      Recombined state
%   - sigmap:   Recombined covariances
%   - alphap:   Recombined rotation angles
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
%   Date: 20 - Sep - 2011

%% Beginning

%% Allocating space in memory for vector 'xip', and cells 'sigmap', 'alphap'
xip    = zeros(n,lambda);
sigmap = cell(1,lambda);
alphap = cell(1,lambda);

%% Recombination of object variables
switch type_obj
  %% No recombination
  case 1
    for i = 1:lambda
      idx       = randsample(1:mu,1);
      xip(:,i)  = xi(:,idx);
    end

  %% Discrete recombination
  case 2
    for i = 1:lambda
      tmp       = randsample(1:mu,2);
      for j = 1:n
        idx       = randsample(tmp,1);
        xip(j,i)  = xi(j,idx);
      end
    end

  %% Panmictic discrete recombination
  case 3
    i = 1;
    while (i <= lambda)
      fixed = randsample(1:mu,1);
      for k = 1:mu
        for j = 1:n
          tmp       = randsample(1:mu,1);
          idx       = randsample([fixed tmp],1);
          xip(j,i)  = xi(j,idx);
        end
      end
      i = i+1;
    end

  %% Intermediate recombination
  case 4
    for i = 1:lambda
      tmp       = randsample(1:mu,2);
      xip(:,i)  = xi(:,tmp(1))  + (xi(:,tmp(2))  - xi(:,tmp(1)))/2;
    end

  %% Panmictic intermediate recombination
  case 5
    i = 1;
    while (i <= lambda)
      fixed = randsample(1:mu,1);
      for k = 1:mu
        for j = 1:n
          tmp       = randsample(1:mu,1);
          xip(j,i)  = xi(j,fixed)  + (xi(j,tmp)  - xi(j,fixed))/2;
        end
      end
      i = i+1;
    end

  %% Generalized intermediate recombination
  case 6
    chi = rand;
    for i = 1:lambda
      tmp       = randsample(1:mu,2);
      xip(:,i)  = xi(:,tmp(1))  + chi*(xi(:,tmp(2))  - xi(:,tmp(1)));
    end

  %% Panmictic generalized intermediate recombination
  case 7
    i = 1;
    while (i <= lambda)
      fixed = randsample(1:mu,1);
      for k = 1:mu
        for j = 1:n
          chi       = rand;
          tmp       = randsample(1:mu,1);
          xip(j,i)  = xi(j,fixed)  + chi*(xi(j,tmp)  - xi(j,fixed));
        end
      end
      i = i+1;
    end

  %% No supported recombination type
  otherwise
    error('Not supported recombination type');
end

%% Recombination of strategy parameters
switch type_str
  %% No recombination
  case 1
    for i = 1:lambda
      idx       = randsample(1:mu,1);
      sigmap{i} = sigma{idx};
      alphap{i} = alpha{idx};
    end

  %% Discrete recombination
  case 2
    for i = 1:lambda
      tmp       = randsample(1:mu,2);
      for j = 1:n
        for jj = j:n
          idx             = randsample(tmp,1);
          sigmap{i}(j,jj) = sigma{idx}(j,jj);
          sigmap{i}(jj,j) = sigma{idx}(j,jj);
          alphap{i}(j,jj) = alpha{idx}(j,jj);
        end
      end
    end

  %% Panmictic discrete recombination
  case 3
    i = 1;
    while (i <= lambda)
      fixed = randsample(1:mu,1);
      for k = 1:mu
        for j = 1:n
          for jj = j:n
            tmp             = randsample(1:mu,1);
            idx             = randsample([fixed tmp],1);
            sigmap{i}(j,jj) = sigma{idx}(j,jj);
            sigmap{i}(jj,j) = sigma{idx}(j,jj);
            alphap{i}(j,jj) = alpha{idx}(j,jj);
          end
        end
      end
      i = i+1;
    end

  %% Intermediate recombination
  case 4
    for i = 1:lambda
      tmp       = randsample(1:mu,2);
      sigmap{i} = sigma{tmp(1)} + (sigma{tmp(2)} - sigma{tmp(1)})/2;
      alphap{i} = alpha{tmp(1)} + (alpha{tmp(2)} - alpha{tmp(1)})/2;

      % validate rotation angles (they must be between [-pi, pi]):
      [p,m]          = find(abs(alphap{i}) > pi);
      alphap{i}(p,m) = alphap{i}(p,m) - 2*pi*(alphap{i}(p,m)/abs(alphap{i}(p,m)));
      
      % validate standard deviations (they must be greater than zero):
      [p,m]          = find(sigmap{i} <= 0);
      sigmap{i}(p,m) = 0.1;
    end

  %% Panmictic intermediate recombination
  case 5
    i = 1;
    while (i <= lambda)
      fixed = randsample(1:mu,1);
      for k = 1:mu
        for j = 1:n
          for jj = j:n
            tmp       = randsample(1:mu,1);
            sigmap{i}(j,jj) = sigma{fixed}(j,jj) + (sigma{tmp}(j,jj) - sigma{fixed}(j,jj))/2;
            sigmap{i}(jj,j) = sigma{fixed}(j,jj) + (sigma{tmp}(j,jj) - sigma{fixed}(j,jj))/2;
            alphap{i}(j,jj) = alpha{fixed}(j,jj) + (alpha{tmp}(j,jj) - alpha{fixed}(j,jj))/2;

            % validate rotation angles (they must be between [-pi, pi]):
            [p,m]          = find(abs(alphap{i}) > pi);
            alphap{i}(p,m) = alphap{i}(p,m) - 2*pi*(alphap{i}(p,m)/abs(alphap{i}(p,m)));
        
            % validate standard deviations (they must be greater than zero):
            [p,m]          = find(sigmap{i} <= 0);
            sigmap{i}(p,m) = 0.1;
          end
        end
      end
      i = i+1;
    end

  %% Generalized intermediate recombination
  case 6
    chi = rand;
    for i = 1:lambda
      tmp       = randsample(1:mu,2);
      sigmap{i} = sigma{tmp(1)} + chi*(sigma{tmp(2)} - sigma{tmp(1)})/2;
      alphap{i} = alpha{tmp(1)} + chi*(alpha{tmp(2)} - alpha{tmp(1)})/2;

      % validate rotation angles (they must be between [-pi, pi]):
      [p,m]          = find(abs(alphap{i}) > pi);
      alphap{i}(p,m) = alphap{i}(p,m) - 2*pi*(alphap{i}(p,m)/abs(alphap{i}(p,m)));
      
      % validate standard deviations (they must be greater than zero):
      [p,m]          = find(sigmap{i} <= 0);
      sigmap{i}(p,m) = 0.1;
    end

  %% Panmictic generalized intermediate recombination
  case 7
    i = 1;
    while (i <= lambda)
      fixed = randsample(1:mu,1);
      for k = 1:mu
        for j = 1:n
          for jj = j:n
            chi       = rand;
            tmp       = randsample(1:mu,1);
            sigmap{i}(j,jj) = sigma{fixed}(j,jj) + chi*(sigma{tmp}(j,jj) - sigma{fixed}(j,jj))/2;
            sigmap{i}(jj,j) = sigma{fixed}(j,jj) + chi*(sigma{tmp}(j,jj) - sigma{fixed}(j,jj))/2;
            alphap{i} = alpha{fixed} + chi*(alpha{tmp} - alpha{fixed})/2;

            % validate rotation angles (they must be between [-pi, pi]):
            [p,m]          = find(abs(alphap{i}) > pi);
            alphap{i}(p,m) = alphap{i}(p,m) - 2*pi*(alphap{i}(p,m)/abs(alphap{i}(p,m)));
        
            % validate standard deviations (they must be greater than zero):
            [p,m]          = find(sigmap{i} <= 0);
            sigmap{i}(p,m) = 0.1;
          end
        end
      end
      i = i+1;
    end

  %% No supported recombination type
  otherwise
    error('Not supported recombination type');
end

end
%% END