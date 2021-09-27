function [min_x, min_sigma, min_alpha] = selection(scheme, mu, lambda, epse, eps, xm, sigmam, x0, sigma, alpham, alpha)
%%  [min_x, min_sigma, min_alpha] = selection(scheme, mu, lambda, epse, eps, xm, sigmam, x0, sigma, alpham, alpha)
%
%   This function carries out the selection process in the "estrategy evolution
%   procedure".
%
%   INPUT DATA:
%
%   - scheme: Selection scheme:
%                   ',' = (mu, lambda)-selection scheme
%                   '+' = (mu + lambda)-selection scheme
%   - mu:     Parent population size (positive integer number)
%   - lambda: Offspring population size (positive integer number)
%   - epse:   Error given by the offspring population (nf x lambda matrix) - 
%             (nf: length handle function 'f')
%   - eps:    Error given by the parents population (nf x lambda matrix) - 
%             (nf: length handle function 'f')
%   - xm:     Mutataed individuals (n x lambda matrix)
%   - sigmam: Mutated standard deviations (nsigma x lambda matrix)
%   - x0:     Parents individuals (n x mu matrix)
%   - sigma:  Parents standard deviations (nsigma x mu matrix)
%   - alpham: Mutated rotation angles (nsigma x lambda matrix)
%   - alpha:  Parents rotation angles (nsigma x mu matrix)
%
%   OUTPUT DATA:
%
%   - min_x:     Individuals selected (n x mu matrix)
%   - min_sigma: Covarinces of the individuals selected (1 x mu cell, each cell
%                contains an nxn symmetric matrix)
%   - min_alpha: Rotation angles of the individuals selected (1 x mu cell, each
%                cell contains an nxn symmetric matrix)
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
switch scheme
  %% (mu, lambda)-selection scheme
  case ','
    if (mu > lambda)
      error('The parents population size is greater than the offspring population size');
    end
    err = epse;
    [xmin, idx] = sort(err);
    min_x       = xm(:,idx(1:mu));
    min_sigma   = sigmam(idx(1:mu));
    min_alpha   = alpham(idx(1:mu));

  %% (mu + lambda)-selection scheme
  case '+'
    err         = [epse eps];    
    xaug        = [xm x0];
    sigmaaug    = [sigmam sigma];
    alphaaug    = [alpham alpha];
    [xmin, idx] = sort(err);
    min_x       = xaug(:,idx(1:mu));
    min_sigma   = sigmaaug(idx(1:mu));
    min_alpha   = alphaaug(idx(1:mu));

  %% no suported selection scheme
  otherwise
    error('not supported selection scheme');
end

end
%% END