%%  EXAMPLES
%
%   Bibliography:
%
%   - BACK, Thomas. "Evolutionary algorithms in theory and practice".
%     Oxford University Press. New York. 1996.
%   - HAUPT, Randy L.; HAUPT, Sue Ellen. "Practical Genetic Algorithms".
%     2nd editon. John Wiley & Sons, Inc. Pennsylvania, USA. 2004.
%   - OLDENHUIS, Rody. Many test functions for global optimizers:
%     "http://www.mathworks.com/matlabcentral/fileexchange/23147-many-
%      testfunctions-for-global-optimizers"
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
clear, clc, close all

fun = 1;                    % function to minimize

switch fun
  case 1
    %%  Rosenbrock's function
    %   Minimum: f(1,1) = 0
    f      = @(x,u) (1-x(1,:)).^2 + 100*(x(2,:)-x(1,:).^2).^2;
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 2
    %%  Himmelblau's function
    %   Minimum:  * f( 3,         2       ) = 0
    %             * f(-2.805118,  3.131312) = 0
    %             * f(-3.779310, -3.283186) = 0
    %             * f( 3.584428, -1.848126) = 0
    f      = @(x,u) (x(1,:).^2 + x(2,:) - 11).^2 + (x(1,:) + x(2,:).^2 - 7).^2;
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)
       
  case 3
    %%  Rastrigin function
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) 20 + (x(1,:).^2 + x(2,:).^2) - 10*(cos(2*pi*x(1,:)) + cos(2*pi*x(2,:)));
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 4
    %%  Sphere function
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) x(1,:).^2 + x(2,:).^2;
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)
    
  case 5
    %%  Ackley's function
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) -20*exp(-0.2*sqrt(0.5*(x(1,:).^2 + x(2,:).^2))) - exp(0.5*(cos(2*pi*x(1,:)) + cos(2*pi*x(2,:)))) + 20 + exp(1);
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)
  
  case 6
    %%  F1 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(0) = 1
    f      = @(x,u) abs(x(1,:)) + cos(x(1,:));
    n_x    = 1;                           % 'n_x' states
    limits = repmat([-20 20], n_x, 1);    % Boundaries
    obj    = 1;                           % objective value (f(x_min) = obj)
    
  case 7
    %%  F2 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(0) = 0
    f      = @(x,u) abs(x(1,:)) + sin(x(1,:));
    n_x    = 1;                           % 'n_x' states
    limits = repmat([-20 20], n_x, 1);    % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 8
    %%  F5 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(zeros(1,n_x)) = -20
    n_x    = 2;                           % 'n_x' states
    f      = @(x,u) sum(abs(x(1:n_x,:)) - 10*cos(sqrt(abs(10*x(1:n_x,:)))));
    limits = repmat([-10 10], n_x, 1);    % Boundaries
    obj    = -20;                         % objective value (f(x_min) = obj)

  case 9
    %%  F6 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(9.6204) = -100.22
    f      = @(x,u) (x(1,:).^2 + x(1,:)).*cos(x(1,:));
    n_x    = 1;                           % 'n_x' states
    limits = repmat([-10 10], n_x, 1);    % Boundaries
    obj    = -100.22;                     % objective value (f(x_min) = obj)
    
  case 10
    %%  F7 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(9.039,8.668) = -18.5547
    f = @(x,u) x(1,:).*sin(4*x(1,:)) + 1.1*x(2,:).*sin(2*x(2,:));
    n_x    = 2;                           % 'n_x' states
    limits = repmat([0 10], n_x, 1);      % Boundaries
    obj    = -18.5547;                    % objective value (f(x_min) = obj)

  case 11
    %%  F8 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(9.039,8.668) = -18.5547
    f = @(x,u) x(2,:).*sin(4*x(1,:)) + 1.1*x(1,:).*sin(2*x(2,:));
    n_x    = 2;                           % 'n_x' states
    limits = repmat([0 10], n_x, 1);      % Boundaries
    obj    = -18.5547;                    % objective value (f(x_min) = obj)
  
  case 12
    %%  F11 (See Appendix I (Haupt, 2004)).
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) 1 + ((x(1,:).^2 + x(2,:).^2)/4000) - cos(x(1,:)).*cos(x(2,:));
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-10 10], n_x, 1);    % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 13
    %%  F12 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(zeros(1,n_x)) = 0
    n_x    = 2;                           % 'n_x' states
    f      = @(x,u) 0.5 + ((sin(sqrt(sum(x(1:n_x,:).^2)))).^2 - 0.5)./(1 + 0.1*(sum(x(1:n_x,:).^2)));
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 14
    %%  F13 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(0,0) = 0
    f      = @(x,u) (x(1,:).^2 + x(2,:).^2).^(0.25).*sin(30*((x(1,:) + 0.5).^2 + x(2,:).^2).^0.1) + abs(x(1,:)) + abs(x(2,:));
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-10 10], n_x, 1);    % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 15
    %%  F15 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(0,0.7688) = -345.36
    f      = @(x,u) -exp(-0.2*sqrt(x(1,:).^2 + x(2,:).^2) + 3*(cos(2*x(1,:)) + sin(2*x(2,:))));
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = -345.36;                     % objective value (f(x_min) = obj)

  case 16
    %%  F16 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(-20,3.15) = -23.59
    f      = @(x,u) -x(1,:).*sin(sqrt(abs(x(1,:) - (x(2,:)+9)))) - (x(2,:)+9).*sin(sqrt(abs(x(2,:) + 0.5*x(1,:) + 9)));
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-20 20], n_x, 1);    % Boundaries
    obj    = -23.59;                      % objective value (f(x_min) = obj)

  case 17
    %%  Other function (no reference)
    %   Minimum: f(0,0) = -20.09
    f      = @(x,u) -exp(-0.2*sqrt(x(1,:).^2 + x(2,:).^2) + 3*(cosd(2*x(1,:)) + sind(2*x(2,:))));
    n_x    = 2;                         % 'n_x' states
    limits = repmat([-5 5], n_x, 1);    % Boundaries
    obj    = -20.09;                    % objective value (f(x_min) = obj)

  case 18
    %%  Other function (No reference)
    %   Minimum: f(zeros(1,n_x)) = 0
    n_x    = 2;                           % 'n_x' states
    f      = @(x,u) 0.5 + ((sqrt(sum(x(1:n_x,:).^2))).^2 - 0.5)./(1 + 0.1*(sum(x(1:n_x,:).^2)));
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 19
    %%  Beale's function
    %   Minimum: f(3,0.5) = 0
    f      = @(x,u) (1.5 - x(1,:) + x(1,:).*x(2,:)).^2 + (2.25 - x(1,:) + x(1,:).*x(2,:).^2).^2 + (2.625 - x(1,:) + x(1,:).*x(2,:).^3).^2;
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-4.5 4.5], n_x, 1);  % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)
   
  case 20
    %%  Rosenbrock's function (Generalization to n_x dimensions)
    %   Minimum: f(1,1,...,1) = 0
    n_x    = 4;                           % 'n_x' states (variable)
    f      = @(x,u) sum(100*(x(2:2:n_x, :) - x(1:2:n_x-1, :).^2).^2 + (1 - x(1:2:n_x-1, :)).^2, 1);
    limits = repmat([-100 100], n_x, 1);  % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)
    
  case 21
    %%  Bird's function
    %   Minimum:  * f( 4.701055751981055,  3.152946019601391) = -106.764537
    %             * f(-1.582142172055011, -3.130246799635430) = -106.764537
    f      = @(x,u) sin(x(1,:)).*exp((1-cos(x(2,:))).^2) + cos(x(2,:)).*exp((1-sin(x(1,:))).^2) + (x(1,:)-x(2,:)).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-2*pi 2*pi], n_x, 1);  % Boundaries
    obj    = -106.764537;                   % objective value (f(x_min) = obj)

  case 22
    %%  Booth's function
    %   Minimum: f(1,3) = 0
    f      = @(x,u) (x(1,:) + 2*x(2,:) - 7).^2 + (2*x(1,:) + x(2,:) - 5).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 23
    %%  Bukin4 function
    %   Minimum: f(-10,0) = 0
    f      = @(x,u) 100*x(2,:).^2 + 0.01*abs(x(1,:) + 10);
    n_x    = 2;                             % 'n_x' states
    limits = [...
                  -15 -5
                  -3   3
             ];                             % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 24
    %%  Bukin6 function
    %   Minimum: f(-10,1) = 0
    f      = @(x,u) 100*sqrt(abs(x(2,:) - 0.01*x(1,:).^2)) + 0.01*abs(x(1,:) + 10);
    n_x    = 2;                             % 'n_x' states
    limits = [...
                  -15 -5
                  -3   3
             ];                             % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 25
    %%  Carrom Table function
    %   Minimum:  * f( 9.646157266348881,  9.646134286497169) = 24.1568155
    %             * f(-9.646157266348881,  9.646134286497169) = 24.1568155
    %             * f( 9.646157266348881, -9.646134286497169) = 24.1568155
    %             * f(-9.646157266348881, -9.646134286497169) = 24.1568155
    f      = @(x,u) -((cos(x(1,:)).*cos(x(2,:)).*exp(abs(1 - sqrt(x(1,:).^2 + x(2,:).^2)/pi))).^2)/30;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = 24.1568155;                    % objective value (f(x_min) = obj)

  case 26
    %%  Chichinadze function
    %   Minimum: f(5.90133, 0.5) = -43.3159
    f      = @(x,u) x(1,:).^2 - 12*x(1,:) + 11 + 10*cos(pi*x(1,:)/2) + 8*sin(5*pi*x(1,:)/2) - 1/sqrt(5)*exp(-((x(2,:) - 0.5).^2)/2);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-30 30], n_x, 1);      % Boundaries
    obj    = -43.3159;                      % objective value (f(x_min) = obj)
    
  case 27
    %%  Cross function
    %   Minimum: f(NaN, NaN) = 0  ---> Too complicated to determine the real solution
    f      = @(x,u) (abs(sin(x(1,:)).*sin(x(2,:)).*exp(abs(100 - sqrt(x(1,:).^2 + x(2,:).^2)/pi))) + 1).^(-0.1);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)
    
  case 28
    %%  Cross-in-tray function
    %   Minimum:  * f( 1.34940668535334, -1.34940668535334) = -2.062611870822739
    %             * f( 1.34940668535334,  1.34940668535334) = -2.062611870822739
    %             * f(-1.34940668535334,  1.34940668535334) = -2.062611870822739
    %             * f(-1.34940668535334, -1.34940668535334) = -2.062611870822739
    f      = @(x,u) -0.0001*(abs(sin(x(1,:)).*sin(x(2,:)).*exp(abs(100 - sqrt(x(1,:).^2 + x(2,:).^2)/pi))) + 1).^(0.1);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = -2.062611870822739;            % objective value (f(x_min) = obj)
    
  case 29
    %%  Cross-leg-table function
    %   Minimum:  f(NaN,NaN) = -1  ---> Too complicated to determine the real solution
    f      = @(x,u) -(abs(sin(x(1,:)).*sin(x(2,:)).*exp(abs(100 - sqrt(x(1,:).^2 + x(2,:).^2)/pi))) + 1).^(-0.1);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = -1;                            % objective value (f(x_min) = obj)

  case 30
    %%  Crowned cross function
    %   Minimum:  f(NaN,NaN) = 0.0001 ---> Too complicated to determine the real solution
    f      = @(x,u) 0.0001*(abs(sin(x(1,:)).*sin(x(2,:)).*exp(abs(100 - sqrt(x(1,:).^2 + x(2,:).^2)/pi))) + 1).^(0.1);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = -1;                            % objective value (f(x_min) = obj)
    
  case 31
    %%  Cube function
    %   Minimum:  f(1,1,...,1) = 0
    n_x    = 4;                             % 'n_x' states
    f      = @(x,u) sum(100*(x(2:n_x,:) - x(1:(n_x-1),:).^3).^2 + (1 - x(1:(n_x-1),:)).^2,1);
    limits = repmat([-20 20], n_x, 1);      % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)
    
  case 32
    %%  Easom function
    %   Minimum:  f(pi,pi) = -1
    f      = @(x,u) -cos(x(1,:)).*cos(x(2,:)).*exp(-((x(1,:)-pi).^2 + (x(2,:)-pi).^2));
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = -1;                            % objective value (f(x_min) = obj)
    
  case 33
    %%  Eggholder function (Just 2D)
    %   Minimum:  f(512, 404.2319) = -959.6407
    f      = @(x,u) sum(-(x(2,:)+47).*sin(sqrt(abs(x(2,:)+x(1,:)/2+47)))-x(1,:).*sin(sqrt(abs(x(1,:)-(x(2,:)+47)))), 1);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-512 512], n_x, 1);    % Boundaries
    obj    = -959.6407;                     % objective value (f(x_min) = obj)
    
  case 34
    %%  Giunta function (Just 2D)
    %   Minimum:  f(0.45834282, 0.45834282) = 0.0602472184
    n_x    = 2;                             % 'n_x' states
    f      = @(x,u) 0.6 + sum(sin(16*x(:,:)/15 - 1) + sin(16*x(:,:)/15 - 1).^2 + sin(4*(16*x(:,:)/15 - 1))/50, 1);
    limits = repmat([-1 1], n_x, 1);        % Boundaries
    obj    = 0.0602472184;                  % objective value (f(x_min) = obj)

  case 35
    %%  Goldstein Price function
    %   Minimum:  f(0,-1) = 3
    f      = @(x,u) (1  + (x(1,:) + x(2,:) + 1).^2.*(19 - 14*x(1,:) +  3*x(1,:).^2 - 14*x(2,:) +  6*x(1,:).*x(2,:) +  3*x(2,:).^2)).*...
                    (30 + (2*x(1,:) - 3*x(2,:)).^2.*(18 - 32*x(1,:) + 12*x(1,:).^2 + 48*x(2,:) - 36*x(1,:).*x(2,:) + 27*x(2,:).^2));
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-2 2], n_x, 1);        % Boundaries
    obj    = 3;                             % objective value (f(x_min) = obj)

  case 36
    %%  Griewank function
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) (x(1,:).^2 + x(2,:).^2)/200 - cos(x(1,:)).*cos(x(2,:)/sqrt(2)) + 1;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 37
    %%  Helical Valley function
    %   Minimum:  f(1,0,0) = 0
    f      = @(x,u) 100*((x(3,:) - 10*atan2(x(2,:), x(1,:))/2/pi).^2 + (sqrt(x(1,:).^2 + x(2,:).^2) - 1).^2) + x(3,:).^2;
    n_x    = 3;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 38
    %%  Holder table function
    %   Minimum:  * f( 8.055023472141116, 9.664590028909654) = -19.2085026
    %             * f(-8.055023472141116, 9.664590028909654) = -19.2085026
    %             * f( 8.055023472141116,-9.664590028909654) = -19.2085026
    %             * f(-8.055023472141116,-9.664590028909654) = -19.2085026
    f      = @(x,u) -abs(sin(x(1,:)).*cos(x(2,:)).*exp(abs(1 - sqrt(x(1,:).^2 + x(2,:).^2)/pi)));
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = -19.2085026;                   % objective value (f(x_min) = obj)

  case 39
    %%  Levi13 function
    %   Minimum:  f(1,1) = 0
    f      = @(x,u) sin(3*pi*x(1,:)).^2 + (x(1,:)-1).^2.*(1 + sin(3*pi*x(2,:)).^2) + (x(2,:)-1).^2.*(1 + sin(2*pi*x(2,:)).^2);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 40
    %%  Matyas function
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) 0.26*(x(1,:).^2 + x(2,:).^2) - 0.48*x(1,:).*x(2,:);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 41
    %%  McCormick function
    %   Minimum:  f(-0.54719, -1.54719) = -1.9133
    f      = @(x,u) sin(x(1,:) + x(2,:)) + (x(1,:)-x(2,:)).^2 - 1.5*x(1,:) + 2.5*x(2,:) + 1;
    n_x    = 2;                             % 'n_x' states
    limits = [...
                -1.5   4
                -3     4
             ];                             % Boundaries
    obj    = -1.9133;                       % objective value (f(x_min) = obj)

  case 42
    %%  Schaffer1 function
    %   Minimum:  f(0, 0) = 0
    f      = @(x,u) 0.5  + (sin(x(1,:).^2 + x(2,:).^2).^2 - 0.5) ./ (1+0.001*(x(1,:).^2 + x(2,:).^2)).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 43
    %%  Schaffer2 function
    %   Minimum:  f(0, 0) = 0
    f      = @(x,u) 0.5  + (sin(x(1,:).^2 - x(2,:).^2).^2 - 0.5) ./ (1+0.001*(x(1,:).^2 + x(2,:).^2)).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 44
    %%  Schaffer3 function
    %   Minimum:  f(0, 1.25313) = 0.00156685
    f      = @(x,u) 0.5  + (sin(cos(abs(x(1,:).^2 - x(2,:).^2))).^2 - 0.5) ./ (1+0.001*(x(1,:).^2 + x(2,:).^2)).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0.00156685;                    % objective value (f(x_min) = obj)

  case 45
    %%  Schaffer4 function
    %   Minimum:  f(0, 1.25313) = 0.292579
    f      = @(x,u) 0.5  + (cos(sin(abs(x(1,:).^2 - x(2,:).^2))).^2 - 0.5) ./ (1+0.001*(x(1,:).^2 + x(2,:).^2)).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0.292579;                      % objective value (f(x_min) = obj)

  case 46
    %%  Pen-Holder function
    %   Minimum:  * f(-9.64617, 9.64617) = -0.96353
    %             * f( 9.64617,-9.64617) = -0.96353
    %             * f( 9.64617, 9.64617) = -0.96353
    %             * f(-9.64617,-9.64617) = -0.96353
    f      = @(x,u) -exp(-(abs(cos(x(1,:)).*cos(x(2,:)).*exp(abs(1 - sqrt(x(1,:).^2 + x(2,:).^2)/pi)))).^(-1));
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-11 11], n_x, 1);      % Boundaries
    obj    = -0.96353;                      % objective value (f(x_min) = obj)

  case 47
    %%  Powell function
    %   Minimum:  f(0,0,0,0) = 0
    f      = @(x,u) (x(1,:) + 10*x(2,:)).^2 + 5*(x(3,:) - x(4,:)).^2 + (x(2,:) - 2*x(3,:)).^4 + 10*(x(1,:) - x(4,:)).^4;
    n_x    = 4;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 48
    %%  Schweffel function
    %   Minimum:  f(420.9687, 420.9687) = -837.9658
    f      = @(x,u) -x(1,:).*sin(sqrt(abs(x(1,:)))) -x(2,:).*sin(sqrt(abs(x(2,:))));
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-500 500], n_x, 1);    % Boundaries
    obj    = -837.9658;                     % objective value (f(x_min) = obj)

  case 49
    %%  Sine Envelope Sine function
    %   Minimum:  f(0,0,...,0) = 0
    n_x    = 4;                             % 'n_x' states (variable)
    f      = @(x,u) sum((sin(sqrt(x(1:n_x-1, :).^2 + x(2:n_x, :).^2)).^2 - 0.5)./(1 + 0.001*(x(1:n_x-1, :).^2 + x(2:n_x, :).^2)).^2 + 0.5, 1);
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 50
    %%  Six Hump Camel Back function
    %   Minimum:  * f( 8.984201368301331e-2,-7.126564032704135e-001) = -1.031628453489877
    %             * f(-8.984201368301331e-2, 7.126564032704135e-001) = -1.031628453489877
    f      = @(x,u) (4 - 2.1*x(1,:).^2 + x(1,:).^4/3).*x(1,:).^2 + x(1,:).*x(2,:) + (4*x(2,:).^2 - 4).*x(2,:).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-5 5], n_x, 1);        % Boundaries
    obj    = -1.031628453489877;            % objective value (f(x_min) = obj)

  case 51
    %%  Styblinski-Tang function
    %   Minimum:  f(-2.903534, -2.903534, ..., -2.903534) = -39.16599*n_x
    n_x    = 2;                             % 'n_x' states (Variable)
    f      = @(x,u) sum(x(:,:).^4 - 16*x(:,:).^2 + 5*x(:,:), 1)/2;
    limits = repmat([-5 5], n_x, 1);        % Boundaries
    obj    = -39.16599 * n_x;               % objective value (f(x_min) = obj)

  case 52
    %%  Test tube Holder function
    %   Minimum:  f(pi/2,0) = -10.8723
    f      = @(x,u) -4*abs(sin(x(1,:)).*cos(x(2,:)).*exp(abs(cos((x(1,:).^2 + x(2,:).^2)/200))));
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = -10.872299901558;              % objective value (f(x_min) = obj)

  case 53
    %%  Three Hump Camel function
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) 2*x(1,:).^2 - 1.05*x(1,:).^4 + x(1,:).^6/6 + x(1,:).*x(2,:) + x(2,:).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-5 5], n_x, 1);        % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 54
    %%  Wood function
    %   Minimum:  f(1,1,1,1) = 0
    f      = @(x,u) 100*(x(1,:).^2 - x(2,:)).^2 + (x(1,:) - 1).^2 + (x(3,:) - 1).^2 + 90*(x(3,:).^2 - x(4,:)).^2 + ...
                    10.1*((x(2,:) - 1).^2 + (x(4,:) - 1).^2)  + 19.8*(x(2,:) - 1).*(x(4,:) - 1);
    n_x    = 4;                             % 'n_x' states
    limits = repmat([-50 50], n_x, 1);      % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 55
    %%  Zetti function
    %   Minimum:  f(-0.0299, 0) = -0.003791
    f      = @(x,u) (x(1,:).^2 + x(2,:).^2 - 2*x(1,:)).^2 + x(1,:)/4;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-5 5], n_x, 1);        % Boundaries
    obj    = -0.003791;                     % objective value (f(x_min) = obj)

  otherwise
    error('Not supported equation');
end

%% Setting initial parameters
nf      = 1;                 % length of the output vector 'f(x,y)'
mu      = 100;               % parent population size
lambda  = 100;               % offspring population size
gen     = 100;               % number of generations
sel     = '+';               % Selection scheme (Pag. 78 in (BACK))
rec_obj = 2;                 % Type of recombination to use on object
                             % variables (Pag. 74 in (BACK))
                             % See 'recombination.m'
rec_str = 4;                 % Type of recombination to use on strategy
                             % parameters (Pag. 74 in (BACK))

u       = 0;                 % external excitation

%% Plot function in 3D
if n_x == 2                  % In case of a 2D function, generate surface
  figure
  [X,Y] = meshgrid(linspace(limits(1,1),limits(1,2),100),linspace(limits(2,1),limits(2,2),100));
  Z = reshape(f([X(:) Y(:)]',[]), 100, 100);
  surf(X,Y,Z)
  shading interp
  title('function','FontSize',18);
end

%% Run "Evolutionary Strategy" (ES):
[min_x, min_f, off, EPS,idx] = evolution_strategy(f, mu, lambda, gen, sel, rec_obj, rec_str, u, obj, nf, n_x, limits);

%% Plot simulation results

%% In case of 1D function
if n_x == 1
  Res    = zeros(1,idx);

  for j = 1:idx
    Res(j) = min_x{j}(1);    % Pass data from cell to matrix
  end

  figure
  plot(Res)
  grid on
  xlabel('Generation','FontSize',16);
  ylabel('x value','FontSize',16);
  title('Evolutionary strategy result','FontSize',18);
  
  % Error shown in semilogarithmic plot
  figure
  semilogy(EPS)
  grid on
  xlabel('Generation','FontSize',16);
  ylabel('log(error)','FontSize',16);
  title('Minimum error at each generation','FontSize',18);

%% In case of 2D function
elseif n_x == 2
  Res    = cell(1,2);
  Res{1} = zeros(1,idx);
  Res{2} = zeros(1,idx);

  for j = 1:idx
    Res{1}(j) = min_x{j}(1);
    Res{2}(j) = min_x{j}(2);
  end

  labels = {'x value', 'y value'};

  for i = 1:2
    figure
    plot(Res{i})
    grid on
    xlabel('Generation','FontSize',16);
    ylabel(labels{i},'FontSize',16);
    title('Evolutionary strategy result','FontSize',18);
  end
  
  % Error shown in semilogarithmic plot
  figure
  semilogy(EPS)
  grid on
  xlabel('Generation','FontSize',16);
  ylabel('log(error)','FontSize',16);
  title('Minimum error at each generation','FontSize',18);

%% In case of functions that cannot be plotted, just the error is shown in
%% semilogarithmic plot
else
  for i = 1:idx
    if (size(min_x{i}) == zeros(1,2))
      idx = idx - 1;
    end
  end
  figure
  semilogy(EPS(1:idx));
  grid on
  xlabel('Generation','FontSize',16);
  ylabel('log(error)','FontSize',16);
  title('Evolutionary strategy result','FontSize',18);
end

%% END