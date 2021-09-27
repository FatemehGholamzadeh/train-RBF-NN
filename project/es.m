function [V,minError,temp] = es(mu,lambda,gen,sel,rec_obj,rec_str,obj,nf,n,X,Ystar,L,limits,mode,model)
%% Initialization:
s = size(X);
dim = s(1,2)
if ((sel ~= ',') && (sel ~= '+'))
  error('not supported selection scheme')
end
e     = 1e-8;  
[x_0, sigma, alpha] = validate_sizes(mu, n, limits);
min_x = cell(1,gen);                  % allocate space in memory for min_x
min_f = cell(1,gen);                  % allocate space in memory for min_f
off   = cell(1,gen);                  % allocate space to store the offspring population
temp = zeros(L,1);
min_x{1} = x_0;                       % first point
value = zeros(nf,mu);                 % allocate space for function evaluation
output =zeros(nf,mu);
for i = 1:1:mu
    x00 = x_0(:,i) ;
   [output(:,i),temp]=func(X,Ystar,L,x00,mode,dim);
  value(:,i) = output(:,i);
end
min_f{1} = value;                     % first approximation
off{1}   = zeros(n,1);

j      = 1;                           % generations counter
jj     = 0;                           % stagnation counter
 eps    = abs(obj - value(1,:));       % initial error
 EPS    = ones(gen,1);                % allocate space in memory for minimum error every generation
 EPS(1) = min(eps);
%% Begin ES
while ((j < gen ) && (min(eps) > e))
  %% Recombine:
  [xr,sigmar,alphar] = recombination(rec_obj,rec_str,n,mu,lambda,min_x{j},sigma,alpha);
  off{j+1}           = xr;            % offspring population
   
  %% Mutation:
  [xm,sigmam,alpham] = mutation(n,lambda,xr,sigmar,alphar,limits);
  %[xm,sigmam,alpham] = mutation(n,lambda,xr,sigmar,alphar,limits);
  %% Evaluation:
  phie = zeros(nf,lambda);
  for i = 1:lambda
      x00 = xm(:,i) ;
      if (model == 1)
          [output(:,i),temp]= func(X,Ystar,L,x00,mode,dim);
      elseif (model == 2 ) 
          [output(:,i),temp]= func(X,Ystar,L,x00,mode,dim);
      end
    phie(:,i) = output(:,i);
  end
  epse = abs(obj - phie(1,:));
  %% Selection:
  [min_x{j+1}, sigma, alpha] = selection(sel, mu, lambda, epse, eps, xm, sigmam, min_x{j}, sigma, alpham, alpha);
  %% Store better results:
  value = zeros(nf,mu);               % allocate space for function evaluation
  for i = 1:mu
      x00 = min_x{j+1}(:,i) ;
      [output(:,i),temp]=func(X,Ystar,L,x00,mode,dim);
    value(:,i) = output(:,i);
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
  j = j+1
  if (jj == 30)
    fprintf('\n\n\tError remains constant for 30 consecutive generations\n\n');
    break
  end
end
minError = min(EPS);
minIndex = find(EPS(:,1) == minError);
 values= min_f{minIndex(1,1)};
 minf = find(values(1,:) == minError);
 vv = min_x{minIndex(1,1)};
 v = vv(:,minf);
 if (model == 1)
      len = n/(dim+1);
 elseif (model ==2 )
      len = n/(dim+2);
 end
V =zeros(len,dim);
 for i = 1:1:dim
     V(:,i) = v(((i-1)*len+1) : (i*len),1);
 end
end
%% END