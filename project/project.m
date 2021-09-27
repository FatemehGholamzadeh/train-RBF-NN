%% inputs are X and Ystar
% 1 is circle and 2 is elliptic
model = 2
m=10;
M = xlsread('regdata2000 (1).xlsx');
s = size(M);
L=s(1,1);
dim = s(1,2);
dim = dim -1 ;
X = M(:,1:(dim) );
Y1 = M(:,dim+1);
minimum = min(Y1);
maximum = max(Y1);
if minimum == -1 
    mode = 2;
elseif (fix(maximum) - maximum ) == 0 
    mode = 3;
else 
    mode = 1;
end  
if mode == 1
   %% regresion mode
   Ystar = Y1;
elseif mode ==2 
    %% 2 classification mode
    Y2 = -1 .* Y1 ;
    Ystar = [Y1  Y2];
elseif mode ==3
    %% multi classfication mode 
    Ystar = zeros(L,maximum);
    for i=1:1:L
        for j = 1:1:maximum
            if Y1(i,1) == j
            Ystar(i,j) = 1;
            else
            Ystar(i,j) = -1;
            end
        end
    end
end

%%   use ES 
maxVector = max(X);
maxElement = max(maxVector.');
N = X ./ maxElement ;
if (model==1 )
   n_x=(dim+1)*m ;
elseif (model ==2)
   n_x=(dim+2)*m; 
end
limits = repmat([0.01 0.99], n_x, 1);
mu = 15;
lambda = 15;
gen =20;
sel = ',';
rec_obj = 2;
rec_str = 4;
obj=0;
nf=1;
if (model==1 )
   n=(dim+1)*m ;
elseif (model ==2)
   n =(dim+2)*m; 
end
X = N;
[V,minError,temp] = es(mu,lambda,gen,sel,rec_obj,rec_str,obj,nf,n,X,Ystar,L,limits,mode)
%% plot final result
%  plot(X(:,1),X(:,2),'*g','MarkerSize',12,'LineWidth',1);
% hold on;
 %plot(V(:,1),V(:,2),'xr','MarkerSize',15,'LineWidth',3); % new centers of clusters
 if (mode ==1)
     counter = zeros(L,1);
     for i=1:1:L
         counter(i,1)=i;
     end
    plot(counter(:,1),Ystar(:,1),'.g','MarkerSize',20,'LineWidth',2);
    hold on
    plot(counter(:,1),temp(:,1),'.b','MarkerSize',20,'LineWidth',2);
    hold on
 else
for i=1:1:L
    if(temp(i,1) == 0)
       plot(X(i,1),X(i,2),'.r','MarkerSize',20,'LineWidth',2); 
       hold on
    else
        plot(X(i,1),X(i,2),'.b','MarkerSize',20,'LineWidth',2);
    end
 end
 end