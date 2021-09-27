function [output] = fun(X,Ystar,L,x00,mode)
sz = size(x00);
len=sz(1,1) / 4 ;
%% calculate V and cigma
V = rand(len,2);
cigma = rand(len,2);
for i=1:1:len 
    V(i,1) =x00(i,1);
end
for i=len+1:1:2*len 
    V(i-len,2) =x00(i,1);
end
for i=(2*len + 1):1:3*len 
    cigma(i - 2*len,1) =x00(i,1);
end
for i=(3*len+1):1:4*len 
     cigma(i - 3*len ,2) =x00(i,1);
end
%% calculate G
G = rand(L,len);
for l=1:1:L
  for i=1:1:len 
      t = (((X(l,1) - V(i,1)) / cigma(i,1)) ^ 2 )+(((X(l,2) - V(i,2)) / cigma(i,2))^ 2 );
      power= exp(-1 * t);
      if (power == 0 )
       G(l,i) = rand(1);
      else 
          G(l,i)  = power;
      end
  end
end        
%% calculate W and Y
GTranspose=G.';
newG=GTranspose*G;
newYstar=GTranspose*Ystar;
yy = size(newG);
for i=1:1:yy(1,1)
    for j = 1:1:yy(1,2)
        if( isnan(newG(i,j)) == 1 )
            newG(i,j) = -255;
        end
    end
end
W=pinv(newG) * newYstar;                
Y=G*W;

%% calcualte output
if mode == 1
 % regresion
 Yh= Ystar - Y;
 output = 1/2 .* (Yh.' * Yh);
else
 % classification
A = Ystar.' ;
B = Y.' ;
for i = 1:1:L
    mm = max(B(:,i));
    index = find(B(:,i) == mm);
    B(index,i) =1;
end
indexA = zeros(L,1);
indexB = zeros(L,1);
for i=1:1:L
       indexB(i,1) = find(B(:,i) == 1);
end
for i=1:1:L
       indexA(i,1) = find(A(:,i) == 1);
end
   temp = sign(abs(indexB - indexA));
   output = sum(temp) / L ;
end
end
