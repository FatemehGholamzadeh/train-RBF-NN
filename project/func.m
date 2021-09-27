function [output,temp] = func(X,Ystar,L,x00,mode,dim)
sz = size(x00);
len = sz(1,1) ./ (dim+1) ;
%% calculate V and cigma
V = rand(len,dim);
    for i = 1:1:dim
        V(:,i) = x00(((i-1)*len+1) : (i*len),1);
    end
    gama = zeros(len,1);
for i = ((dim )*len +1) : 1 :(dim +1) * len
    gama((i - dim*len),1)=x00(i,1);
end
%% calculate G
G = rand(L,len);
for l=1:1:L
  for u=1:1:len 
        G(l,u) = Gi(X(l,:),V(u,:),gama(u,1));
  end
end  
%% calculate W and Y
GTranspose=G.';
invs = inv(GTranspose * G + 5*eye(length(GTranspose * G)));
W = invs * GTranspose * Ystar ;
Y=G*W;
%% calcualte output
if mode == 1
 %% regresion
 Yh= Ystar - Y;
 output = 1/2 .* (Yh.' * Yh);
 temp = Y;
else
%% classification
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
   output = sum(temp) / L 
end
end