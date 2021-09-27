function [gi] =Gi(X,Vi,gama_i)
temp = (X - Vi) ;
powere = temp * temp.' ;
gi = exp(-1 * gama_i * powere);
end