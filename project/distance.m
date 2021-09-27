function[d] = distance(x,v,dim)
 d =0 ;
  for i =1:1:dim
     d = d + (x(1,i) -v(1,i) ) .^2 ;
  end
  d = sqrt(d);
end