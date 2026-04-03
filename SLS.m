function s = SLS(v,X,Y,kernfun)
%SLS implementation of the semiparametric least square (SLS) criterion for
%the single index model
   beta = v(1:length(v)-1);
   h = v(length(v));
   n = size(X,1);
   s = 0;
   for i = 1:n
       s = s + 1/n * (Y(i) - LOOLC(X,Y,beta,h,kernfun,i))^2;
   end

end

