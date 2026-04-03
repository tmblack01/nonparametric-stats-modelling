function g = LOOLC(X,Y,beta,h,kernfun,i)
%LOOLC implementation of the leave one out local costant estimator which is
% used for estimating the parameters of the single index model
n = size(X,1);
    num = 0;
    den = 0;
    for j = 1:n
        if i ~= j
            num = num + 1/h * kernfun((X(i,:) * beta' - X(j,:) * beta')/h) * Y(i);
            den = den + 1/h * kernfun((X(i,:) * beta' - X(j,:) * beta')/h);
        end
    end
    g = num/den;
end

