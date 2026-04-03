function f = LC(X,Y,pts,h,kerfun)
    %LC implementation of local constant regression
    f = sum(kerfun((pts - X')/h).*Y',2)./sum(kerfun((pts - X')/h),2);
end

