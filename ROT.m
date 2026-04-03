function h = ROT(X)
    %ROT find bandwidth of KDE using rule of thumb estimator
    std_est = min([std(X),iqr(X)/1.34]);
    h = 1.06 * std_est * size(X,1) ^ (-1/5);
end

