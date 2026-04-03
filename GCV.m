function lambda = GCV(X,Y,L,a,b,D2)
%GCV perform generalised cross-validation to select the smoothing parameter
% lambda for the 2-dimensional cubic penalised spline regression
    
    % create the omega matrix
    knots1 = [a(1),quantile(X(:,1),L(1)),b(1)];
    knots2 = [a(2),quantile(X(:,2),L(2)),b(2)];
    B1 = Bspline(X(:,1),a(1),b(1),knots1,3);
    B2 = Bspline(X(:,2),a(2),b(2),knots2,3);
    B = row_kron(B1,B2);
    
    % loop through each candidate lambda value and find the one that
    % minimises the CV criterion
    cv_min = Inf;
    n = size(X,1);
    for lambda_can = 0:0.1:8
        omega = B*pinv(B'*B + lambda_can*D2)*B';
        m = PSpline2(X,a,b,X,Y,L,lambda_can,D2);
        cv = n/((n - trace(omega))^2)*sum((Y-m).^2);

        if cv < cv_min
            cv_min = cv;
            lambda = lambda_can;
        end

    end
                

end


