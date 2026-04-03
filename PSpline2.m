function m = PSpline2(pts, a, b, X, Y, L, lambda, D2)
    %PSpline2: implement 2-dimensional tensor-product cubic penalised
    % spline regression.
    
    % Calculate tensor product on observations
    knot1 = [a(1), quantile(X(:,1),L(1)), b(1)];
    knot2 = [a(2), quantile(X(:,2),L(2)), b(2)];
    B1 = Bspline(X(:,1), a(1), b(1), knot1, 3);
    B2 = Bspline(X(:,2), a(2), b(2), knot2, 3);
    B = row_kron(B1,B2);
    
    % Calculate tensor-product on points
    B1_pts = Bspline(pts(:,1), a(1), b(1), knot1, 3);
    B2_pts = Bspline(pts(:,2), a(2), b(2), knot2, 3);
    B_pts = row_kron(B1_pts,B2_pts);
    
    m = B_pts*pinv(B'*B + lambda*D2)*B'*Y;
end