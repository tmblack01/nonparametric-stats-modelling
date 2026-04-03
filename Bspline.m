function B = Bspline(x, a, b, Knots, d)
% x: values of covariate where the B-spline functions are evaluated
% a: lower bound of the support of x
% b: upper bound of the support of x
% Knots: [a, interior knots, b], an L+2 array. E.g., a=0, b=1, interior
% knots = 0.05:0.05:0.95, then Knots = 0:0.05:1
% d: degree of the spline

L = length(Knots) - 2;
Knots = reshape(Knots, 1, L+2);
allknots = [repmat(a,1,d),Knots,repmat(b,1,d)]; 

B = fnval(spmak(allknots,eye(L+d+1)),x)';

end