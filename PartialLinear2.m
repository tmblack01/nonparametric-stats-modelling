function m = PartialLinear2(pts,X,Y,kernfun)
%PARTIALLINEAR2 Implementation of the partial linear model where X is 
% 2-dimensional and Y is 1-dimensional
T = X(:,1);
U = X(:,2);
T_pts = pts(:,1);
U_pts = pts(:,2);

h = ROT(U);
S = kernfun((T - T')/h);
den = sum(kernfun((T - T')/h),2);
for i = 1:size(S,1)
    S(i,:) = S(i,:)/den(i);
end

U2 = (eye(size(S,1)) - S) * U;
Y2 = (eye(size(S,1)) - S) * Y;

beta = inv(U2'*U2)*U2'*Y2;

S_pts = kernfun((T_pts - T')/h);
den = sum(kernfun((T_pts - T')/h),2);
for i = 1:size(S_pts,1)
    S_pts(i,:) = S_pts(i,:)/den(i);
end

g = S_pts*(Y-U*beta');

m = U_pts*beta' + g;

end