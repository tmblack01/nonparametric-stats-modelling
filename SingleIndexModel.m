function m = SingleIndexModel(pts,X,Y)
%SINGLEINDEXMODEL Implementation of the single index model where X is 
% 2-dimensional and Y is 1-dimensional
    beta = inv(X'*X)*X'*Y;
    h = 0.5*(ROT(X(:,1)) + ROT(X(:,2)));

    % find beta and bandwidth that minimises SLS under constraints
    v0 = [beta', h];
    v0(2) = 1;
    A = [0,0,-1];
    b = 0;
    Aeq = [0,1,0];
    beq = 1;
    v = fmincon(@(v)SLS(v,X,Y,@normpdf),v0,A,b,Aeq,beq);
    
    beta = v(1:length(v)-1);
    h = v(length(v));
    eta = X*beta';
    eta_pts = pts*beta';
    m = LC(eta,Y,eta_pts,h,@normpdf);
end

