%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of Non-Parametric Statistical Models
%
% Author: Thomas Black
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set random seed to repeat results
rng(90111)

% initialise values
J = 200;
n = 100;
mu = [0 0];
Sigma = [1 0.2;0.2 1];
mse1 = zeros([J,3]);
mse2 = zeros([J,3]);

for j = 1:J
    % generate X from a bivariate normal distribution
    X = mvnrnd(mu,Sigma,n);
    a = [min(X(:,1)),min(X(:,2))];
    b = [max(X(:,1)),max(X(:,2))];
    
    % DGP1:
    beta_0 = [1/sqrt(3), sqrt(2/3)];
    m1 = 0.1.*(X*beta_0').^2 .* exp(X*beta_0');
    Y1 = m1 + normrnd(0,1,[n,1]);
    
    % DGP2:
    m2 = cos(X(:,1)) + X(:,2);
    Y2 = m2 + normrnd(0,1,[n,1]);
    
    % Perform 2-dimensional tensor-product spline regression
    L = [min(floor(sqrt(n)/1.5),35),min(floor(sqrt(n)/1.5),35)];

    % get a grid of values for X1 and X2
    knots1 = [a(1),quantile(X(:,1),L(1)),b(1)];
    knots2 = [a(2),quantile(X(:,2),L(2)),b(2)];
    dx1 = (b(1) - a(1))/99;
    dx2 = (b(2) - a(2))/99;
    xseq1 = a(1):dx1:b(1);
    xseq2 = a(2):dx2:b(2);
    [xgrid1,xgrid2] = meshgrid(xseq1,xseq2);
    xgrid = [xgrid1(:),xgrid2(:)];
    
    % Calculate D2 using Riemann sum
    allknots1 = [repmat(a(1),1,3),knots1,repmat(b(1),1,3)];
    allknots2 = [repmat(a(2),1,3),knots2,repmat(b(2),1,3)];
    Derv2_B1 = fnval(fnder(spmak(allknots1,eye(L(1)+4)),2),xgrid(:,1))';
    Derv2_B2 = fnval(fnder(spmak(allknots2,eye(L(2)+4)),2),xgrid(:,2))';
    Derv2_B = row_kron(Derv2_B1,Derv2_B2);
    D2 = Derv2_B'*Derv2_B*dx1*dx2;

    lambda1 = GCV(X,Y1,L,a,b,D2);
    lambda2 = GCV(X,Y2,L,a,b,D2);
    m1_hat1 = PSpline2(X,a,b,X,Y1,L,lambda1,D2);
    m2_hat1 = PSpline2(X,a,b,X,Y2,L,lambda2,D2);
    
    % Perform Single Index Model Regression 
    m1_hat2 = SingleIndexModel(X,X,Y1);
    m2_hat2 = SingleIndexModel(X,X,Y2);
    
    % Perform Partial Linear Regression 
    m1_hat3 = PartialLinear2(X,X,Y1,@normpdf);
    m2_hat3 = PartialLinear2(X,X,Y2,@normpdf);
    
    % Calculate MSE for each model fitted on DGP-1
    mse1(j,1) = (1/n)*sum((m1_hat1 - m1).^2);
    mse1(j,2) = (1/n)*sum((m1_hat2 - m1).^2);
    mse1(j,3) = (1/n)*sum((m1_hat3 - m1).^2);
    
    % Calculate MSE for each model fitted on DGP-2
    mse2(j,1) = (1/n)*sum((m2_hat1 - m2).^2);
    mse2(j,2) = (1/n)*sum((m2_hat2 - m2).^2);
    mse2(j,3) = (1/n)*sum((m2_hat3 - m2).^2);
end

%% Plot Boxplot for DGP 1

figure
boxplot(mse1,'Labels',["penalised splines", "single index","partial linear"],DataLim=[-0.1,6])
xlabel("model")
ylabel('mean square error (mse)')
title("Mean Square Errors of Statistical Models (DGP-1)")

%% Plot Boxplot for DGP 2

figure
boxplot(mse2,'Labels',["penalised splines", "single index","partial linear"])
xlabel("model")
ylabel('mean square error (mse)')
title("Mean Square Errors of Statistical Models (DGP-2)")

%% Plot Spline Regression Model fitted to DGP-1

% points for plotting
pts1 = min(X(:,1)):0.05:max(X(:,1));
pts2 = min(X(:,2)):0.05:max(X(:,2));
[pts1,pts2] = meshgrid(pts1,pts2);
pts = [pts1(:),pts2(:)];

m1_hat1_plot = PSpline2(pts,a,b,X,Y1,L,lambda1,D2);
m1_hat1_plot = reshape(m1_hat1_plot,[length(unique(pts(:,2))),length(unique(pts(:,1)))]);
figure
surf(unique(pts(:,1)),unique(pts(:,2)),m1_hat1_plot)
hold on
scatter3(X(:,1),X(:,2),Y1,'filled',MarkerFaceColor="red")

%% Plot Single Index Model fitted to DGP-1

m1_hat2_plot = SingleIndexModel(pts,X,Y1);
m1_hat2_plot = reshape(m1_hat2_plot,[length(unique(pts(:,2))),length(unique(pts(:,1)))]);
figure
surf(unique(pts(:,1)),unique(pts(:,2)),m1_hat2_plot)
hold on
scatter3(X(:,1),X(:,2),Y1,'filled',MarkerFaceColor="red")

%% Plot Partial Linear Model fitted to DGP-1

m1_hat3_plot = PartialLinear2(pts,X,Y1,@normpdf);
m1_hat3_plot = reshape(m1_hat3_plot,[length(unique(pts(:,2))),length(unique(pts(:,1)))]);
figure
surf(unique(pts(:,1)),unique(pts(:,2)),m1_hat3_plot)
hold on
scatter3(X(:,1),X(:,2),Y1,'filled',MarkerFaceColor="red")

%% Plot Spline Model fitted to DGP-2

m2_hat1_plot = PSpline2(pts,a,b,X,Y2,L,lambda2,D2);
m2_hat1_plot = reshape(m2_hat1_plot,[length(unique(pts(:,2))),length(unique(pts(:,1)))]);
figure
surf(unique(pts(:,1)),unique(pts(:,2)),m2_hat1_plot)
hold on
scatter3(X(:,1),X(:,2),Y2,'filled',MarkerFaceColor="red")

%% Plot Single Index Model fitted to DGP-2

m2_hat2_plot = SingleIndexModel(pts,X,Y2);
m2_hat2_plot = reshape(m2_hat2_plot,[length(unique(pts(:,2))),length(unique(pts(:,1)))]);
figure
surf(unique(pts(:,1)),unique(pts(:,2)),m2_hat2_plot)
hold on
scatter3(X(:,1),X(:,2),Y2,'filled',MarkerFaceColor="red")

%% Plot Partial Linear Model fitted to DGP-2

m2_hat3_plot = PartialLinear2(pts,X,Y2,@normpdf);
m2_hat3_plot = reshape(m2_hat3_plot,[length(unique(pts(:,2))),length(unique(pts(:,1)))]);
figure
surf(unique(pts(:,1)),unique(pts(:,2)),m2_hat3_plot)
hold on
scatter3(X(:,1),X(:,2),Y2,'filled',MarkerFaceColor="red")