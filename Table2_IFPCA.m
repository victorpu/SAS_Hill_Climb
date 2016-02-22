rand_if = zeros(1,50);
sdiff_if = zeros(1,50);
set0 = 1:50;
c0 = [ones(1,30),2*ones(1,30),3*ones(1,30)];
Sigma1 = diag(linspace(1,5,500));
Sigma2 = Sigma1;
Sigma3 = Sigma1;
for b = 1:50
    U1 = orth(randn(500));
    Sigma1t = transpose(U1) * Sigma1 * U1;

    U2 = orth(randn(500));
    Sigma2t = transpose(U2) * Sigma2 * U2;

    U3 = orth(randn(500));
    Sigma3t = transpose(U3) * Sigma3 * U3;

    mu = zeros(1,500);
    Y1 = mvnrnd(mu, Sigma1t, 30);
    Y2 = mvnrnd(mu, Sigma2t, 30);
    Y3 = mvnrnd(mu, Sigma3t, 30);
    Y = [Y1; Y2; Y3];

    mu1 = [linspace(1.02,2,50), zeros(1, 450)];
    mu2 = [linspace(2.02,3,50), zeros(1, 450)];
    mu3 = [linspace(3.02,4,50), zeros(1, 450)];
    X1 = repmat(mu1, 30, 1);
    X2 = repmat(mu2, 30, 1);
    X3 = repmat(mu3, 30, 1);
    X = [X1; X2; X3];

    X0 = X + Y;
    Data = X0';
    [IFlabel, stats, L]  = ifpca(Data, 3);
    set = stats.ranking(1:L);
    sdiff_if(b) = length(setxor(set,set0));
    [~,RI,~,~]= RandIndex(c0,IFlabel'); 
    rand_if(b) = RI;    
end
sprintf('mean(rand_if) is %f', mean(rand_if))
sprintf('std(rand_if) is %f', std(rand_if))
sprintf('mean(sdiff_if) is %f', mean(sdiff_if))
sprintf('mean(sdiff_if) is %f', mean(sdiff_if))
