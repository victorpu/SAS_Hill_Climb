c0 = [ones(1,20),2*ones(1,20),3*ones(1,20)];
P = [100,200,500,1000];
u = [0.6,0.7,0.8,0.9,1.0];
fid = fopen('simulation1.txt', 'w');
for i = 1:4
    p = P(i);
        fprintf(fid, 'current loop: p = ');
        fprintf(fid, ',%d', P(i));
        fprintf(fid, ';\n');
    for j = 1:5
        rand_if = zeros(1,50);
        mu = u(j)*ones(1,50);
        for b = 1:50
            X1 = mvnrnd(mu, eye(50), 20);
            X2 = mvnrnd(zeros(1,50),eye(50),20);
            X3 = mvnrnd(-mu, eye(50),20);
            Y = [X1;X2;X3];
            X0 = mvnrnd(zeros(1,p-50),eye(p-50),60);
            X = [Y,X0];
            Data = zscore(X)';
            [IFlabel, ~, ~] = ifpca(Data, 3);
            [~,RI,~,~]= RandIndex(c0,IFlabel');
            rand_if(b) = RI;
        end
        fprintf(fid, 'Rand index for IF-PCA:');
        for k = 1:50
            fprintf(fid, ',%f', rand_if(k));
        end
        fprintf(fid, ';\n');
    end
end
fclose(fid);
