function EE = get_distribution(M,K,gm2d,U,s,v,VectorImage,plot_flag,NonLiner_flag)
% この関数は，一つのベクトルイメージから二次元発光分布を再構成します．

Z=U'*VectorImage.';

gmin=-15;gmax=15;
lg_gamma=gmin:1:gmax;
l_g = numel(lg_gamma);
gamma=10.^(lg_gamma);

V1 = zeros(1,l_g);
V2 = zeros(1,l_g);
Vgamma = zeros(1,21);

for n=1:l_g
    rho = M*gamma(n)./(s.^2+M*gamma(n));

    V11 = rho.*(Z.');
    V1(n)=M*sum(V11.^2);
    V2(n)=(sum(rho))^2;
    Vgamma(n)=V1(n)/V2(n);
end

if plot_flag
    figure;plot(lg_gamma,Vgamma,'*');
    xlabel('logγ');
    ylabel('GCV');
end

[~,gamma_index]=min(Vgamma);
E = zeros(1,K);

% for i=1:K
%     if M>K
%         v_1 = [v(i,:) zeros(1,M-K)];
%     else
%         v_1 = v(i,:);
%     end
%     E1 = (s./(s.^2+M*10^(lg_gamma(gamma_index)))).*v_1.*(Z.');
%     E(i)=sum(E1);
% end
% EE = reshape(E,sqrt(K),sqrt(K)); %ここで縦がr、横がzで左下が最小になる

if NonLiner_flag
    gamma = 10^(lg_gamma(gamma_index));
    C = Laplacian(sqrt(K)-1);
    H = gm2d;
    G = VectorImage;
    W = eye(size(C, 1));
    diag_idx = find(W);
    for i=1:K
        % v_1 = v(i,:);
        if M>K
            v_1 = [v(i,:) zeros(1,M-K)];
        else
            v_1 = v(i,:);
        end
        E1 = (s./(s.^2+M*10^(lg_gamma(gamma_index)))).*v_1.*(Z.');
        E(i)=sum(E1);
    end
    EE = E;
    W(diag_idx) = 1./EE;
    W(W==Inf) = -1;
    W(W<0) = max(W, [], 'all');
    EE = (H' * H + (M * gamma) .* (C'* W * C))^(-1) * H' * G'; 

    EE = reshape(EE, sqrt(K), sqrt(K));
else
    for i=1:K
        if M>K
            v_1 = [v(i,:) zeros(1,M-K)];
        else
            v_1 = v(i,:);
        end
        E1 = (s./(s.^2+M*10^(lg_gamma(gamma_index)))).*v_1.*(Z.');
        E(i)=sum(E1);
    end
    EE = reshape(E,sqrt(K),sqrt(K)); %ここで縦がr、横がzで左下が最小になる
end

% contourfで単調増加する軸から生成されたmeshgridを使ってプロットすると上下が反転する
EE = flipud(EE);
% figure;imagesc(EE);

end

function C = Laplacian(N_grid)
k=N_grid+1;
K=k*k;
C=zeros(K);
for i=1:1:k
    for j=1:1:k
           C((i-1)*k+j,(i-1)*k+j)=-4;
        if j+1<=k
            C((i-1)*k+j,(i-1)*k+j+1)=1;
        end
        
        if j-1>=1
            C((i-1)*k+j,(i-1)*k+j-1)=1;
        end
        
        if i-1-1>=0
            C((i-1)*k+j,(i-1-1)*k+j)=1;
        end
        
        if i-1+1<=k-1
            C((i-1)*k+j,(i-1+1)*k+j)=1;
        end
    end
end
end