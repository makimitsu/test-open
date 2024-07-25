function EE = get_distribution(M,K,gm2d,U,s,v,VectorImage,plot_flag,NL)

Z=U'*VectorImage.';

gmin=-15;gmax=15;
% gmin=-30;gmax=0;
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
    figure;semilogy(lg_gamma,Vgamma,'*');
    % figure;plot(lg_gamma,Vgamma,'*');
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

if NL %最大エントロピー
    gamma = 10^(lg_gamma(gamma_index));
    L = gm2d;
    f = zeros(K, 1);
    df = zeros(K, 1);
    eps = 0.1; % Newton法やめる時ようのε
    iter = 0;
    while (df.'*df > eps*(f.'*f) || iter == 0) && (iter < 15) % iter15以下で大丈夫なのかな
        disp(iter);
        Ex = diag(exp(f)); 
        invEx = Ex^(-1);%逆行列計算を楽にする
        Dx = M*gamma/2*Ex;
        invDx = Dx^(-1);%逆行列計算を楽にする
        Phif = M*gamma/2*(f+1)+ L.'*L*exp(f)-L.'*s.';
        A = eye(M)+L*invDx*L.'; %正定値対称行列になっているはずなんだよね
        b = L*invDx*Phif;
        xi = cholesky(A,b);
        
        df = -invEx*invDx*(Phif-L.'*xi);
        m = 0.5;
        df(df > m) = 1+m - m ./ df(df > m); % dfが大きくなりすぎるのを阻止。エラーの原因になる。
        df(df < -m) = -(1+m) - m ./df(df<-m); % dfが小さくなりすぎるのを阻止。エラーの原因になる。
        f = f+ df;
        iter = iter+1;
    end
    E = exp(f);

    EE = reshape(E.',sqrt(K),sqrt(K));
    
%{
if NL % 最小フィッシャー
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
    EE(EE<1e-5) = -1;
    % ここで10^-5未満を0に設定？
    W(diag_idx) = 1./EE;
    % W(W==Inf) = -1;
    W(W<0) = max(W, [], 'all');
    % if plot_flag
    %     figure;histogram(W(diag_idx));
    % end
    EE = (H' * H + (M * gamma) .* (C'* W * C))^(-1) * H' * G'; 

    EE = reshape(EE, sqrt(K), sqrt(K));
%}
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
    disp(max(E));
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
