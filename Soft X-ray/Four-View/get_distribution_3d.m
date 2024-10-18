function EE = get_distribution_3d(M3d, K3d, U3d, s3d, v3d, VectorImage, N_grid)
%TP
Z=U3d'*VectorImage.';
gmin=-15;gmax=15;
% gmin=-30;gmax=0;
lg_gamma=gmin:1:gmax;
l_g = numel(lg_gamma);
gamma=10.^(lg_gamma);
V1 = zeros(1,l_g);
V2 = zeros(1,l_g);
Vgamma = zeros(1,21);
for n=1:l_g
    rho = M3d*gamma(n)./(s3d.^2+M3d*gamma(n));

    V11 = rho.*(Z.');
    V1(n)=M3d*sum(V11.^2);
    V2(n)=(sum(rho))^2;
    Vgamma(n)=V1(n)/V2(n);
end
[~,gamma_index]=min(Vgamma);
E = zeros(1,K3d);
for i=1:K3d
        if M3d>K3d
            v_1 = [v3d(i,:) zeros(1,M3d-K3d)];
        else
            v_1 = v3d(i,:);
        end
        E1 = (s3d./(s3d.^2+M3d*10^(lg_gamma(gamma_index)))).*v_1.*(Z.');
        E(i)=sum(E1);
end
EE = reshape(E,N_grid+1,N_grid+1,N_grid+1);
negativeEE = find(EE<0);
EE(negativeEE) = zeros(size(negativeEE));

end