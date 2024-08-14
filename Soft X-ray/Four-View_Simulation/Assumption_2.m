function [EE,Iwgn] = Assumption_2(N_projection,gm2d,plot_flag)
% 実際の画像を利用したテスト

[~, N_g] = size(gm2d);
% N_projection = sqrt(N_p);
N_grid = sqrt(N_g);
m=N_grid;
n=N_grid;
% z_0=0;
% r_0=-0.3;
% z=linspace(-1,1,m);
% r=linspace(-1,1,n);

% z_grid = linspace(200,-200,m);
z_grid = linspace(-200,200,m);
r_grid = linspace(330,70,n);

% [r_space,z_space] = meshgrid(r,z); %rが横、zが縦の座標系（左上最小）
% r0_space = sqrt((z_space-z_0).^2+(r_space-r_0).^2);
% r1_space = abs(0.5*(z_space-z_0)+(r_space-r_0))/sqrt(1.25);
% EE = exp(-0.5*r0_space.^2).*exp(-5*r1_space)+1.5*exp(-r0_space.^2);
% 
% % EE = zeros(m,n);
% % for i=1:m
% %     for j=1:n
% %         r0 = sqrt((z(i)-z_0)^2+(r(j)-r_0)^2);
% %         r1 = abs(0.5*(z(i)-z_0)+(r(j)-r_0))/sqrt(1.25);
% %         EE(i,j) = 1*exp(-0.5*r0^2)*exp(-5*r1)+1.5*exp(-r0^2);
% %     end
% % end
% % size(E)
% EE = EE./max(max(EE));
% EE = fliplr(rot90(EE)); %rが縦、zが横、右下最小

loadpath = '4_high.txt';
EE = readmatrix(loadpath); %縦がr、横がzで左上が最小になる
% EE = rot90(EE,2); %rが縦、zが横、右下最小
EE = flipud(EE);

% なぜか反転、おそらくファントムの生成に問題
% 要検証

% 2D matrix is transformed to 1D transversal vector
E = reshape(EE,1,[]);
% E = fliplr(reshape(EE,1,[]));
% whos gm2d
% whos E
I=gm2d*(E)';
Iwgn=awgn(I,10*log10(10),'measured'); % 5 related to 20%; 10 related to 10%;
Iwgn(Iwgn<0)=0;


% 1D column vector is transformed to 2D matrix
n_p = N_projection;
II = zeros(n_p);
IIwgn = zeros(n_p);
k=FindCircle(n_p/2);
II(k) = I;
IIwgn(k) = Iwgn;
Iwgn = Iwgn.';

% I_check = II(75,:);
% j = 1:numel(I_check);
% figure;plot(j,I_check);

if plot_flag
    [mesh_z,mesh_r] = meshgrid(z_grid,r_grid);
    % [~,h] = contourf(mesh_z,mesh_r,fliplr(EE),20);
    [~,h] = contourf(mesh_z,mesh_r,EE,20);
    h.LineStyle = 'none';
    c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
    % figure;imagesc(z_grid,r_grid,EE);c=colorbar('Ticks',[0.1,0.5,1]);
    % c.Label.String='Assumed Intensity [a.u.]';xlabel('Z [mm]');ylabel('R [mm]');
    % axis xy
    % ax = gca;
    % ax.XDir = 'reverse';
    figure;imagesc(II);c=colorbar('Ticks',[0,20,40]);
    c.Label.String='Assumed Intensity [a.u.]';xlabel('Z Pixels');ylabel('R Pixels');
    figure;imagesc(IIwgn);c=colorbar('Ticks',[0,20,40]);
    c.Label.String='Assumed Intensity [a.u.]';xlabel('Z Pixels');ylabel('R Pixels');
end

end

function k = FindCircle(L)
R = zeros(2*L);
for i = 1:2*L
    for j = 1:2*L
        R(i,j) = sqrt((L-i+0.5)^2+(j-L-0.5)^2);
    end
end
% figure;imagesc(R)
k = find(R<L);
end