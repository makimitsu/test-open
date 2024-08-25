function plot_SXR_test()
addpath '/Users/shohgookazaki/Documents/GitHub/test-open/Soft X-ray/Four-view';
NL = false;
plot_flag = true;

% 再構成条件の定義
newProjectionNumber = 50;%80; %投影数＝視線数の平方根
newGridNumber = 90;%100; %グリッド数（再構成結果の画素数の平方根）

[gm2d1, gm2d2, gm2d3, gm2d4, U1, U2, U3, U4, ...
          s1, s2, s3, s4, v1, v2, v3, v4, M, K, range, N_projection, N_grid] = parametercheck(newProjectionNumber, newGridNumber);


% number = (t-start)/interval+1;
% ファントムテスト用の画像を準備（4視点分）
[~,Iwgn1] = Assumption(N_projection,gm2d1,plot_flag);
%[~,Iwgn2] = Assumption(N_projection,gm2d2,false);
%[~,Iwgn3] = Assumption(N_projection,gm2d3,false);
%[~,Iwgn4] = Assumption(N_projection,gm2d4,false);
% % こっちを使う時は N_projection_new = 80, N_grid_new = 100
% [~,Iwgn1] = Assumption_2(N_projection,gm2d1,true);
% [~,Iwgn2] = Assumption_2(N_projection,gm2d2,false);
% [~,Iwgn3] = Assumption_2(N_projection,gm2d3,false);
% [~,Iwgn4] = Assumption_2(N_projection,gm2d4,false);


EE1 = get_distribution(M,K,gm2d1,U1,s1,v1,Iwgn1,plot_flag,NL);
%EE2 = get_distribution(M,K,gm2d2,U2,s2,v2,Iwgn2,plot_flag,NL);
%EE3 = get_distribution(M,K,gm2d3,U3,s3,v3,Iwgn3,plot_flag,NL);
%EE4 = get_distribution(M,K,gm2d4,U4,s4,v4,Iwgn4,plot_flag,NL);
%{
EE1 = get_MFI_reconstruction(Iwgn1, gm2d1);
EE2 = get_MFI_reconstruction(Iwgn2, gm2d2);
EE3 = get_MFI_reconstruction(Iwgn3, gm2d3);
EE4 = get_MFI_reconstruction(Iwgn4, gm2d4);
%}
% f = figure;
% f.Units = 'normalized';
% f.Position = [0.1,0.2,0.8,0.4];


% 表示範囲の設定に使うパラメータを取得
range = range./1000;
zmin = range(1);
zmax = range(2);
rmin = range(5);
rmax = range(6);
r_space_SXR = linspace(rmin,rmax,size(EE1,1));
z_space_SXR = linspace(zmin,zmax,size(EE1,2));

r_range = find(0.060<=r_space_SXR & r_space_SXR<=0.330);
r_space_SXR = r_space_SXR(r_range);
z_range = find(-0.17<=z_space_SXR & z_space_SXR<=0.17);

z_space_SXR = z_space_SXR(z_range);

EE1 = EE1(r_range,z_range);
%EE2 = EE2(r_range,z_range);
%EE3 = EE3(r_range,z_range);
%EE4 = EE4(r_range,z_range);

subplot(2,2,1);
[SXR_mesh_z,SXR_mesh_r] = meshgrid(z_space_SXR,r_space_SXR);
[~,h1] = contourf(SXR_mesh_z,SXR_mesh_r,EE1,20);
h1.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('1');
%{
subplot(2,2,2);
[~,h2] = contourf(SXR_mesh_z,SXR_mesh_r,EE2,20);
h2.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('2');

subplot(2,2,3);
[~,h3] = contourf(SXR_mesh_z,SXR_mesh_r,EE3,20);
h3.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('3');

subplot(2,2,4);
[~,h4] = contourf(SXR_mesh_z,SXR_mesh_r,EE4,20);
h4.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('4');
%}

end