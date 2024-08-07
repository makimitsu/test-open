% function plot_SXR_test()

NL = false;

% Definition of reconstruction condition
N_projection_new = 80; %square root of the projection number
N_grid_new = 100; %square root of the grid number

% load parameters for the reconstruction
filepath = 'parameters.mat';

% if the parameter file does not exist, calculate the parameters
if isfile(filepath)
    load(filepath, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
        's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');
    if N_projection_new ~= N_projection || N_grid_new ~= N_grid
        disp('Different parameters - Start calculation!');
        clc_parameters(N_projection_new,N_grid_new,filepath);
        load(filepath, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
            's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');
    end
else
    disp('No parameter - Start calculation!');
    clc_parameters(N_projection_new,N_grid_new,filepath);
    load(filepath, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
        's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');   
end
    
% number = (t-start)/interval+1;
plot_flag = false;
% % ファントムテスト用の画像を準備（4視点分）
% [~,Iwgn1] = Assumption(N_projection,gm2d1,true);
% [~,Iwgn2] = Assumption(N_projection,gm2d2,false);
% [~,Iwgn3] = Assumption(N_projection,gm2d3,false);
% [~,Iwgn4] = Assumption(N_projection,gm2d4,false);
% こっちを使う時は N_projection_new = 80, N_grid_new = 100

% prepare the phantom images for reconstruction test
[~,Iwgn_phantom1] = Assumption_2(N_projection,gm2d1,true);
[~,Iwgn_phantom2] = Assumption_2(N_projection,gm2d2,false);
[~,Iwgn_phantom3] = Assumption_2(N_projection,gm2d3,false);
[~,Iwgn_phantom4] = Assumption_2(N_projection,gm2d4,false);

% reconstruct original images from noisy signal
EE1 = clc_distribution(M,K,gm2d1,U1,s1,v1,Iwgn_phantom1,plot_flag,NL);
EE2 = clc_distribution(M,K,gm2d2,U2,s2,v2,Iwgn_phantom2,plot_flag,NL);
EE3 = clc_distribution(M,K,gm2d3,U3,s3,v3,Iwgn_phantom3,plot_flag,NL);
EE4 = clc_distribution(M,K,gm2d4,U4,s4,v4,Iwgn_phantom4,plot_flag,NL);


E1 = reshape(flipud(EE1),[],1);
E2 = reshape(flipud(EE2),[],1);
E3 = reshape(flipud(EE3),[],1);
E4 = reshape(flipud(EE4),[],1);

Iwgn1 = (gm2d1*E1).';
Iwgn2 = (gm2d2*E2).';
Iwgn3 = (gm2d3*E3).';
Iwgn4 = (gm2d4*E4).';

figure;
subplot(2,2,1);plot(Iwgn_phantom1);hold on;plot(Iwgn1);
subplot(2,2,2);plot(Iwgn_phantom2);hold on;plot(Iwgn2);
subplot(2,2,3);plot(Iwgn_phantom3);hold on;plot(Iwgn3);
subplot(2,2,4);plot(Iwgn_phantom4);hold on;plot(Iwgn4);

EE_new1 = clc_distribution(M,K,gm2d1,U1,s1,v1,Iwgn_phantom1,plot_flag,NL);
EE_new2 = clc_distribution(M,K,gm2d2,U2,s2,v2,Iwgn_phantom2,plot_flag,NL);
EE_new3 = clc_distribution(M,K,gm2d2,U2,s2,v2,Iwgn_phantom3,plot_flag,NL);
EE_new4 = clc_distribution(M,K,gm2d2,U2,s2,v2,Iwgn_phantom4,plot_flag,NL);


% f = figure;
% f.Units = 'normalized';
% f.Position = [0.1,0.2,0.8,0.4];


% plot the reconstructed images
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
EE2 = EE2(r_range,z_range);
EE3 = EE3(r_range,z_range);
EE4 = EE4(r_range,z_range);

subplot(2,2,1);
[SXR_mesh_z,SXR_mesh_r] = meshgrid(z_space_SXR,r_space_SXR);
[~,h1] = contourf(SXR_mesh_z,SXR_mesh_r,EE1,20);
h1.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('1');

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

EE_new1 = EE_new1(r_range,z_range);
EE_new2 = EE_new2(r_range,z_range);
EE_new3 = EE_new3(r_range,z_range);
EE_new4 = EE_new4(r_range,z_range);

subplot(2,2,1);
[SXR_mesh_z,SXR_mesh_r] = meshgrid(z_space_SXR,r_space_SXR);
[~,h1] = contourf(SXR_mesh_z,SXR_mesh_r,EE_new1,20);
h1.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('1');

subplot(2,2,2);
[~,h2] = contourf(SXR_mesh_z,SXR_mesh_r,EE_new2,20);
h2.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('2');

subplot(2,2,3);
[~,h3] = contourf(SXR_mesh_z,SXR_mesh_r,EE_new3,20);
h3.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('3');

subplot(2,2,4);
[~,h4] = contourf(SXR_mesh_z,SXR_mesh_r,EE_new4,20);
h4.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('4');

error1 = sum((EE_new1-EE1).^2/max(EE1,[],'all'),'all')/numel(EE1);
error2 = sum((EE_new2-EE2).^2/max(EE2,[],'all'),'all')/numel(EE2);
error3 = sum((EE_new3-EE3).^2/max(EE3,[],'all'),'all')/numel(EE3);
error4 = sum((EE_new4-EE4).^2/max(EE4,[],'all'),'all')/numel(EE4);

disp(error1);
disp(error2);
disp(error3);
disp(error4);



% end