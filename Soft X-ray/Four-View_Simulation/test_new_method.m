addpath('/Users/shinjirotakeda/Documents/GitHub/test-open/Soft X-ray/Four-View');
SXRfilename = '/Users/shinjirotakeda/Library/CloudStorage/OneDrive-TheUniversityofTokyo/Documents/SXR_Images/230920/shot056.tif';
date = 230920;
number = 6;
doPlot = false;
doNLR = true;

parameterFile = 'parameters.mat';
load(parameterFile, 'gm2d1', 'gm2d2', 'gm2d3', 'gm2d4', 'U1', 'U2', 'U3', 'U4', ...
            's1', 's2', 's3', 's4', 'v1', 'v2', 'v3', 'v4', 'M', 'K', 'range','N_projection', 'N_grid');
im = imread(SXRfilename);
[f1,f2,f3,f4] = get_sxr_image(date,number,N_projection,im);

EE1 = get_distribution(M,K,gm2d1,U1,s1,v1,f1,doPlot,doNLR);
g_old1 = reshape(EE1,[],1);
g_new1 = get_MFI_reconstruction(f1,gm2d1);

EE_new1 = reshape(g_new1,sqrt(K),sqrt(K));
EE_new1 = flipud(EE_new1);
EE_new1(EE_new1<0) = 0;

figure;
plot(g_old1);
hold on
plot(g_new1);

figure;
subplot(1,2,1);imagesc(EE1);colorbar;clim([0 0.5]);
subplot(1,2,2);imagesc(EE_new1);colorbar;clim([0 0.5]);

