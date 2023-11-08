load('shot15.mat','Vx','Vy','drReF')
drReF15 = drReF;
load('shot16.mat','drReF')
drReF16 = drReF;

load('shot19.mat','drReF')
drReF19 = drReF;
load('shot21.mat','drReF')
drReF21 = drReF;

levels = 10;
% levels = [-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05];%ìôçÇê¸ñ{êî
figure('Position',[700 150 550 500])
contourf(Vx, Vy, drReF15+drReF16-drReF19-drReF21, levels)%2éüå≥ìôçÇê¸ê}
colorbar
title('Reconstructed distribution function','FontSize',20);
xlabel('Vx [km/s]','FontSize',15);
ylabel('Vy [km/s]','FontSize',15);

