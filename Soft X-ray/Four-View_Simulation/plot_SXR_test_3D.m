function plot_SXR_test_3D()
addpath '/Users/shohgookazaki/Documents/GitHub/test-open/Soft X-ray/Four-view';

newProjectionNumber = 30; %投影数＝視線数の平方根
newGridNumber = 20; %グリッド数

if evalin('base', 'exist(''N_projection'', ''var'')')
    NP = evalin('base', 'N_projection');
    if NP ~= newProjectionNumber
        [N_projection, N_grid3d, gm3d,U3d,s3d, v3d, M3d, K3d] = parametercheck_3d(newProjectionNumber, newGridNumber);
    end
else
    [N_projection, N_grid3d, gm3d,U3d,s3d, v3d, M3d, K3d] = parametercheck_3d(newProjectionNumber, newGridNumber);
end

proj = Assumption3d(N_projection,gm3d,N_grid3d,true);

EE = get_distribution_3d(M3d, K3d, U3d, s3d, v3d, proj.', N_grid3d);

EE = permute(EE, [2, 1, 3]);

rmin3d = 0; 
rmax3d = 375;
dmin3d = -375;
dmax3d = 375;
zmin3d = -300;
zmax3d = 300;
r = linspace(rmin3d, rmax3d, newGridNumber+1);
z = linspace(zmin3d, zmax3d, newGridNumber+1);
d = linspace(dmin3d, dmax3d, newGridNumber+1);

% Create a grid for the data
[Z, R, D] = meshgrid(z, r, d);

figure;
slice(Z, R, D, EE, [-100,100], [100,200], [-100,100]); % Example slice positions along Z

xlabel('z軸');
ylabel('r軸');
zlabel('depth');
title(strcat('3D visualization of EE'));

% Adjust color and viewing angle
colormap(jet); % Use jet colormap
colorbar;      % Add a color bar
view(3);       % Set to 3D view
end