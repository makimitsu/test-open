function plot_SXR_test_3D()
addpath '/Users/shohgookazaki/Documents/GitHub/test-open/Soft X-ray/Four-view';

newProjectionNumber = 50;%80; %投影数＝視線数の平方根
newGridNumber = 90;%100; %グリッド数（再構成結果の画素数の平方根）

if evalin('base', 'exist(''N_projection'', ''var'')')
    NP = evalin('base', 'N_projection');
    if NP ~= newProjectionNumber
        [gm2d1, ~, gm2d3, gm2d4, U1, U2, U3, U4, ...
                  s1, s2, s3, s4, v1, v2, v3, v4, M, K, range, N_projection, N_grid,gm3d,U3d,s3d, v3d, M3d, K3d] = parametercheck(newProjectionNumber, newGridNumber);
    end
else
    [gm2d1, gm2d2, gm2d3, gm2d4, U1, U2, U3, U4, ...
              s1, s2, s3, s4, v1, v2, v3, v4, M, K, range, N_projection, N_grid, gm3d,U3d,s3d, v3d, M3d, K3d] = parametercheck(newProjectionNumber, newGridNumber);
end
N_grid = 15;
proj = Assumption3d(N_projection,gm3d,true);

EE = get_distribution_3d(M3d, K3d, U3d, s3d, v3d, proj.', N_grid);

x = 1:16;
y = 1:16;
z = 1:31;

% Create a grid for the data
[X, Y, Z] = meshgrid(x, y, z);

% Plot slices at specific positions along the Z-axis
figure;
slice(X, Y, Z, EE, [5,10], [5,10], [10,20]); % Example slice positions along Z

% Add labels and title
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Visualization of EE');

% Adjust color and viewing angle
colormap(jet); % Use jet colormap
colorbar;      % Add a color bar
view(3);       % Set to 3D view
end