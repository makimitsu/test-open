function data_set(N_projection,gm2d,name, num_images)
% Set the number of images to generate

% Get size from gm2d matrix
[~, N_g] = size(gm2d);

% Calculate grid size
N_grid = sqrt(N_g);
m = N_grid;
n = N_grid;

% Grid for z and r coordinates
z = linspace(-1, 1, m);
r = linspace(-1, 1, n);

% Grids for plotting
z_grid = linspace(-200, 200, m);
r_grid = linspace(330, 70, n);

% Create directories for saving images if not exist
output_dir = ['output_images_', name];
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
initial_dir = [output_dir, '/initial_images_', name];
if ~exist(initial_dir, 'dir')
    mkdir(initial_dir);
end
projected_dir = [output_dir, '/projected_images_', name];
if ~exist(projected_dir, 'dir')
    mkdir(projected_dir);
end
withnoise_dir = [output_dir, '/withnoise_images_', name];
if ~exist(withnoise_dir, 'dir')
    mkdir(withnoise_dir);
end


% Metadata structure to store information about each image
metadata = struct('image_id', {}, 'num_sources', {}, 'source_positions', {}, 'source_sizes', {}, 'source_intensities', {});

for i = 1:num_images
    % Initialize the intensity matrix EE
    EE = zeros(m, n);

    % Randomize the number of light sources (e.g., between 1 and 10)
    num_sources = randi([1, 100]);

    % Preallocate arrays to store source properties
    source_positions = zeros(num_sources, 2); % [z_0, r_0]
    source_sizes = zeros(num_sources, 1);
    source_intensities = zeros(num_sources, 1);
    
    div = 4; %div*divに分割する
    zdiv = floor(rand()*div+1)/div*2-1;
    rdiv = floor(rand()*div+1)/div*2-1;
    for j = 1:num_sources
        % Randomize parameters for each light source
        if j < num_sources*1/2 %ある部分に発光を集中させる
            
            z_0 = zdiv - 2/div * rand();        % Random z-position
            r_0 = rdiv - 2/div * rand();        % Random r-position
        else
            z_0 = 1 - 2 * rand();        % Random z-position
            r_0 = 1 - 2 * rand();        % Random r-position
        end
        intensity = rand()*0.1; % Random intensity multiplier
        sigma = rand() * 0.1;     % Random size of the source

        % Store properties for metadata
        source_positions(j, :) = [z_0, r_0];
        source_sizes(j) = sigma;
        source_intensities(j) = intensity;

        % Calculate r_space and z_space
        [r_space, z_space] = meshgrid(r, z);

        % Calculate distances for the exponential terms
        r0_space = sqrt((z_space - z_0).^2 + (r_space - r_0).^2);

        % Update the intensity matrix EE with random parameters
        EE = EE + intensity * exp(-0.5 * (r0_space / sigma).^2);
    end

    % 背景の光を追加する
    background_intensity = 0.01;
    EE = EE + background_intensity;
    
    %正規化
    EE = EE./max(max(EE));
    EE = fliplr(rot90(EE)); %rが縦、zが横、右下最小
    
    %2D matrix is transformed to 1D transversal vector
    E = reshape(EE,1,[]);
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
    
    %非線形フィルター
    %[IIwgn,~] = imnlmfilt(IIwgn,'SearchWindowSize',25,'ComparisonWindowSize',15);

    %Save the image
    save(fullfile(initial_dir, sprintf('/image_%04d.mat', i)), 'EE');
    save(fullfile(projected_dir, sprintf('/image_%04d.mat', i)), 'II');
    save(fullfile(withnoise_dir, sprintf('/image_%04d.mat', i)), 'IIwgn');

    % Save metadata
    metadata(i).image_id = i;
    metadata(i).num_sources = num_sources;
    metadata(i).source_positions = source_positions;
    metadata(i).source_sizes = source_sizes;
    metadata(i).source_intensities = source_intensities;
end

% Save metadata to a MAT-file for future use
save([output_dir, '/metadata.mat'], 'metadata');


%可視化ファイル
convert_mat_to_png(initial_dir, [initial_dir, '/converted'], z_grid, r_grid, 'EE');
convert_mat_to_png(projected_dir, [projected_dir, '/converted'], z_grid, r_grid, 'II');
convert_mat_to_png(withnoise_dir, [withnoise_dir, '/converted'], z_grid, r_grid, 'IIwgn');

disp('Image generation and metadata creation completed.');
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