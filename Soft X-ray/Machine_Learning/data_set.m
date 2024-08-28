addpath '/Users/shohgookazaki/Documents/GitHub/test-open/Soft X-ray/Four-view';
Newdata = true;

num_images = 5; % 決める

% 再構成条件の定義
newProjectionNumber = 50; %投影数＝視線数の平方根
newGridNumber = 90; %グリッド数（再構成結果の画素数の平方根）

if evalin('base', 'exist(''N_projection'', ''var'')')
    NP = evalin('base', 'N_projection');
    if NP ~= newProjectionNumber
        [gm2d1, gm2d2, gm2d3, gm2d4, U1, U2, U3, U4, ...
                  s1, s2, s3, s4, v1, v2, v3, v4, M, K, range, N_projection, N_grid] = parametercheck(newProjectionNumber, newGridNumber);
    end
else
    [gm2d1, gm2d2, gm2d3, gm2d4, U1, U2, U3, U4, ...
              s1, s2, s3, s4, v1, v2, v3, v4, M, K, range, N_projection, N_grid] = parametercheck(newProjectionNumber, newGridNumber);
end

%ここからファントム生成

% Get size from gm2d matrix
[~, N_g] = size(gm2d1);

% Calculate grid size
N_grid = sqrt(N_g);
m=N_grid;
n=N_grid;
z_0=0;
r_0=-0.3;
z=linspace(-1,1,m);
r=linspace(-1,1,n);
z_grid = linspace(-200,200,m);
r_grid = linspace(330,70,n);

[r_space,z_space] = meshgrid(r,z);

% Grid for z and r coordinates
z = linspace(-1, 1, m);
r = linspace(-1, 1, n);

% Grids for plotting
z_grid = linspace(-200, 200, m);
r_grid = linspace(330, 70, n);

% Create directories for saving images if not exist
output_dir = 'test_output_images';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
initial_dir = [output_dir, '/initial_images'];
if ~exist(initial_dir, 'dir')
    mkdir(initial_dir);
end
projected_dir = [output_dir, '/projected_images'];
if ~exist(projected_dir, 'dir')
    mkdir(projected_dir);
end
withnoise_dir = [output_dir, '/withnoise_images'];
if ~exist(withnoise_dir, 'dir')
    mkdir(withnoise_dir);
end

EETikhonov_dir = [output_dir, '/EETikhonov_images'];
if ~exist(EETikhonov_dir, 'dir')
    mkdir(EETikhonov_dir);
end



% Metadata structure to store information about each image
metadata = struct('image_id', {}, 'num_sources', {}, 'source_positions', {}, 'source_sizes', {}, 'source_intensities', {});

for i = 1:num_images
    % Initialize the intensity matrix EE
    EE = zeros(size(r_space));

    % Randomize the number of light sources (e.g., between 1 and 10)
    num_sources = randi([1, 50]);

    % Preallocate arrays to store source properties
    source_positions = zeros(num_sources, 2); % [z_0, r_0]
    source_sizes = zeros(num_sources, 1);
    source_intensities = zeros(num_sources, 1);
    
    div = rand()*8; %div*divに分割する
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
        
        if mod(num_images, 2) == 0
            intensity = rand()*10;
            sigma = rand();
        else
            intensity = rand()*0.1; % Random intensity multiplier
            sigma = rand() * 0.1;     % Random size of the source
        end

        % Store properties for metadata
        source_positions(j, :) = [z_0, r_0];
        source_sizes(j) = sigma;
        source_intensities(j) = intensity;
        
        % Calculate distances for the exponential terms
        r0_space = sqrt((z_space - z_0).^2 + (r_space - r_0).^2);

        % Update the intensity matrix EE with random parameters
        EE = EE + intensity * exp(-0.5 * (r0_space / sigma).^2);
    end

    % 背景の光を追加する
    %background_intensity = 0.01;
    %EE = EE + background_intensity;
    
    %正規化
    EE = EE./max(max(EE));
    %EE = fliplr(rot90(EE)); %rが縦、zが横、右下最小
    
    %2D matrix is transformed to 1D transversal vector
    E = reshape(EE,1,[]);
    I1=gm2d1*(E)';
    I2=gm2d2*(E)';
    I3=gm2d3*(E)';
    I4=gm2d4*(E)';

    Iwgn1=awgn(I1,10*log10(10),'measured'); % 5 related to 20%; 10 related to 10%;
    Iwgn2=awgn(I2,10*log10(10),'measured');
    Iwgn3=awgn(I3,10*log10(10),'measured');
    Iwgn4=awgn(I4,10*log10(10),'measured');

    Iwgn1(Iwgn1<0)=0;
    Iwgn2(Iwgn2<0)=0;
    Iwgn3(Iwgn3<0)=0;
    Iwgn4(Iwgn4<0)=0;

    
    
    % 1D column vector is transformed to 2D matrix
    n_p = N_projection;

    II1 = zeros(n_p);
    II2 = zeros(n_p);
    II3 = zeros(n_p);
    II4 = zeros(n_p);

    IIwgn1 = zeros(n_p);
    IIwgn2 = zeros(n_p);
    IIwgn3 = zeros(n_p);
    IIwgn4 = zeros(n_p);

    k=FindCircle(n_p/2);
    II1(k) = I1;
    II2(k) = I2;
    II3(k) = I3;
    II4(k) = I4;

    IIwgn1(k) = Iwgn1;
    IIwgn2(k) = Iwgn2;
    IIwgn3(k) = Iwgn3;
    IIwgn4(k) = Iwgn4;

    Iwgn1 = Iwgn1.';
    Iwgn2 = Iwgn2.';
    Iwgn3 = Iwgn3.';
    Iwgn4 = Iwgn4.';
    
    %非線形フィルター
    %[IIwgn,~] = imnlmfilt(IIwgn,'SearchWindowSize',25,'ComparisonWindowSize',15);
    
    EEt1 = get_distribution(M,K,gm2d1,U1,s1,v1,Iwgn1,false,false);
    EEt2 = get_distribution(M,K,gm2d2,U2,s2,v2,Iwgn2,false,false);
    EEt3 = get_distribution(M,K,gm2d3,U3,s3,v3,Iwgn3,false,false);
    EEt4 = get_distribution(M,K,gm2d4,U4,s4,v4,Iwgn4,false,false);

    EEt1 = flipud(EEt1);
    EEt2 = flipud(EEt2);
    EEt3 = flipud(EEt3);
    EEt4 = flipud(EEt4);
    
    %Save
    EE1 = EE; %学習時の名前の統一。全部同じデータ入っているけど
    EE2 = EE;
    EE3 = EE;
    EE4 = EE;
    
    sxr1 = IIwgn1; % 学習時の名前の統一
    sxr2 = IIwgn2;
    sxr3 = IIwgn3;
    sxr4 = IIwgn4;

    save(fullfile(initial_dir, sprintf('/image_%04d.mat', i)), 'EE1','EE2','EE3','EE4');
    save(fullfile(projected_dir, sprintf('/image_%04d.mat', i)), 'II1','II2','II3','II4');
    save(fullfile(withnoise_dir, sprintf('/image_%04d.mat', i)), 'sxr1','sxr2','sxr3','sxr4');
    save(fullfile(EETikhonov_dir, sprintf('/image_%04d.mat', i)), 'EEt1','EEt2','EEt3','EEt4');

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
%disp('converting to png...');
convert_mat_to_png(initial_dir, [initial_dir, '/converted'], z_grid, r_grid, 'EE1');
%convert_mat_to_png(projected_dir, [projected_dir, '/converted'], z_grid, r_grid, 'II1');
%convert_mat_to_png(withnoise_dir, [withnoise_dir, '/converted'], z_grid, r_grid, 'sxr1');
convert_mat_to_png(EETikhonov_dir, [EETikhonov_dir, '/converted'], z_grid, r_grid, 'EEt1');

disp('Image generation and metadata creation completed.');


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