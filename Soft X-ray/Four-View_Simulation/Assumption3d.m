function proj = Assumption3d(N_projection, gm3d, N_grid3d, plotflag)
    m = N_grid3d+1; % z dimension
    n = N_grid3d+1; % r dimension
    p = N_grid3d+1; % y dimension
    
    z = linspace(-1, 1, m);
    r = linspace(-1, 1, n);
    y = linspace(-1, 1, p); % New y dimension
    
    [r_space, z_space, y_space] = meshgrid(r, z, y); % 3D meshgrid
    
    % Define the centers and properties for the light sources in 3D
    % centers = [% Center of 3D light sources, with z_center, r_center, y_center, semi_major_axis, semi_minor_axis, intensity
    %     0, -0.7, -0.3, 0.2, 0.2, 0.8;% Example source
    %     -0.7,-0.7,0,0.2,0.2,0.8;
    %     ];

    centers = [
        0,-0.8,0,0.5,0.5,0.7;
        0,0.8,0,0.3,0.3,0.4;
        ];
    % Add more configurations as needed
        
    % Initialize EE matrix for 3D
    EE = zeros(size(r_space));
        
    % Add multiple light sources in 3D (ellipsoidal)
    for i = 1:size(centers, 1)
        z_center = centers(i, 1);
        r_center = centers(i, 2);
        y_center = centers(i, 3);
        a = centers(i, 4); % semi-major axis
        b = centers(i, 5); % semi-minor axis
        intensity = centers(i, 6);
            
        % Define the ellipsoidal light source in 3D
        EE = EE + intensity * exp(-((z_space - z_center).^2 / (2 * a^2) + ...
                                     (r_space - r_center).^2 / (2 * b^2) + ...
                                     (y_space - y_center).^2 / (2 * b^2)));
    end
        
    % Normalize
    EE = EE ./ max(max(max(EE)));
    
    E = reshape(EE,1,[]);
    I = gm3d*(E)';
    n = 10;
    Iwgn = awgn(I,10*log10(100/n), 'measured');
    Iwgn(Iwgn<0) = 0;
    
    I = reshape(I,[],4);
    Iwgn = reshape(Iwgn,[],4);
    n_p = N_projection;
    II = zeros(n_p,n_p,4);
    IIwgn = zeros(n_p,n_p,4);
    k = FindCircle(n_p/2);
    for i = 1:4
        III = II(:,:,i);
        IIIwgn = IIwgn(:,:,i);
        III(k) = I(:,i);
        IIIwgn(k) = Iwgn(:,i);
        II(:,:,i) = III;
        IIwgn(:,:,i) = IIIwgn;
    end
    
    if plotflag
        figure;
        for i = 1:4
            subplot(2, 2, i);       % 2行2列のグリッドにプロット
            imagesc(IIwgn(:, :, i));    % i番目のスライスをプロット
            colorbar;               % カラーバーを追加
            title(['Slice ', num2str(i)]);
        end
    end

    if plotflag
        EE = permute(EE, [2, 1, 3]);
        x = 1:21;
        y = 1:21;
        z = 1:21;
        [X, Y, Z] = meshgrid(x, y, z);
        figure;
        slice(X, Y, Z, EE, [5,10], [10,20], [5,10]); % Example slice positions along Z
        
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

    proj = cat(1, Iwgn(:,1),Iwgn(:,2),Iwgn(:,3),Iwgn(:,4));

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