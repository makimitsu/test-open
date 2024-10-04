function [] = get_parameters(N_projection,N_grid,filepath)

% filepath = '/Users/shinjirotakeda/Documents/GitHub/test-open/Soft X-ray/Four-View_Simulation/parameters.mat';

% 視線の分布、重み行列の作成
zhole1=40;zhole2=-40;                                  
zmin1=-240;zmax1=320;zmin2=-320;zmax2=240;             
rmin=55;rmax=375;
range = [zmin1,zmax1,zmin2,zmax2,rmin,rmax];            
l1 = MCPLine_up(N_projection,zhole2,false);
gm2d1 = LineProjection(l1,N_grid,zmin2,zmax2,rmin,rmax,false,true); 
l2 = MCPLine_down(N_projection,zhole2,false);
gm2d2 = LineProjection(l2,N_grid,zmin2,zmax2,rmin,rmax,false,false);
l3 = MCPLine_up(N_projection,zhole1,false);
gm2d3 = LineProjection(l3,N_grid,zmin1,zmax1,rmin,rmax,false,true);
l4 = MCPLine_down(N_projection,zhole1,false);
gm2d4 = LineProjection(l4,N_grid,zmin1,zmax1,rmin,rmax,false,false);

% plot_overlapping_rays(l1,l2,l3,l4,5);
l = {l1,l2,l3,l4};
N_grid3d = 25;
N_grid_x = N_grid3d; N_grid_y = 2*N_grid3d; N_grid_z = N_grid3d;
zmin3d = -400; zmax3d = 400; ymin3d = -400; ymax3d = 100; xmin3d = -100; xmax3d = 50;
gm3d = LineProjection3D(l, N_grid_x, N_grid_y, N_grid_z, xmin3d, xmax3d, ymin3d, ymax3d, zmin3d, zmax3d, false);

disp('Calculation for TP starting...')
% ラプラシアン行列の計算と特異値分解
C = Laplacian(N_grid);
C3d = Laplacian_3D(N_grid3d);
[U1,S1,V1]=svd(gm2d1*(C^(-1)),'econ');
[U2,S2,V2]=svd(gm2d2*(C^(-1)),'econ');
[U3,S3,V3]=svd(gm2d3*(C^(-1)),'econ');
[U4,S4,V4]=svd(gm2d4*(C^(-1)),'econ');
[U3d,S3d,V3d] = svd(gm3d*(C3d^(-1)), 'econ');
% [U1,S1,V1]=svd(gm2d1*(C^(-1)));
% [U2,S2,V2]=svd(gm2d2*(C^(-1)));
% [U3,S3,V3]=svd(gm2d3*(C^(-1)));
% [U4,S4,V4]=svd(gm2d4*(C^(-1)));
v1=(C^(-1)*V1);
v2=(C^(-1)*V2);
v3=(C^(-1)*V3);
v4=(C^(-1)*V4);
v3d = (C3d^(-1)*V3d);

[M,K] = size(gm2d1);
[M3d,K3d] = size(gm3d);
if K>M
    v1 = v1(:,1:M);
    v2 = v2(:,1:M);
    v3 = v3(:,1:M);
    v4 = v4(:,1:M);
    v3d = v3d(:,1,M3d);
end
s1 = (diag(S1)).';
s2 = (diag(S2)).';
s3 = (diag(S3)).';
s4 = (diag(S4)).';
s3d = (diag(S3d)).';
if M>K
    s1 = [s1 zeros(1,M-K)];
    s2 = [s2 zeros(1,M-K)];
    s3 = [s3 zeros(1,M-K)];
    s4 = [s4 zeros(1,M-K)];
    s3d = [s3d zeros(1,M3d-K3d)];
end

save(filepath,'gm2d1','gm2d2','gm2d3','gm2d4', ...
    'U1','U2','U3','U4','s1','s2','s3','s4', ...
    'v1','v2','v3','v4','M','K','range','N_projection','N_grid','gm3d','U3d','s3d','v3d','M3d','K3d','-v7.3');

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

function l = MCPLine_up(N_projection,Z_hole,plot_flag)
d_hole = 24.4; % distance between the hole and the MCP
r_mcp=10;  % radius of the MCP plate

% Y_hole = 427.85+12; 
Y_hole = 413.24+12; %?

% X_hole = 209.62;
X_hole = 208.13; %?

Y_initial=Y_hole+d_hole;  X_initial=X_hole-r_mcp;  Z_initial=Z_hole+r_mcp;    %CCD position
X_end=X_hole+r_mcp;      Z_end=Z_hole-r_mcp;

Nh=N_projection-1; 
Dhx=(X_end-X_initial)/Nh;
Dhz=(Z_end-Z_initial)/Nh;

X=X_initial:Dhx:X_end;
Z=Z_initial:Dhz:Z_end;
Y=repelem(Y_initial,N_projection);

r_center = 55;
r_device = 375;
ll(N_projection,N_projection) = struct('x',[],'y',[],'z',[]);
if plot_flag
    f1=figure;
    f1.Name = "Upper MCP Line - Lightlines";
    f2=figure;
    f2.Name = "Upper MCP Line - 視線の行列";
end
for i=1:N_projection
    for j=1:N_projection
        ll(i,j).y=Y(j):-10:-400;
        length = numel(ll(i,j).y);
        ll(i,j).x=(ll(i,j).y-Y_hole)*(X(i)-X_hole)/(Y(j)-Y_hole)+X_hole;
        ll(i,j).z=(ll(i,j).y-Y_hole)*(Z(j)-Z_hole)/(Y(j)-Y_hole)+Z_hole;
        %中心軸で視線が遮られることを考慮
        r = sqrt(ll(i,j).y.^2+ll(i,j).x.^2);
        A = find(r<=r_center);
        if isempty(A) == 0
            obs1 = A(1)-1;
            ll(i,j).x = [repelem(ll(i,j).x(1),length-obs1) ll(i,j).x(1:obs1)];
            ll(i,j).y = [repelem(ll(i,j).y(1),length-obs1) ll(i,j).y(1:obs1)];
            ll(i,j).z = [repelem(ll(i,j).z(1),length-obs1) ll(i,j).z(1:obs1)];
        end
        B = find(r>=r_device & ll(i,j).y<=0);
        if isempty(B) == 0
            obs2 = B(1)-1;
            ll(i,j).x = [repelem(ll(i,j).x(1),length-obs2) ll(i,j).x(1:obs2)];
            ll(i,j).y = [repelem(ll(i,j).y(1),length-obs2) ll(i,j).y(1:obs2)];
            ll(i,j).z = [repelem(ll(i,j).z(1),length-obs2) ll(i,j).z(1:obs2)];
        end
        %plot the lightline
        if plot_flag
            figure(f1);plot3(X(i),Y(j),Z(j),'*',ll(i,j).x,ll(i,j).y,ll(i,j).z);  
            hold on;grid on; 
        end
    end
end

%視線の行列のうち円の内部に含まれるものだけをベクトル化
L = N_projection/2;
k = FindCircle(L);
l = ll(k);
if plot_flag
    for i = 1:numel(l)
        figure(f2);plot3(l(i).x,l(i).y,l(i).z);
        hold on;grid on;
    end
end
end

function l = MCPLine_down(N_projection,Z_hole,plot_flag)
d_hole =24.4; % distance between the hole and the MCP
r_mcp=10;  %radius of the MCP plate
% Y_hole=425.24;  X_hole=195.13;          % position of the hole

% Y_hole = 427.85+12; %
Y_hole = 413.24+12; %?

% X_hole = 145.62; 
X_hole = 144.13; %?

Y_initial=Y_hole+d_hole;  X_initial=X_hole-r_mcp;  Z_initial=Z_hole+r_mcp;    %CCD position
X_end=X_hole+r_mcp;      Z_end=Z_hole-r_mcp;

Nh=N_projection-1; 
Dhx=(X_end-X_initial)/Nh;
Dhz=(Z_end-Z_initial)/Nh;

X=X_initial:Dhx:X_end;
Z=Z_initial:Dhz:Z_end;
Y=repelem(Y_initial,N_projection);

r_center = 55;
r_device = 375;
ll(N_projection,N_projection) = struct('x',[],'y',[],'z',[]);
if plot_flag
    f1=figure;
    f1.Name = "Lower MCP Line - Lightlines";
    f2=figure;
    f2.Name = "Lower MCP Line - 視線の行列";
end
for i=1:N_projection
    for j=1:N_projection
        ll(i,j).y=Y(j):-10:-400; %大から小
        length = numel(ll(i,j).y);
        ll(i,j).x=(ll(i,j).y-Y_hole)*(X(i)-X_hole)/(Y(j)-Y_hole)+X_hole;
        ll(i,j).z=(ll(i,j).y-Y_hole)*(Z(j)-Z_hole)/(Y(j)-Y_hole)+Z_hole;
        %中心軸で視線が遮られることを考慮
        r = sqrt(ll(i,j).y.^2+ll(i,j).x.^2);
        A = find(r<=r_center);
        if isempty(A) == 0
            obs1 = A(1)-1;
            ll(i,j).x = [repelem(ll(i,j).x(1),length-obs1) ll(i,j).x(1:obs1)];
            ll(i,j).y = [repelem(ll(i,j).y(1),length-obs1) ll(i,j).y(1:obs1)];
            ll(i,j).z = [repelem(ll(i,j).z(1),length-obs1) ll(i,j).z(1:obs1)];
        end
        B = find(r>=r_device & ll(i,j).y<=0);
        if isempty(B) == 0
            obs2 = B(1)-1;
            ll(i,j).x = [repelem(ll(i,j).x(1),length-obs2) ll(i,j).x(1:obs2)];
            ll(i,j).y = [repelem(ll(i,j).y(1),length-obs2) ll(i,j).y(1:obs2)];
            ll(i,j).z = [repelem(ll(i,j).z(1),length-obs2) ll(i,j).z(1:obs2)];
        end
        %plot the lightline
        if plot_flag
            figure(f1);plot3(X(i),Y(j),Z(j),'*',ll(i,j).x,ll(i,j).y,ll(i,j).z);  
            hold on;grid on; 
        end
    end
end

%視線の行列のうち円の内部に含まれるものだけをベクトル化
L = N_projection/2;
k = FindCircle(L);
l = ll(k);
if plot_flag
    for i = 1:numel(l)
        figure(f2);plot3(l(i).x,l(i).y,l(i).z);
        hold on;grid on;
    end
end
end

function l = MCPLine(N_projection,Z_hole,plot_flag)
d_hole =24.4; % distance between the hole and the MCP
r_mcp=10;  %radius of the MCP plate
Y_hole=425.24;  X_hole=195.13;          % position of the hole
% Y_hole_new = 427.85+12; 
% X_hole_new_up = 209.62; X_hole_new_down = X_hole_new_up - 64;
Y_initial=Y_hole+d_hole;  X_initial=X_hole-r_mcp;  Z_initial=Z_hole+r_mcp;    %CCD position
X_end=X_hole+r_mcp;      Z_end=Z_hole-r_mcp;

Nh=N_projection-1; 
Dhx=(X_end-X_initial)/Nh;
Dhz=(Z_end-Z_initial)/Nh;

X=X_initial:Dhx:X_end;
Z=Z_initial:Dhz:Z_end;
Y=repelem(Y_initial,N_projection);

r_center = 55;
r_device = 375;
ll(N_projection,N_projection) = struct('x',[],'y',[],'z',[]);
if plot_flag
    f1=figure;
    f2=figure;
end
for i=1:N_projection
    for j=1:N_projection
        ll(i,j).y=Y(j):-10:-400;
        length = numel(ll(i,j).y);
        ll(i,j).x=(ll(i,j).y-Y_hole)*(X(i)-X_hole)/(Y(j)-Y_hole)+X_hole;
        ll(i,j).z=(ll(i,j).y-Y_hole)*(Z(j)-Z_hole)/(Y(j)-Y_hole)+Z_hole;
        %中心軸で視線が遮られることを考慮
        r = sqrt(ll(i,j).y.^2+ll(i,j).x.^2);
        A = find(r<=r_center);
        if isempty(A) == 0
            obs1 = A(1)-1;
            ll(i,j).x = [repelem(ll(i,j).x(1),length-obs1) ll(i,j).x(1:obs1)];
            ll(i,j).y = [repelem(ll(i,j).y(1),length-obs1) ll(i,j).y(1:obs1)];
            ll(i,j).z = [repelem(ll(i,j).z(1),length-obs1) ll(i,j).z(1:obs1)];
        end
        B = find(r>=r_device & ll(i,j).y<=0);
        if isempty(B) == 0
            obs2 = B(1)-1;
            ll(i,j).x = [repelem(ll(i,j).x(1),length-obs2) ll(i,j).x(1:obs2)];
            ll(i,j).y = [repelem(ll(i,j).y(1),length-obs2) ll(i,j).y(1:obs2)];
            ll(i,j).z = [repelem(ll(i,j).z(1),length-obs2) ll(i,j).z(1:obs2)];
        end
        %plot the lightline
        if plot_flag
            figure(f1);plot3(X(i),Y(j),Z(j),'*',ll(i,j).x,ll(i,j).y,ll(i,j).z);  
            hold on;grid on; 
        end
    end
end

%視線の行列のうち円の内部に含まれるものだけをベクトル化
L = N_projection/2;
k = FindCircle(L);
l = ll(k);
if plot_flag
    for i = 1:numel(l)
        figure(f2);plot3(l(i).x,l(i).y,l(i).z);
        hold on;grid on;
    end
end
end

function gm2d = LineProjection(l,N_grid,zmin,zmax,rmin,rmax,plot_flag,up_flag)
    % rmin=70;rmax=330;
    % rmin=55;rmax=375;
    % rmin=70;rmax=280;
    
    N_p = numel(l);
    N_g = N_grid+1;
    DR=(rmax-rmin)/N_grid;
    DZ=(zmax-zmin)/N_grid;
    gm2d = zeros(N_p,N_g^2);
    
    if plot_flag
        f1=figure;
        f2=figure;
        if up_flag
            f1.Name = "Upper Line Projection";
            f2.Name = "Upper Line Projection";
        else
            f1.Name = "Lower Line Projection";
            f2.Name = "Lower Line Projection";
        end
    end
    
    for i = 1:N_p
        %各視線の座標からrz座標を計算、プロット
        x=l(i).x;
        y=l(i).y;
        z=l(i).z;
        r=sqrt(x.^2+y.^2);
        if plot_flag
            figure(f1);plot(z,r);
            hold on;grid on;
            xlabel('Z [mm]');ylabel('R [mm]');
        end
        
        %再構成対象の領域内の視線を抽出
        k=find(r>=rmin&r<=rmax&z>=zmin&z<=zmax);
        pl_r=r(k);
        pl_z=z(k);
        if plot_flag
            figure(f2);plot(pl_z,pl_r,'.');
            hold on;grid on;
            xlabel('Z [mm]');ylabel('R [mm]');
        end
        
        %各点のグリッド座標を求め、各グリッド毎に含まれる点の数を数え上げる
        r_grid = fix((pl_r-rmin)./DR)+1;%グリッド座標
        z_grid = fix((pl_z-zmin)./DZ)+1;
        num_p = numel(r_grid);
        gm_temp = zeros(N_g);
        for ct = 1:num_p
            gm_temp(r_grid(ct),z_grid(ct)) = gm_temp(r_grid(ct),z_grid(ct))+1;
        end
        gm2d(i,:) = reshape(flipud(gm_temp),1,[]);
    end
end

function C = Laplacian(N_grid)
    k = N_grid+1;
    K = k*k;
    C=zeros(K);
    for i=1:1:k
        for j=1:1:k
               C((i-1)*k+j,(i-1)*k+j)=-4;
            if j+1<=k
                C((i-1)*k+j,(i-1)*k+j+1)=1;
            end
            
            if j-1>=1
                C((i-1)*k+j,(i-1)*k+j-1)=1;
            end
            
            if i-1-1>=0
                C((i-1)*k+j,(i-1-1)*k+j)=1;
            end
            
            if i-1+1<=k-1
                C((i-1)*k+j,(i-1+1)*k+j)=1;
            end
        end
    end
end

function plot_overlapping_rays(l1, l2, l3, l4, threshold)
    % 光線データのセル配列と総光線数
    rays = {l1, l2, l3, l4};
    N_rays = length(rays);
    
    % 全ての光線の合計点数を計算
    total_points = 0;
    for i = 1:N_rays
        for j = 1:numel(rays{i})
            total_points = total_points + numel(rays{i}(j).x);
        end
    end
    
    % 全ての点を格納する配列を準備
    all_points = zeros(total_points, 3);
    point_idx = 1;
    
    % 各光線の点を1つの配列にまとめる
    for i = 1:N_rays
        for j = 1:numel(rays{i})
            x = rays{i}(j).x;
            y = rays{i}(j).y;
            z = rays{i}(j).z;
            num_points = numel(x);
            all_points(point_idx:point_idx+num_points-1, :) = [x', y', z'];
            point_idx = point_idx + num_points;
        end
    end
    
    % KD-treeを使って最近傍探索
    MdlKDT = KDTreeSearcher(all_points);
    idx_pairs = rangesearch(MdlKDT, all_points, threshold);
    
    % 重なり度合いをカウントする配列
    overlap_count = zeros(total_points, 1);
    
    % 各点の重複数をカウント
    for i = 1:total_points
        overlap_count(i) = numel(idx_pairs{i});
    end
    
    % プロット設定
    figure; hold on; grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    
    % % 重複数が2以上の点だけを対象にする
    % idx_2_or_more = find(overlap_count > 1);
    
    % % 重複数が2の点を青色でプロット
    % idx_2 = find(overlap_count == 2);
    % if ~isempty(idx_2)
    %     plot3(all_points(idx_2, 1), all_points(idx_2, 2), all_points(idx_2, 3), 'bo', 'DisplayName', 'Overlap 2');
    % end
    % 
    % % 重複数が3の点を緑色でプロット
    % idx_3 = find(overlap_count == 3);
    % if ~isempty(idx_3)
    %     plot3(all_points(idx_3, 1), all_points(idx_3, 2), all_points(idx_3, 3), 'go', 'DisplayName', 'Overlap 3');
    % end
    
    % 重複数が4の点を赤色でプロット
    idx_4 = find(overlap_count == 4);
    if ~isempty(idx_4)
        plot3(all_points(idx_4, 1), all_points(idx_4, 2), all_points(idx_4, 3), 'ro', 'DisplayName', 'Overlap 4');
    end
    
    legend;
end

function gm3d = LineProjection3D(l, N_grid_x, N_grid_y, N_grid_z, xmin, xmax, ymin, ymax, zmin, zmax, plot_flag)
    % l : 光線データのセル配列 {l1, l2, l3, l4}
    % N_grid_x, N_grid_y, N_grid_z : X, Y, Z 方向のグリッド数
    % xmin, xmax, ymin, ymax, zmin, zmax : それぞれの座標範囲
    % plot_flag : プロットフラグ (trueでプロットを表示)

    % 各セルに含まれる光線の総数をカウントする
    N_p = sum(cellfun(@numel, l));  % l1, l2, l3, l4 のすべての光線データを合算
    N_gx = N_grid_x + 1;
    N_gy = N_grid_y + 1;
    N_gz = N_grid_z + 1;

    Dx = (xmax - xmin) / N_grid_x;
    Dy = (ymax - ymin) / N_grid_y;
    Dz = (zmax - zmin) / N_grid_z;

    gm3d = zeros(N_p, N_gx * N_gy * N_gz);  % 出力行列の初期化

    if plot_flag
        figure;
        hold on;
        grid on;
        xlabel('X [mm]');
        ylabel('Y [mm]');
        zlabel('Z [mm]');
    end

    % 光線データのインデックスを追跡
    idx = 1;  % 全体で何番目の光線かを追跡するインデックス

    % 各セル配列に対して処理を行う
    for n = 1:length(l)
        current_l_set = l{n};  % 各セル配列の光線データを取得

        % current_l_set 内の各光線に対して処理を行う
        for i = 1:numel(current_l_set)
            current_l = current_l_set(i);

            % 各座標を取得
            x = current_l.x;
            y = current_l.y;
            z = current_l.z;

            if plot_flag
                plot3(x, y, z, '-');
            end

            % 再構成領域内の視線を抽出
            k = find(x >= xmin & x <= xmax & y >= ymin & y <= ymax & z >= zmin & z <= zmax);
            pl_x = x(k);
            pl_y = y(k);
            pl_z = z(k);

            % 各点のグリッド座標を計算
            x_grid = fix((pl_x - xmin) ./ Dx) + 1;
            y_grid = fix((pl_y - ymin) ./ Dy) + 1;
            z_grid = fix((pl_z - zmin) ./ Dz) + 1;

            % 各グリッドに含まれる点の数を数え上げる
            gm_temp = zeros(N_gx, N_gy, N_gz);  % 3次元グリッドのテンポラリ
            num_p = numel(x_grid);
            for ct = 1:num_p
                gm_temp(x_grid(ct), y_grid(ct), z_grid(ct)) = gm_temp(x_grid(ct), y_grid(ct), z_grid(ct)) + 1;
            end

            % 3次元グリッドを1次元に変換して `gm3d` に格納
            gm3d(idx, :) = reshape(flipud(gm_temp), 1, []);

            % インデックスを次の光線に進める
            idx = idx + 1;
        end
    end
end

function C = Laplacian_3D(N_grid)
    nx = N_grid + 1;      % x-dimension
    ny = 2 * N_grid + 1;  % y-dimension
    nz = N_grid + 1;      % z-dimension
    K = nx * ny * nz;     % Total number of grid points
    C = sparse(K, K);     % Initialize sparse Laplacian matrix

    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                idx = (i-1)*ny*nz + (j-1)*nz + k;  % Linear index for the current point
                
                % Set the center value (-6 for 3D Laplacian)
                C(idx, idx) = -6;
                
                % x-direction neighbors
                if i+1 <= nx
                    C(idx, (i)*ny*nz + (j-1)*nz + k) = 1;  % Forward x-direction
                end
                if i-1 >= 1
                    C(idx, (i-2)*ny*nz + (j-1)*nz + k) = 1;  % Backward x-direction
                end
                
                % y-direction neighbors
                if j+1 <= ny
                    C(idx, (i-1)*ny*nz + (j)*nz + k) = 1;  % Forward y-direction
                end
                if j-1 >= 1
                    C(idx, (i-1)*ny*nz + (j-2)*nz + k) = 1;  % Backward y-direction
                end
                
                % z-direction neighbors
                if k+1 <= nz
                    C(idx, (i-1)*ny*nz + (j-1)*nz + k + 1) = 1;  % Forward z-direction
                end
                if k-1 >= 1
                    C(idx, (i-1)*ny*nz + (j-1)*nz + k - 1) = 1;  % Backward z-direction
                end
            end
        end
    end
end
