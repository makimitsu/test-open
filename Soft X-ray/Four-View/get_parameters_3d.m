function [] = get_parameters_3d(N_projection,N_grid3d,filepath)%UNTITLED Summary of this function goes here

% 視線の分布、重み行列の作成
zhole1=40;zhole2=-40;                                  
% zmin1=-240;zmax1=320;zmin2=-320;zmax2=240;             
% rmin=55;rmax=375;
% range = [zmin1,zmax1,zmin2,zmax2,rmin,rmax];            
l1 = MCPLine_up(N_projection,zhole2,false);
l2 = MCPLine_down(N_projection,zhole2,false);
l3 = MCPLine_up(N_projection,zhole1,false);
l4 = MCPLine_down(N_projection,zhole1,false);

% plot_overlapping_rays(l1,l2,l3,l4,5);
%l = {l1,l2,l3,l4};
l = {l1,l2,l4};
rmin3d = 0; rmax3d = 375; dmin3d = -375; dmax3d = 375; zmin3d = -300; zmax3d = 300;

gm3d = LineProjection3D(l, N_grid3d, zmin3d, zmax3d, dmin3d, dmax3d, rmin3d, rmax3d, true);
C3d = Laplacian_3D(N_grid3d);invC3d = C3d^(-1);
[U3d,S3d,V3d] = svd(gm3d*invC3d, 'econ');
v3d = (invC3d*V3d);
[M3d,K3d] = size(gm3d);
if K3d>M3d
    v3d = v3d(:,1:M3d);
end
s3d = (diag(S3d)).';
if M3d>K3d
    s3d = [s3d zeros(1,M3d-K3d)];
end
save(filepath,'N_projection','N_grid3d','gm3d','U3d','s3d','v3d','M3d','K3d','-v7.3');

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
        ll(i,j).y=Y(j):-1:-400; 
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
        ll(i,j).y=Y(j):-1:-400; %大から小 視線の画素数をコントロール
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

function gm3d = LineProjection3D(l, N_grid3d, zmin, zmax, dmin, dmax, rmin, rmax, plot_flag)
    % l : 光線データのセル配列 {l1, l2, l3, l4}
    % N_grid_x, N_grid_y, N_grid_z : X, Y, Z 方向のグリッド数
    % xmin, xmax, ymin, ymax, zmin, zmax : それぞれの座標範囲
    % plot_flag : プロットフラグ (trueでプロットを表示)

    % 各セルに含まれる光線の総数をカウントする
    N_p = sum(cellfun(@numel, l));  % l1, l2, l3, l4 のすべての光線データを合算
    N_gz = N_grid3d + 1;
    N_gd = N_grid3d + 1;
    N_gr = N_grid3d + 1;

    Dz = (zmax - zmin) / N_grid3d;
    Dd = (dmax - dmin) / N_grid3d;
    Dr = (rmax - rmin) / N_grid3d;

    gm3d = zeros(N_p, N_gz * N_gd * N_gr);  % 出力行列の初期化

    if plot_flag
        figure;
        hold on;
        grid on;
        xlabel('Z [mm]');
        ylabel('d [mm]');
        zlabel('r [mm]');
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
            z = current_l.z;
            d = current_l.y;
            r = current_l.x;

            if plot_flag
                plot3(z, d, r, '-');
            end

            % 再構成領域内の視線を抽出
            k = find(z >= zmin & z <= zmax & d >= dmin & d <= dmax & r >= rmin & r <= rmax);
            pl_z = z(k);
            pl_d = d(k);
            pl_r = r(k);

            % 各点のグリッド座標を計算
            z_grid = fix((pl_z - zmin) ./ Dz) + 1;
            d_grid = fix((pl_d - dmin) ./ Dd) + 1;
            r_grid = fix((pl_r - rmin) ./ Dr) + 1;

            % 各グリッドに含まれる点の数を数え上げる
            gm_temp = zeros(N_gz, N_gd, N_gr);  % 3次元グリッドのテンポラリ
            num_p = numel(z_grid);
            for ct = 1:num_p
                gm_temp(r_grid(ct), z_grid(ct), d_grid(ct)) = gm_temp(r_grid(ct), z_grid(ct), d_grid(ct)) + 1;
            end
            % 3次元グリッドを1次元に変換して `gm3d` に格納
            gm3d(idx, :) = reshape(gm_temp, 1, []);

            % インデックスを次の光線に進める
            idx = idx + 1;
        end
    end
end

function C = Laplacian_3D(N_grid)
    nz = N_grid + 1;      % z-dimension
    nd = N_grid + 1;      % d-dimension
    nr = N_grid + 1;      % r-dimension
    K = nz * nd * nr;     % Total number of grid points
    
    % Estimate the total number of non-zero elements (7 per grid point)
    num_nonzeros = 7 * K;
    
    % Preallocate arrays for sparse matrix construction
    row_idx = zeros(1, num_nonzeros);
    col_idx = zeros(1, num_nonzeros);
    values = zeros(1, num_nonzeros);
    
    count = 0;
    
    for i = 1:nz
        for j = 1:nd
            for k = 1:nr
                idx = (i-1)*nd*nr + (j-1)*nr + k;  % Linear index for the current point
                
                % Center value (-6 for 3D Laplacian)
                count = count + 1;
                row_idx(count) = idx;
                col_idx(count) = idx;
                values(count) = -6;
                
                % z-direction neighbors
                if i+1 <= nz
                    count = count + 1;
                    row_idx(count) = idx;
                    col_idx(count) = (i)*nd*nr + (j-1)*nr + k;
                    values(count) = 1;
                end
                if i-1 >= 1
                    count = count + 1;
                    row_idx(count) = idx;
                    col_idx(count) = (i-2)*nd*nr + (j-1)*nr + k;
                    values(count) = 1;
                end
                
                % d-direction neighbors
                if j+1 <= nd
                    count = count + 1;
                    row_idx(count) = idx;
                    col_idx(count) = (i-1)*nd*nr + (j)*nr + k;
                    values(count) = 1;
                end
                if j-1 >= 1
                    count = count + 1;
                    row_idx(count) = idx;
                    col_idx(count) = (i-1)*nd*nr + (j-2)*nr + k;
                    values(count) = 1;
                end
                
                % r-direction neighbors
                if k+1 <= nr
                    count = count + 1;
                    row_idx(count) = idx;
                    col_idx(count) = (i-1)*nd*nr + (j-1)*nr + k + 1;
                    values(count) = 1;
                end
                if k-1 >= 1
                    count = count + 1;
                    row_idx(count) = idx;
                    col_idx(count) = (i-1)*nd*nr + (j-1)*nr + k - 1;
                    values(count) = 1;
                end
            end
        end
    end
    
    % Create sparse matrix all at once
    row_idx = row_idx(1:count);  % Trim the unused part of the preallocated arrays
    col_idx = col_idx(1:count);
    values = values(1:count);
    C = sparse(row_idx, col_idx, values, K, K);
end
