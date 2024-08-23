function [] = get_parameters_new(N_projection,N_grid,filepath)

% filepath = '/Users/shinjirotakeda/Documents/GitHub/test-open/Soft X-ray/Four-View_Simulation/parameters.mat';

plotFlagProjection = false;
plotFlagLines = false;

% 視線の分布、重み行列の作成
zhole1=40;zhole2=-40;                                  
zmin1=-100;zmax1=180;zmin2=-180;zmax2=100;             
rmin=70;rmax=330;
range = [zmin1,zmax1,zmin2,zmax2,rmin,rmax];            
l1 = MCPLine(N_projection,zhole2,plotFlagLines,true);
gm2d1 = LineProjection(l1,N_grid,zmin2,zmax2,rmin,rmax,plotFlagProjection,true); 
l2 = MCPLine(N_projection,zhole2,plotFlagLines,false);
gm2d2 = LineProjection(l2,N_grid,zmin2,zmax2,rmin,rmax,plotFlagProjection,false);
l3 = MCPLine(N_projection,zhole1,plotFlagLines,true);
gm2d3 = LineProjection(l3,N_grid,zmin1,zmax1,rmin,rmax,plotFlagProjection,true);
l4 = MCPLine(N_projection,zhole1,plotFlagLines,false);
gm2d4 = LineProjection(l4,N_grid,zmin1,zmax1,rmin,rmax,plotFlagProjection,false);

% ラプラシアン行列の計算と特異値分解
C = Laplacian(N_grid);
[U1,S1,V1]=svd(gm2d1*(C^(-1)),'econ');
[U2,S2,V2]=svd(gm2d2*(C^(-1)),'econ');
[U3,S3,V3]=svd(gm2d3*(C^(-1)),'econ');
[U4,S4,V4]=svd(gm2d4*(C^(-1)),'econ');
v1=(C^(-1)*V1);
v2=(C^(-1)*V2);
v3=(C^(-1)*V3);
v4=(C^(-1)*V4);

[M,K] = size(gm2d1);
if K>M
    v1 = v1(:,1:M);
    v2 = v2(:,1:M);
    v3 = v3(:,1:M);
    v4 = v4(:,1:M);
end
s1 = (diag(S1)).';
s2 = (diag(S2)).';
s3 = (diag(S3)).';
s4 = (diag(S4)).';
if M>K
    s1 = [s1 zeros(1,M-K)];
    s2 = [s2 zeros(1,M-K)];
    s3 = [s3 zeros(1,M-K)];
    s4 = [s4 zeros(1,M-K)];
end

save(filepath,'gm2d1','gm2d2','gm2d3','gm2d4', ...
    'U1','U2','U3','U4','s1','s2','s3','s4', ...
    'v1','v2','v3','v4','M','K','range','N_projection','N_grid');

end


function l = MCPLine(N_projection,Z_hole,plot_flag,up_flag)
d_hole = 24.4; % distance between the hole and the MCP
r_mcp=10;  % radius of the MCP plate

Y_hole = 413.24+12;

if up_flag
    X_hole = 208.13;
else
    X_hole = 144.13;
end

Y_initial=Y_hole+d_hole;  X_initial=X_hole-r_mcp;  Z_initial=Z_hole+r_mcp;    %CCD position
X_end=X_hole+r_mcp;      Z_end=Z_hole-r_mcp;

Nh=N_projection-1; 
Dhx=(X_end-X_initial)/Nh;
Dhz=(Z_end-Z_initial)/Nh;

x_mcp=X_initial:Dhx:X_end;
z_mcp=Z_initial:Dhz:Z_end;
[X,Z] = meshgrid(x_mcp,z_mcp); %Xは縦に一様で右に向かって増える、Zは横に一様で上に向かって増える
X=flipud(rot90(X));Z=flipud(rot90(Z)); %Xは横に一様で下に向かって増える、Zは縦に一様で左に向かって増える
Y=repelem(Y_initial,N_projection,N_projection);

r_center = 55;
r_device = 375;

ll_y = repmat(reshape(Y_initial:-10:-400,1,1,[]),N_projection,N_projection,1);
if plot_flag
    f1=figure;
    if up_flag
        f1.Name = "Upper MCP Line - Lightlines";
    else
        f1.Name = "Lower MCP Line - Lightlines";
    end
    f2=figure;
    if up_flag
        f2.Name = "Upper MCP Line - 視線の行列";
    else
        f2.Name = "Lower MCP Line - 視線の行列";
    end
end

ll_x = (ll_y-Y_hole).*(X-X_hole)./(Y-Y_hole)+X_hole;
ll_z = (ll_y-Y_hole).*(Z-Z_hole)./(Y-Y_hole)+Z_hole;

nLine = numel(reshape(Y_initial:-10:-400,1,1,[]));

r = sqrt(ll_y.^2+ll_x.^2);
A = find(r<=r_center | r>= r_device);
ll_x(A) = 0;
ll_y(A) = 0;
ll_z(A) = 0;

if plot_flag
    figure(f1);plot3(X,Y,Z,'*',reshape(ll_x,[],nLine),reshape(ll_y,[],nLine),reshape(ll_z,[],nLine),'.');
end


%視線の行列のうち円の内部に含まれるものだけをベクトル化
L = N_projection/2;
k = find_circle(L); %左上から下方向にカウントする線型インデックス→X最小、Z最大からZをまず減らし、その後Xを増やしている
[row,col] = ind2sub([N_projection N_projection],k);

Np = numel(k);
l.x = zeros(Np,nLine);
l.y = zeros(Np,nLine);
l.z = zeros(Np,nLine);
for i = 1:Np
    l.x(i,:) = ll_x(row(i),col(i),:);
    l.y(i,:) = ll_y(row(i),col(i),:);
    l.z(i,:) = ll_z(row(i),col(i),:);
end

if plot_flag
    figure(f2);plot3(reshape(l.x,[],nLine),reshape(l.y,[],nLine),reshape(l.z,[],nLine),'.');
end
end



function gm2d = LineProjection(l,N_grid,zmin,zmax,rmin,rmax,plot_flag,up_flag)
% rmin=70;rmax=330;
% rmin=55;rmax=375;
% rmin=70;rmax=280;

% N_p = numel(l);
N_p = size(l.z,1);
N_g = N_grid;
DR=(rmax-rmin)/N_grid;
DZ=(zmax-zmin)/N_grid;
gm2d = zeros(N_p,N_g^2);

if plot_flag
    f1=figure;
    % f2=figure;
    if up_flag
        f1.Name = "Upper Line Projection";
        % f2.Name = "Upper Line Projection";
    else
        f1.Name = "Lower Line Projection";
        % f2.Name = "Lower Line Projection";
    end
end

l_z = l.z;
l_r = sqrt((l.x).^2+(l.y).^2);
nLine = size(l_z,2);
if plot_flag
    figure(f1);plot(reshape(l_z,[],nLine),reshape(l_r,[],nLine),'.');
    hold on;grid on;
    xlabel('Z [mm]');ylabel('R [mm]');
end

l_z = reshape(l_z,N_p,nLine); %線型インデックスと同じ順番
z_grid = ceil((l_z-zmin)./DZ);
l_r = reshape(l_r,N_p,nLine);
r_grid = ceil((l_r-rmin)./DR);
l_grid_rz = cat(3,r_grid,z_grid);

FOV = zeros(N_g);

for i = 1:N_p
    A = squeeze(l_grid_rz(i,:,:));
    [uniqueRows, ~, rowIndices] = unique(A, 'rows', 'stable');%左がr,右がz
    counts = histcounts(rowIndices, 'BinMethod', 'integers', 'BinLimits', [1, size(uniqueRows, 1)]);
    gm_temp = zeros(N_g); %左上が最小
    I = find(min(uniqueRows,[],2)>=1&max(uniqueRows,[],2)<=N_g);
    uniqueRows = uniqueRows(I,:);
    counts = counts(I);
    gm_temp(sub2ind(size(gm_temp),uniqueRows(:,1),uniqueRows(:,2))) = counts.';%縦がr、横がzの配列（左上最小）
    gm_temp = flipud(gm_temp);%縦がr、横がzの配列（左下最小）
    % gm_temp = rot90(gm_temp); %縦がr、横がzの配列（左下最小） なんか向き違うかも　縦がzになってる
    FOV = FOV + gm_temp;
    gm2d(i,:) = reshape(gm_temp,1,[]); %RZ最小からまずzを増やし、その後rを増やす
end

if plot_flag
    figure;imagesc(FOV);clim([0 1]);
end

end

function C = Laplacian(N_grid)
k=N_grid;
K=k*k;
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
