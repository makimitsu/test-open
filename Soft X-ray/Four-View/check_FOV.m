clear
N_projection = 50; % projection number
N_grid = 90; % grid number
z_plot_lim = 120; % plot
% z_plot_lim = 320; % plot areaarea

zhole1=40;zhole2=-40;                                  
zmin1=-240;zmax1=320;zmin2=-320;zmax2=240;             
rmin=55;rmax=375;

f = figure;
f.Units = 'normalized';
f.Position = [0.1,0.2,0.8,0.8];

l1 = MCPLine_up(N_projection,zhole2,false);
gm2d1 = LineProjection(l1,N_grid,zmin2,zmax2,rmin,rmax,false,true); 
l2 = MCPLine_down(N_projection,zhole2,false);
gm2d2 = LineProjection(l2,N_grid,zmin2,zmax2,rmin,rmax,false,false);
l3 = MCPLine_up(N_projection,zhole1,false);
gm2d3 = LineProjection(l3,N_grid,zmin1,zmax1,rmin,rmax,false,true);
l4 = MCPLine_down(N_projection,zhole1,false);
gm2d4 = LineProjection(l4,N_grid,zmin1,zmax1,rmin,rmax,false,false);

FOV = zeros(N_grid+1,N_grid+1,4);
FOV(:,:,1) = reshape(sum(gm2d1),N_grid+1,[]); % 縦軸r, 横軸zで左下最小
FOV(:,:,2) = reshape(sum(gm2d2),N_grid+1,[]);
FOV(:,:,3) = reshape(sum(gm2d3),N_grid+1,[]);
FOV(:,:,4) = reshape(sum(gm2d4),N_grid+1,[]);

r_space_SXR = linspace(rmin,rmax,N_grid+1);
z_space_SXR1 = linspace(zmin1,zmax1,N_grid+1);
z_space_SXR2 = linspace(zmin2,zmax2,N_grid+1);

r_range = find(60<=r_space_SXR & r_space_SXR<=330);
r_space_SXR = r_space_SXR(r_range);
z_range1 = find(-1*z_plot_lim<=z_space_SXR1 & z_space_SXR1<=z_plot_lim);
z_range2 = find(-1*z_plot_lim<=z_space_SXR2 & z_space_SXR2<=z_plot_lim);

z_space_SXR1 = z_space_SXR1(z_range1);
z_space_SXR2 = z_space_SXR2(z_range2);

[SXR_mesh_z1,SXR_mesh_r] = meshgrid(z_space_SXR1,r_space_SXR);
[SXR_mesh_z2,~] = meshgrid(z_space_SXR2,r_space_SXR);

positionList = [2,4,1,3];

for i = 1:4
    if i <= 2
        z_range = z_range2;
        SXR_mesh_z = SXR_mesh_z2;
    else
        z_range = z_range1;
        SXR_mesh_z = SXR_mesh_z1;
    end
    FOV_plot = flipud(FOV(r_range,z_range,i)); % 縦軸r, 横軸zで左上最小
    p = positionList(i);
    subplot(2,2,p);
    % [~,h] = contourf(SXR_mesh_z,SXR_mesh_r,FOV_plot,20);
    % h.LineStyle = 'none';
    pcolor(SXR_mesh_z,SXR_mesh_r,FOV_plot);
    clim([0 1]);
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

