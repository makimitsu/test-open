clear all
close all
% this script reads:
% processed Ti data calculated using Photron_passive and magnetic field data
% this script calculates and saves:
% thermal speed and gyroradius and s parameter

log_cell = readcell('/Users/yunhancai/Google Drive/Data/log2022-2023.xlsx','Sheet','TOTAL');
a039_shotlist_target=2746:2920;%【input】dtacqの保存番号
a039_shot = cell2mat(log_cell(2:end,3));
[a039_list,idx] = intersect(a039_shot,a039_shotlist_target,'stable');
a040_list = cell2mat(log_cell(idx+1,5));
a039_tf_list = cell2mat(log_cell(idx+1,4));
a040_tf_list = cell2mat(log_cell(idx+1,6));
gas_list = string(log_cell(idx+1,25));
shotlist = [a039_list,a040_list];
tfshotlist = [a039_tf_list,a040_tf_list];
i_cai_list = cell2mat(log_cell(idx+1,29));
i_cai_tf_list = zeros(size(idx));
date_list = cell2mat(log_cell(idx+1,2));%【input】計測日
n_data = numel(a039_list);

z_measurement = 0; % measurement at z=0
Z = 1; % electric charge
kb = 1.380649 * 10^-23; % boltzmann const
u = 1.6605402 * 10^-27; % atomic mass
qc = 1.602176634 * 10^-19; % coulumbd charge
c = 299792458; % speed of light
ev2kelvin = 11604.525;

for index = 1:n_data
shot_a039_num = a039_list(index);
filename_doppler = strcat('/Users/yunhancai/Google Drive/Data/Doppler/Photron/processed/Ti_a039_',num2str(shot_a039_num),'.mat');
filename_magnetic = strcat('/Users/yunhancai/Google Drive/Data/pcb_probe/pre_processed/a039_',num2str(shot_a039_num),'.mat');
s_kinetic_path = '/Users/yunhancai/Google Drive/Data/Doppler/Photron/ti_related_calculation/';
if ~exist(filename_doppler,"file") || ~exist(filename_magnetic,"file")
    disp('doppler or magnetic data not found.')
    continue
end
load(filename_doppler),load(filename_magnetic);
gas = gas_list(index);
if gas == "H"
    mu_gas = 1.00784; % m_gas/m_proton
elseif gas == "Ar"
    mu_gas = 39.948;
else
    disp("Gas not specified.");
    continue;
end

[~,z_index] = min(abs(grid2D.zprobepcb-z_measurement));
B_total = (data2D.Bz(:,z_index,:).^2+data2D.Br(:,z_index,:).^2+data2D.Bt(:,z_index,:).^2).^0.5;
B_total = squeeze(B_total);
B_data_index = round(time2)-data2D.trange(1)+1;
[t_mesh,r_mesh] = meshgrid(data2D.trange(B_data_index),grid2D.rprobepcb); % new mesh
B_total_interp = griddata(data2D.trange,grid2D.rprobepcb,B_total,t_mesh,r_mesh); % interpolate B to new mesh
Ti_local_interp = smoothdata(griddata(time2,yy,Ti_local2,t_mesh,r_mesh) * ev2kelvin); % interpolate Ti to new mesh
Ti_local_interp_max = smoothdata(griddata(time2,yy,Ti_local_max2,t_mesh,r_mesh) * ev2kelvin);
Ti_local_interp_min = smoothdata(griddata(time2,yy,Ti_local_min2,t_mesh,r_mesh) * ev2kelvin);

% calculates minor radius
minor_radius = zeros(size(time2));
separatrix_index = zeros(size(time2));
for i = 1:numel(time2)
    pp = spline(grid2D.rprobepcb,data2D.psi(:,z_index,B_data_index(i)));
    pp_zeros = fnzeros(pp,[0.05 0.4]);
    if numel(pp_zeros(1,:)) == 2
        minor_radius(i) = (pp_zeros(1,2)-pp_zeros(1,1))/2;
        [~,separatrix_index(i)] = min(abs(grid2D.rprobepcb-pp_zeros(1,2)));
    else
        minor_radius(i) = NaN;
    end
end

% calculates v_Ti and iongyro profile
v_Ti = (kb*Ti_local_interp/(u*mu_gas)).^0.5; % v_ti = (2*kb*Ti/m_ti)^0.5
v_Ti_max = (kb*Ti_local_interp_max/(u*mu_gas)).^0.5; % v_ti = (2*kb*Ti/m_ti)^0.5
v_Ti_min = (kb*Ti_local_interp_min/(u*mu_gas)).^0.5; % v_ti = (2*kb*Ti/m_ti)^0.5
iongyro = mu_gas*u * v_Ti./(qc*Z*B_total_interp);iongyro(iongyro==0) = NaN;
iongyro_max = mu_gas*u * v_Ti_max./(qc*Z*B_total_interp);iongyro(iongyro==0) = NaN;
iongyro_min = mu_gas*u * v_Ti_min./(qc*Z*B_total_interp);iongyro(iongyro==0) = NaN;

% calculates s by integration: s =
% integrate(r*dr/(R_separatrix*iongyro),R_OH to R_separatrix),belova2003
dr = abs(grid2D.rprobepcb(1)-grid2D.rprobepcb(2));
s_integrant = repmat(grid2D.rprobepcb',1,size(iongyro,2))./(iongyro)*dr;
s_integrant_max = repmat(grid2D.rprobepcb',1,size(iongyro_min,2))./(iongyro_min)*dr;
s_integrant_min = repmat(grid2D.rprobepcb',1,size(iongyro_max,2))./(iongyro_max)*dr;
s_integrant(isnan(s_integrant)) = 0;s_integrant(isinf(s_integrant)) = 0;
s_integrant_max(isnan(s_integrant_max)) = 0;s_integrant_max(isinf(s_integrant_max)) = 0;
s_integrant_min(isnan(s_integrant_min)) = 0;s_integrant_min(isinf(s_integrant_min)) = 0;
s = zeros(size(time2));s_max = zeros(size(time2));s_min = zeros(size(time2));
for i = 1:numel(time2)
    if separatrix_index(i) == 0 
        continue
    end
    r_separatrix = grid2D.rprobepcb(separatrix_index(i));
    s(i) = sum(s_integrant(1:separatrix_index(i),i)/r_separatrix,1);
    s_max(i) = sum(s_integrant_max(1:separatrix_index(i),i)/r_separatrix,1);
    s_min(i) = sum(s_integrant_min(1:separatrix_index(i),i)/r_separatrix,1);
end

% calculates iongyro and s by averaging over entire r
iongyro_avg = mean(iongyro,1,"omitnan");
iongyro_avg_max = mean(iongyro_max,1,"omitnan");
iongyro_avg_min = mean(iongyro_min,1,"omitnan");
s_avg = 2*minor_radius./iongyro_avg;
s_avg_max = 2*minor_radius./iongyro_avg_min;
s_avg_min = 2*minor_radius./iongyro_avg_max;

% calculates v_Ti, iongyro and s using averaged line-integrated Ti
v_Ti_int = (kb*mean(Ti_line_integrated*ev2kelvin,1,"omitnan")/(u*mu_gas)).^0.5;
v_Ti_int_max = (kb*mean(Ti_line_integrated_max*ev2kelvin,1,"omitnan")/(u*mu_gas)).^0.5;
v_Ti_int_min = (kb*mean(Ti_line_integrated_min*ev2kelvin,1,"omitnan")/(u*mu_gas)).^0.5;
iongyro_int = zeros(size(time2));
iongyro_int_max = zeros(size(time2));
iongyro_int_min = zeros(size(time2));
for i = 1:numel(time2)
    iongyro_int(i) = mu_gas*u * v_Ti_int(i)/(qc*Z*mean(B_total_interp(1:separatrix_index(i),i),1));
    iongyro_int_max(i) = mu_gas*u * v_Ti_int_max(i)/(qc*Z*mean(B_total_interp(1:separatrix_index(i),i),1));
    iongyro_int_min(i) = mu_gas*u * v_Ti_int_min(i)/(qc*Z*mean(B_total_interp(1:separatrix_index(i),i),1));
end
s_int = 2*minor_radius./iongyro_int;
s_int_max = 2*minor_radius./iongyro_int_min;
s_int_min = 2*minor_radius./iongyro_int_max;

% h = figure('Position', [0 0 300 300],'visible','off');
% hold on;
% errorbar(time2,s,s-s_min,s_max-s,'k');
% errorbar(time2,s_avg,s_avg-s_avg_min,s_avg_max-s_avg,'b');
% errorbar(time2,s_int,s_int-s_int_min,s_int_max-s_int,'r');
% legend("integrated","averaged","line-integrated");xlabel("Time(s)");ylabel("s")
% hold off

h = figure('Position', [0 0 300 300],'visible','off');
idxs = logical(~isnan(s_int).*~isnan(s_int_min).*~isnan(s_int_min));
errorbar(time2(idxs),s_int(idxs),s_int(idxs)-s_int_min(idxs),s_int_max(idxs)-s_int(idxs),'k');
xlim([time2(1) time2(end)]);xlabel("Time(s)");ylabel("s");

Ti_related =struct(...
    's',s,'s_max',s_max,'s_min',s_min,...
    's_avg',s_avg,'s_avg_max',s_avg_max,'s_avg_min',s_avg_min,...
    's_int',s_int,'s_int_max',s_int_max,'s_int_min',s_int_min,...
    'iongyro',iongyro,'iongyro_max',iongyro_max,'iongyro_min',iongyro_min,...
    'iongyro_int',iongyro_int,'iongyro_int_max',iongyro_int_max,'iongyro_int_min',iongyro_int_min,...
    'v_Ti',v_Ti,'v_Ti_min',v_Ti_min,'v_Ti_max',v_Ti_max,...
    'v_Ti_int',v_Ti_int,'v_Ti_int_max',v_Ti_int_max,'v_Ti_int_min',v_Ti_int_min,...
    'Ti',Ti_local_interp/ev2kelvin,'Ti_max',Ti_local_interp_max/ev2kelvin,'Ti_min',Ti_local_interp_min/ev2kelvin,...
    'Ti_line_integrated',Ti_line_integrated,'Ti_line_integrated_max',Ti_line_integrated_max,'Ti_line_integrated_min',Ti_line_integrated_min,...
    'B_total',B_total_interp,...
    'minor_radius',minor_radius,...
    'z_measurement',z_measurement,...
    't_mesh',t_mesh,...
    'r_mesh',r_mesh,...
    'time',time2,...
    'r',grid2D.rprobepcb);

saveas(h,strcat(s_kinetic_path,'s_',num2str(shot_a039_num),'.png'));
save([s_kinetic_path,'a039_',num2str(shot_a039_num),'.mat'],"Ti_related");

close(h);
disp(['Finished ',num2str(index),'/',num2str(n_data)])
end