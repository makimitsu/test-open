clear all

% ********* ALWAYS RUN THIS FUNCTION FIRST ***********
% parameters:(date,TF_shot,shot,offset_TF,offset_EF)
%[B_z,r_probe,z_probe,ch_dist,data] = get_B_z(200202,8,14,true,true);
[B_z,r_probe,z_probe,ch_dist,data] = get_B_z(200430,12,13,true,true);

% run this as well (probes at one edge is not used)
B_z = B_z(:,2:end,:);
z_probe = z_probe(2:end);
ch_dist = ch_dist(:,2:end);

% ************* PLOTTING FUNCTIONS *******************
% parameters:(B_z,ch_dist,start_time,end_time)
plot_B_z_in_time(B_z,ch_dist,350,600);

% parameters:(B_z,r_probe,z_probe,t,fitting,fill,fixed_Clayer,show_probe)
plot_psi_at_t(B_z,r_probe,z_probe,469,true,true,false,true);