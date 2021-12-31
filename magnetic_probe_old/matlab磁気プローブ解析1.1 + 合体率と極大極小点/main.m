clear all

% ********* ALWAYS RUN THIS FUNCTION FIRST ***********
[B_z,r_probe,z_probe,ch_dist] = get_B_z('190315010.mrd','190315001.mrd');

% ************* PLOTTING FUNCTIONS *******************

%plot_B_z_in_time(B_z,ch_dist,400,600);

%plot_B_z_in_rz(B_z,r_probe,z_probe);

%plot_psi(B_z,r_probe,z_probe);

plot_fitrate(B_z,r_probe);

% merging animated
%{
[psi_mesh_z,psi_mesh_r] = meshgrid(linspace(z_probe(25),z_probe(1),25),...
                                    linspace(r_probe(7),r_probe(1),25));
figure
for i = 450:2:550
    psi = get_psi(B_z,r_probe,i);
    psi_interp = griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r);
    contourf(psi_interp,30);
    pause(0.05);
end
%}
