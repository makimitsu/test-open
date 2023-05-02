function [] = plot_pcb_SXR_at_t(psi_all,rq,zq,date,shot,t,layer,area,start,interval,ifsave,SXRfilename)
% plot SXR emission on psi in rz plane
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   integer: date, date of experiment
%   integer: shot, number of shot
%   integer: t, time of interest (us)
%   boolean: layer, option for changing the contour property
%   boolean: area, option for narrowing the reconstruction area
%   integer: start, start time (us)
%   integer: interval, interval time of the framing camera (us)
%   boolean: save, option for saving the reconstruction result
%   string: SXRfilename, name of the SXR image file

    filename = strcat('SXR_',num2str(date),num2str(shot,'%03i'),'_',num2str(t),'us');
    filemat=strcat(filename,'.mat');

    if isfile(filemat)
        load(filemat)
        psi_mesh_z=zq;
        psi_mesh_r=rq;
        f = figure;
        % f.WindowState = 'maximized';
        f.Units = 'normalized';
        f.Position = [0.1,0.2,0.8,0.6];
        pos1 = [0.07,0.2,0.35,0.6];
        pos2 = [0.58,0.2,0.35,0.6];
        
        % subplot(2,1,1);
        subplot('Position',pos1);
        [~,h1] = contourf(SXR_mesh_z1,SXR_mesh_r1,EE1,40);
        h1.LineStyle = 'none';
        % caxis([0,0.5]);
        % caxis([0,0.2]);
        c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
        hold on
        [~,hp1]=contourf(psi_mesh_z,psi_mesh_r,psi,30,'-k','Fill','off');
    
        title(strcat(num2str(t),' us'));
        xlabel('z [m]');
        ylabel('r [m]');
        ax = gca;
        ax.FontSize = 18; 
        hold off
            axis image
            axis tight manual
        % subplot(2,1,2);
        subplot('Position',pos2);
    
        [~,h2] = contourf(SXR_mesh_z2,SXR_mesh_r2,EE2,40);
        h2.LineStyle = 'none';
        % caxis([0,0.15]);
        % caxis([0,0.2]);
        c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
        hold on
        [~,hp2]=contourf(psi_mesh_z,psi_mesh_r,psi,30,'-k','Fill','off');
        title(strcat(num2str(t),' us'));
        xlabel('z [m]');
        ylabel('r [m]');
        ax = gca;
        ax.FontSize = 18; 
            axis image
            axis tight manual
        hold off
    %save(filemat,'psi','t','SXR_mesh_z1','SXR_mesh_r1','EE1','SXR_mesh_z2','SXR_mesh_r2','EE2')
    
%     filename = strcat('/',num2str(t),'us.png');
%    saveas(gcf,strcat(foldername,filename));
    else
        number = (t-start)/interval+1;
        psi=psi_all(:,:,number);
        if layer
            psi_min = min(min(psi));
            psi_max = max(max(psi));
            contour_layer = linspace(psi_min,psi_max,25);
        else
            max_psi_color = 0.0400;
            min_psi_color = -0.0200;
            layer_resolution = 0.0002;
            contour_layer =  min_psi_color:layer_resolution:max_psi_color;
        end
        if date == 210309
            if shot <= 15
                number = (t-430)/10;
            else
                number = (t-445)/5;
            end
        elseif date == 210310
            number = (t-445)/5;
        elseif date == 210319
            if shot == 4
                number = (t-445)/5;
            elseif shot <= 6
                number = (t-440)/10;
            else
                number = (t-415)/10; %ここいじった
            end
        end
        
        zmin1=-240;zmax1=320;zmin2=-320;zmax2=240;
        zhole1=40;zhole2=-40;
        
        [VectorImage1,VectorImage2] = get_SXRImage(date,shot,number,SXRfilename);
        % [EE1,EE2] = get_SXR(VectorImage1,VectorImage2,false);
        EE1 = clc_SXR(VectorImage1,zmin1,zmax1,zhole1,false);
        EE2 = clc_SXR(VectorImage2,zmin2,zmax2,zhole2,false);
        
        % % pathname = '/Users/shinjirotakeda/Dropbox/ReconstructionResults/';
        % pathname = '/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/ReconstructionResults/';
        % if area
        %     foldername = strcat(pathname,num2str(date),'/shot',num2str(shot));
        % else
        %     foldername = strcat(pathname,num2str(date),'/shot',num2str(shot),'_wide');
        % end
        % if exist(foldername,'dir') == 0
        %     mkdir(foldername);
        % end
        
        r_space_SXR = linspace(0.055,0.375,size(EE1,1));
        r_range = 1:numel(r_space_SXR);
        z_space_SXR1 = linspace(zmin1/1000,zmax1/1000,size(EE1,2));
        z_space_SXR2 = linspace(zmin2/1000,zmax2/1000,size(EE2,2));
        z_range1 = find(-0.12<=z_space_SXR1 & z_space_SXR1<=0.20);
        z_range2 = find(-0.20<=z_space_SXR2 & z_space_SXR2<=0.12);
        z_space_SXR1 = z_space_SXR1(z_range1);
        z_space_SXR2 = z_space_SXR2(z_range2);
        if area
            r_range = find(rq(end,1)<=r_space_SXR & r_space_SXR<=rq(1,1));
            r_space_SXR = r_space_SXR(r_range);
        % else
        %     z_space_SXR1 = linspace(zmin1/1000,zmax1/1000,size(EE1,2));
        %     z_space_SXR2 = linspace(zmin2/1000,zmax2/1000,size(EE1,2));
        end
        % figure;imagesc(EE1);
        EE1 = EE1(r_range,z_range1);
        EE2 = EE2(r_range,z_range2);
        % figure;imagesc(EE1);
        % EE1 = flipud(EE1);
        % EE2 = flipud(EE2);
        
        f = figure;
        % f.WindowState = 'maximized';
        f.Units = 'normalized';
        f.Position = [0.1,0.2,0.8,0.6];
        pos1 = [0.07,0.2,0.35,0.6];
        pos2 = [0.58,0.2,0.35,0.6];
        
        % subplot(2,1,1);
        subplot('Position',pos1);
        [SXR_mesh_z1,SXR_mesh_r1] = meshgrid(z_space_SXR1,r_space_SXR);
        [~,h1] = contourf(SXR_mesh_z1,SXR_mesh_r1,EE1,40);
        h1.LineStyle = 'none';
        % caxis([0,0.5]);
        % caxis([0,0.2]);
        c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
        hold on
        [~,hp1]=contourf(zq,rq,psi,30,'-k','Fill','off');
        if layer
            hp1.LineWidth = 1;
        end
        title(strcat(num2str(t),' us'));
        xlabel('z [m]');
        ylabel('r [m]');
        ax = gca;
        ax.FontSize = 18; 
        hold off
            axis image
            axis tight manual
            xlim([min(zq,[],'all') max(zq,[],'all')])
            ylim([min(rq,[],'all') max(rq,[],'all')])


        % subplot(2,1,2);
        subplot('Position',pos2);
        [SXR_mesh_z2,SXR_mesh_r2] = meshgrid(z_space_SXR2,r_space_SXR);
        [~,h2] = contourf(SXR_mesh_z2,SXR_mesh_r2,EE2,40);
        h2.LineStyle = 'none';
        % caxis([0,0.15]);
        % caxis([0,0.2]);
        c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
        hold on
        [~,hp2]=contourf(zq,rq,psi,30,'-k','Fill','off');
        if layer
            hp2.LineWidth = 1;
        end
        title(strcat(num2str(t),' us'));
        xlabel('z [m]');
        ylabel('r [m]');
        ax = gca;
        ax.FontSize = 18; 
            axis image
            axis tight manual
            xlim([min(zq,[],'all') max(zq,[],'all')])
            ylim([min(rq,[],'all') max(rq,[],'all')])
        hold off
    end
if ifsave
    % pathname = '/Users/shinjirotakeda/Dropbox/ReconstructionResults/';
    %cd 'I:\makimitsu\211223';
    cd 'C:\Users\Moe Akimitsu\Desktop\dronkaiseki'
%     if area
%         foldername = strcat(pathname,num2str(date),'/shot',num2str(shot));
%     else
%         foldername = strcat(pathname,num2str(date),'/shot',num2str(shot),'_wide');
%     end
%     if exist(foldername,'dir') == 0
%         mkdir(foldername);
%     end
%    filename = strcat('/shot',num2str(shot),'_',num2str(t),'us.png');
     save(strcat(filename,'.mat'),'psi','zq','rq','SXR_mesh_z1','SXR_mesh_r1','EE1','SXR_mesh_z2','SXR_mesh_r2','EE2')
     saveas(gcf,strcat(filename,'.png'))
     close
end

end