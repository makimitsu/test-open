function [] = plot_ionflow(V_i,absV,T_i,date,ICCD,factor,mpoints,overlay_plot,save_fig)
%V_i/absV/T_i/������/ICCD�ϐ�/���T�C�Y(���l:0.05�Ȃ�)/�v���_�z��/���C�ʂɏd�˂�/fig��ۑ�

time = round(ICCD.trg+ICCD.exp_w/2);%�v������
if overlay_plot%���C�ʂƏd�˂�
    unit = 1E-2;%�P�ʂ�cm����m�ɕϊ�
else
    unit = 1;%�P�ʕϊ����Ȃ�
    figure('Position',[600 150 300 600])
end
mpoints.z = mpoints.z * unit;
mpoints.r = mpoints.r * unit;
factor = factor * unit;

T_icon = repmat(T_i(:,1),1,5);%�������}�p���σC�I�����x
zcon = mpoints.z(1,1);
s = pcolor([zcon-2*unit zcon-1*unit zcon zcon+1*unit zcon+2*unit],mpoints.r,T_icon);
s.FaceColor = 'interp';
s.EdgeAlpha = 0;
colormap('jet')
colorbar
hold on
plot(mpoints.z,mpoints.r,'xr');%�h�b�v���[�v���[�u�v���_��\��
hold on
if mpoints.n_z == 1
    q = quiver(mpoints.z,mpoints.r,V_i(:,1)*factor,V_i(:,2)*factor);
end
if mpoints.n_z == 2
    q = quiver(mpoints.z,mpoints.r,[V_i(:,1),V_i(:,3)]*factor,[V_i(:,2),V_i(:,4)]*factor);
end
q.LineWidth = 2;
q.MaxHeadSize = 10;
q.AutoScale = 'off';
q.Color = 'k';
absV = round(absV,1);
for j = 1:mpoints.n_z
    for i = 1:mpoints.n_r
        txt = text(mpoints.z(i,j)+0.5*unit,mpoints.r(i,j)-0.2*unit,num2str(absV(i,j)));
        txt.FontSize = 18;
        txt.Color = 'k';
        txt.FontWeight = 'bold';
    end
end
if overlay_plot
    s.FaceAlpha = 0.3;%�����x(0�Ŋ��S�ɓ���)
    xlim([-0.17 0.17])
    ylim([0.06 0.33])
    ax = gca;
    ax.FontSize = 14;
else
    xlim([min(mpoints.z,[],'all')-2.5 max(mpoints.z,[],'all')+2.5])
    ylim([min(mpoints.r,[],'all')-2.5 max(mpoints.r,[],'all')+2.5])
    title([num2str(time) '[us] - shot' num2str(ICCD.shot) newline 'Ion Flow [km/s]'],'Color','black','FontWeight','bold')
    xlabel('Z [cm]')
    ylabel('R [cm]')
    ax = gca;
    ax.FontSize = 14;
    grid on
    daspect([1 1 1])
    if save_fig
        if not(exist(['ionflow_fig/',num2str(date)],'dir'))
            mkdir(sprintf("ionflow_fig/%s", num2str(date)));
        end
        saveas(gcf,['ionflow_fig/',num2str(date),'/',num2str(time),'us_shot',num2str(ICCD.shot),'.png'])
        hold off
        close
    else
        hold off
    end
end
end