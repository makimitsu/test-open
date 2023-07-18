NofCH = 24;%�`�����l����
NT = 9;%�v���g���K������
NData = 5;%����v���g���K�ł̌v���V���b�g��
nz = 1;%z�����v���_��
factor = 0.3;%�t���[���̒����{��
r_measured = zeros(NofCH/4,nz);%�x�N�g���v���b�gr���W1���data1�A2���data2
z_measured = zeros(NofCH/4,nz);%�x�N�g���v���b�gz���W1���data1�A2���data2
z_measured(:,1) = 2.1;%data1��z���W
r_measured(:,1) = [22.5, 25, 27.5, 30, 32.5, 35];%data1��r���W

fig = figure;
frames(NT) = struct('cdata', [], 'colormap', []); % �e�t���[���̉摜�f�[�^���i�[����z��

for t = 1:NT
    Ti_multi = zeros(NofCH/4,NData);
    Vz_multi = zeros(NofCH/4,NData);
    Vr_multi = zeros(NofCH/4,NData);
    for ndata = 1:NData
        load(['magprobe/mat/save_',num2str(463+2*t),'us_',num2str(ndata),'.mat'],'V','Ti')
        Ti_multi(:,ndata) = Ti(1:NofCH/4,1);
        Vz_multi(:,ndata) = V(1:NofCH/4,1);
        Vr_multi(:,ndata) = V(1:NofCH/4,2);
    end
    Ti_mean = trimmean(Ti_multi,60,2);
    Vz_mean = trimmean(Vz_multi,60,2);
    Vr_mean = trimmean(Vr_multi,60,2);
    absV_mean = zeros(NofCH/4,1);
    if t == 1
        Ti_offset = Ti_mean;
    end
    for i = 1:NofCH/4
        absV_mean(i,1) = sqrt(Vz_mean(i,1)^2 + Vr_mean(i,1)^2);
    end
    
    %���ϗ����A���ω��x���v���b�g
    %     Ticon = repmat(Ti_mean,1,2);%�������}�p���σC�I�����x
    Ticon = repmat(Ti_mean - Ti_offset,1,2);%�������}�p���σC�I�����x�㏸
    zcon = z_measured(1,1);
    s = pcolor([zcon-0.8 zcon+0.8],r_measured,Ticon);
    s.FaceColor = 'interp';
    s.FaceAlpha = 0.6;%�J���[�}�b�v�s�����x(0~1)
    s.EdgeAlpha = 0;
    colormap('jet')
    %     caxis([55 110])
    caxis([0 50])
    hold on
    plot(z_measured,r_measured,'xr','MarkerSize',5,'LineWidth',2);
    hold on
    if nz == 1
        q = quiver(z_measured,r_measured,Vz_mean*factor,Vr_mean*factor);
    end
    if nz == 2
        q = quiver(z_measured,r_measured,[V(:,1),V(:,3)]*factor,[V(:,2),V(:,4)]*factor);
    end
    q.LineWidth = 3;
    q.MaxHeadSize = 15;
    q.AutoScale = 'off';
    q.Color = 'black';
    xlim([min(z_measured,[],'all')-2.5 max(z_measured,[],'all')+2.5])
    ylim([min(r_measured,[],'all')-1 max(r_measured,[],'all')+1])
    %     title({[num2str(463+2*t),'us'];'Ion Flow [km/s]'},'Color','black','FontWeight','bold')
    title([num2str(463+2*t),'us'],'Color','black','FontWeight','bold')
    absV_mean = round(absV_mean,1);
    for j = 1:nz
        for i = 1:NofCH/4
            txt = text(z_measured(i,j)+0.9,r_measured(i,j)-0.2,num2str(absV_mean(i,j)));
            txt.FontSize = 16;
            txt.Color = 'k';
            txt.FontWeight = 'bold';
        end
    end
    xlabel('Z [cm]')
    ylabel('R [cm]')
    ax = gca;
    ax.FontSize = 12;
    grid on
    daspect([1 1 1])
    c = colorbar;
    % c.Label.String = 'Ion Temperature [eV]';
    c.Label.String = 'Increase of Ion Temperature [eV]';
    c.FontSize = 12;
    drawnow; % �`����m���Ɏ��s������
    frames(t) = getframe(fig); % �}���摜�f�[�^�Ƃ��ē���
    hold off
end

%gif�t�@�C����ۑ�
filename = 'movie_flow.gif'; % �t�@�C����
for t = 1:NT
    [A, map] = rgb2ind(frame2im(frames(t)), 256); % �摜�`���ϊ�
    if t == 1
        imwrite(A, map, filename, 'gif', 'DelayTime', 1); % �o�͌`��(FPS)��ݒ�
    else
        imwrite(A, map, filename, 'gif', 'DelayTime', 1, 'WriteMode', 'append'); % 2�t���[���ڈȍ~��"�ǋL"�̐ݒ���K�v
    end
end
