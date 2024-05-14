clear all
close all
addpath '/Users/rsomeya/Documents/lab/matlab/common';
run define_path.m

%--�yInput�z----
date = 240313;%�Z��������
calib_CHlist = [49];%�Z��CH�ԍ����X�g
plot_cont = true;%�������}��`��
cal_CH1 = true;%1st�pCH�ʒu����
cal_1st = true;%1st�s�[�N����
cal_CH2 = true;%2nd�pCH�ʒu����
cal_2nd = true;%2nd�s�[�N����
cal_CH3 = true;%2nd�pCH�ʒu����
cal_3rd = true;%2nd�s�[�N����
save_fit = false;%�t�B�b�e�B���O��png�ŕۑ�
save_cal = false;%�Z�����ʂ�mat�ŕۑ�
N_CH = 1;%�Z��CH����(��{�I�ɂ�1)
width = 6;%�`�����l���؂��蕝
th_ratio = 0.8;
l_mov = 15;

filename = [pathname.TE,'/Andor/20',num2str(date),'/CVI_filter_49CH_Neon_lamp.asc'];

for m = 1:size(calib_CHlist,2)
    calib_CH = calib_CHlist(m);%�Z��CH�ԍ�
    cal_filename = [pathname.TE,'/Andor/20',num2str(date),'/CVI_filter_49CH_Neon_lamp.asc'];%ICCD�t�@�C����
    lambda_1st = 529.8189;%��1�s�[�N�g��
    lambda_2nd = 528.0086;%��2�s�[�N�g��
    lambda_3rd = 530.4758;%��3�s�[�N�g��

    Min_CH_CH1 = 440;
    % Min_CH_CH1 = round(855.1-21.45*mod(calib_CH-1,32))-12;%��1�s�[�N�p�`�����l�����؂���ŏ��l(1~1024)
    Min_L_CH1 = 370;%��1�s�[�N�p�g�����؂���ŏ��l(1~1024)
    Min_L_CH2 = 520;%��2�s�[�N�p�g�����؂���ŏ��l(1~1024)
    Min_L_CH3 = 330;%��2�s�[�N�p�g�����؂���ŏ��l(1~1024)
    Max_CH_CH1 = Min_CH_CH1+80;%��1�s�[�N�p�`�����l�����؂���ő�l(1~1024)
    Min_CH_CH2 = Min_CH_CH1;%��2�s�[�N�p�`�����l�����؂���ŏ��l(1~1024)
    Min_CH_CH3 = Min_CH_CH1;%��2�s�[�N�p�`�����l�����؂���ŏ��l(1~1024)

    Max_CH_CH2 = Max_CH_CH1;%��2�s�[�N�p�`�����l�����؂���ő�l(1~1024)
    Max_CH_CH3 = Max_CH_CH1;%��3�s�[�N�p�`�����l�����؂���ő�l(1~1024)
    Max_L_CH1 = Min_L_CH1+100;%��1�s�[�N�p�g�����؂���ő�l(1~1024)
    Max_L_CH2 = Min_L_CH2+100;%��2�s�[�N�p�g�����؂���ő�l(1~1024)
    Max_L_CH3 = Min_L_CH3+50;%��3�s�[�N�p�g�����؂���ő�l(1~1024)

    Min_L_L1 = Min_L_CH1;%��1�s�[�N�p�g�����؂���ŏ��l(1~1024)
    Max_L_L1 = Max_L_CH1;%��1�s�[�N�p�g�����؂���ő�l(1~1024)
    Min_L_L2 = Min_L_CH2;%��2�s�[�N�p�g�����؂���ŏ��l(1~1024)
    Max_L_L2 = Max_L_CH2;%��2�s�[�N�p�g�����؂���ő�l(1~1024)
    Min_L_L3 = Min_L_CH3;%��3�s�[�N�p�g�����؂���ŏ��l(1~1024)
    Max_L_L3 = Max_L_CH3;%��3�s�[�N�p�g�����؂���ő�l(1~1024)

    %-----�Z���t�@�C���ǂݍ���----
    cal_data = importdata(cal_filename);
    cal_data = cal_data(:,2:1025);%�f�[�^1��ڂ͒ʂ��ԍ��Ȃ̂Ő؂�̂Ă�

    %----�z���`----
    ax_pixel = transpose(linspace(1,1024,1024));%1~1024�̐�����
    cal_result = zeros(3,4);%�Z������

    %-----ICCD���f�[�^���v���b�g-----
    if plot_cont
        figure
        contour(ax_pixel,ax_pixel,cal_data)
        xlabel('Pixel number in CH direction')
        ylabel('Pixel number in Lambda direction')
        hold off
    end

    figure('Position',[300 50 1200 800],'visible','on')
    %-------��1�s�[�N�pCH�ʒu����------
    if cal_CH1
        spectrum_CH1 = ...
            transpose(sum(cal_data(Min_L_CH1:Max_L_CH1,Min_CH_CH1:Max_CH_CH1),1));
        % Offset_CH1 = mean(spectrum_CH1(1:3,1));
        Offset_CH1 = ...
            mean(cal_data([Min_L_CH1-5:Min_L_CH1, Max_L_CH1:Max_L_CH1+5],...
            [Min_CH_CH1-5:Min_CH_CH1, Max_CH_CH1:Max_CH_CH1+5]),'all') * (Max_L_CH1 - Min_L_CH1 + 1);
        Y_CH1 = spectrum_CH1 - Offset_CH1;
        f_CH1 = fit(ax_pixel(Min_CH_CH1:Max_CH_CH1),Y_CH1,'gauss1');
        % figure
        subplot(2,3,1)
        plot(f_CH1,ax_pixel(Min_CH_CH1:Max_CH_CH1),Y_CH1)
        title(['Detecting CH position' newline 'for the 1st peak'])
        legend('off')
        xlabel('Pixel number in CH direction')
        ylabel('Intensity [cnt]')
        Coef_CH1 = coeffvalues(f_CH1);
        cal_result(1,1) = sort(Coef_CH1(1,2),'descend');%�c�ʒu��傫�����Ŏ擾
    end

    %-------��2�s�[�N�pCH�ʒu����------
    if cal_CH2
        spectrum_CH2 = ...
            transpose(sum(cal_data(Min_L_CH2:Max_L_CH2,Min_CH_CH2:Max_CH_CH2),1));
        % Offset_CH2 = mean(spectrum_CH2(1:5,1));
        Offset_CH2 = ...
            mean(cal_data([Min_L_CH2-5:Min_L_CH2, Max_L_CH2:Max_L_CH2+5],...
            [Min_CH_CH2-5:Min_CH_CH2, Max_CH_CH2:Max_CH_CH2+5]),'all') * (Max_L_CH2 - Min_L_CH2 + 1);
        Y_CH2 = spectrum_CH2 - Offset_CH2;
        % Y_CH2 = movmean(Y_CH2,15);
        f_CH2 = fit(ax_pixel(Min_CH_CH2:Max_CH_CH2),Y_CH2,'gauss1');
        % figure
        subplot(2,3,2)
        plot(f_CH2,ax_pixel(Min_CH_CH2:Max_CH_CH2),Y_CH2)
        title(['Detecting CH position' newline 'for the 2nd peak'])
        legend('off')
        xlabel('Pixel number in CH direction')
        ylabel('Intensity [cnt]')
        Coef_CH2 = coeffvalues(f_CH2);
        cal_result(2,1) = sort(Coef_CH2(1,2),'descend');%�c�ʒu��傫�����Ŏ擾
    end

    %-------��3�s�[�N�pCH�ʒu����------
    if cal_CH3
        spectrum_CH3 = ...
            transpose(sum(cal_data(Min_L_CH3:Max_L_CH3,Min_CH_CH3:Max_CH_CH3),1));
        % Offset_CH3 = mean(spectrum_CH3(1:5,1));
        Offset_CH3 = ...
            mean(cal_data([Min_L_CH3-5:Min_L_CH3, Max_L_CH3:Max_L_CH3+5],...
            [Min_CH_CH3-5:Min_CH_CH3, Max_CH_CH3:Max_CH_CH3+5]),'all') * (Max_L_CH3 - Min_L_CH3 + 1);
        Y_CH3 = spectrum_CH3 - Offset_CH3;
        % Y_CH3 = movmean(Y_CH3,15);
        f_CH3 = fit(ax_pixel(Min_CH_CH3:Max_CH_CH3),Y_CH3,'gauss1');
        % figure
        subplot(2,3,3)
        plot(f_CH3,ax_pixel(Min_CH_CH3:Max_CH_CH3),Y_CH3)
        title(['Detecting CH position' newline 'for the 3rd peak'])
        legend('off')
        xlabel('Pixel number in CH direction')
        ylabel('Intensity [cnt]')
        Coef_CH3 = coeffvalues(f_CH3);
        cal_result(3,1) = sort(Coef_CH3(1,2),'descend');%�c�ʒu��傫�����Ŏ擾
    end

    %-------�g���ʒu����--------
    %��1Ne�s�[�N���o
    if cal_1st
        % figure
        subplot(2,3,4)
        Min_CH_L1 = round(cal_result(1,1)-width);%�`�����l�����؂���ŏ��l
        Max_CH_L1 = round(cal_result(1,1)+width);%�`�����l�����؂���ő�l
        spectrum_L1 = ...
            sum(cal_data(Min_L_L1:Max_L_L1,Min_CH_L1:Max_CH_L1),2);
        % Offset_L1 = mean(spectrum_L1(1:10,1));
        Offset_L1 = ...
            mean(cal_data([Min_L_L1-30:Min_L_L1, Max_L_L1:Max_L_L1+30],...
            [Min_CH_L1-30:Min_CH_L1, Max_CH_L1:Max_CH_L1+30]),'all') * (Max_CH_L1 - Min_CH_L1 + 1);
        Y_L1 = spectrum_L1 - Offset_L1;
        Y_L1 = movmean(Y_L1,l_mov);
        MAX1 = max(Y_L1);
        S1 = [ax_pixel(Min_L_L1:Max_L_L1) Y_L1]; %[�g��,���x]
        s1 = size(S1);
        deleted_S1 = zeros(1,2);
        ori_S1 = S1;
        j = 1;
        while j < s1(1)+1 %SN�̈����f�[�^������
            if S1(j,2) < MAX1*th_ratio
                deleted_S1 = cat(1,deleted_S1,S1(j,:));
                S1(j,:) = [];
            else
                j = j+1;
            end
            s1 = size(S1);
        end
        % f_L1 = fit(ax_pixel(Min_L_L1:Max_L_L1),Y_L1,'gauss1');
        % plot(f_L1,ax_pixel(Min_L_L1:Max_L_L1),Y_L1)
        f_L1 = fit(S1(:,1),S1(:,2),'gauss1');
        plot(f_L1,'r-',ori_S1(:,1),ori_S1(:,2),'.w');
        hold on
        plot(S1(:,1),S1(:,2),'bo');
        hold on
        plot(deleted_S1(:,1),deleted_S1(:,2),'kx');
        hold on
        yline(MAX1*th_ratio,'g','LineWidth',3);
        title(['Detecting Lambda position' newline 'for the 1st peak'])
        legend('off')
        xlabel('Pixel number in Lambda direction')
        ylabel('Intensity [cnt]')
        xlim([Min_L_L1 Max_L_L1])
        Coef_L1 = coeffvalues(f_L1);
        cal_result(1,2) = Coef_L1(1,2);%��1�s�[�N���擾
        cal_result(1,3) = Coef_L1(1,3);%��1�s�[�N��(���u�֐�)���擾
        cal_result(1,4) = Coef_L1(1,1)/1e5;%���Ί��x���擾
    end

    %��2Ne�s�[�N���o
    if cal_2nd
        % figure
        subplot(2,3,5)
        Min_CH_L2 = round(cal_result(2,1)-width);%�`�����l�����؂���ŏ��l
        Max_CH_L2 = round(cal_result(2,1)+width);%�`�����l�����؂���ő�l
        spectrum_L2 = ...
            sum(cal_data(Min_L_L2:Max_L_L2,Min_CH_L2:Max_CH_L2),2);
        % Offset_L2 = mean(spectrum_L2(1:10,1));
        Offset_L2 = ...
            mean(cal_data([Min_L_L2-30:Min_L_L2, Max_L_L2:Max_L_L2+30],...
            [Min_CH_L2-30:Min_CH_L2, Max_CH_L2:Max_CH_L2+30]),'all') * (Max_CH_L2 - Min_CH_L2 + 1);
        Y_L2 = spectrum_L2 - Offset_L2;
        Y_L2 = movmean(Y_L2,l_mov);
        MAX2 = max(Y_L2);
        S2 = [ax_pixel(Min_L_L2:Max_L_L2) Y_L2]; %[�g��,���x]
        s1 = size(S2);
        deleted_S2 = zeros(1,2);
        ori_S2 = S2;
        j = 1;
        while j < s1(1)+1 %SN�̈����f�[�^������
            if S2(j,2) < MAX2*th_ratio
                deleted_S2 = cat(1,deleted_S2,S2(j,:));
                S2(j,:) = [];
            else
                j = j+1;
            end
            s1 = size(S2);
        end
        % f_L2 = fit(ax_pixel(Min_L_L2:Max_L_L2),Y_L2,'gauss1');
        % plot(f_L2,ax_pixel(Min_L_L2:Max_L_L2),Y_L2)
        f_L2 = fit(S2(:,1),S2(:,2),'gauss1');
        plot(f_L2,'r-',ori_S2(:,1),ori_S2(:,2),'.w');
        hold on
        plot(S2(:,1),S2(:,2),'bo');
        hold on
        plot(deleted_S2(:,1),deleted_S2(:,2),'kx');
        hold on
        yline(MAX2*th_ratio,'g','LineWidth',3);
        title(['Detecting Lambda position' newline 'for the 2nd peak'])
        legend('off')
        xlabel('Pixel number in Lambda direction')
        ylabel('Intensity [cnt]')
        xlim([Min_L_L2 Max_L_L2])
        Coef_L2 = coeffvalues(f_L2);
        cal_result(2,2) = Coef_L2(1,2);%��2�s�[�N���擾
        cal_result(2,3) = Coef_L2(1,3);%��2�s�[�N��(���u�֐�)���擾
        cal_result(2,4) = Coef_L2(1,1)/1e5;%���Ί��x���擾
        px2nm = -(lambda_1st - lambda_2nd)/(cal_result(1,2) - cal_result(2,2));%px2nm���v�Z
    end

    %��3Ne�s�[�N���o
    if cal_3rd
        % figure
        subplot(2,3,6)
        Min_CH_L3 = round(cal_result(3,1)-width);%�`�����l�����؂���ŏ��l
        Max_CH_L3 = round(cal_result(3,1)+width);%�`�����l�����؂���ő�l
        spectrum_L3 = ...
            sum(cal_data(Min_L_L3:Max_L_L3,Min_CH_L3:Max_CH_L3),2);
        % Offset_L3 = mean(spectrum_L3(1:10,1));
        Offset_L3 = ...
            mean(cal_data([Min_L_L3-30:Min_L_L3, Max_L_L3:Max_L_L3+30],...
            [Min_CH_L3-30:Min_CH_L3, Max_CH_L3:Max_CH_L3+30]),'all') * (Max_CH_L3 - Min_CH_L3 + 1);
        Y_L3 = spectrum_L3 - Offset_L3;
        Y_L3 = movmean(Y_L3,l_mov);
        MAX3 = max(Y_L3);
        S3 = [ax_pixel(Min_L_L3:Max_L_L3) Y_L3]; %[�g��,���x]
        s1 = size(S3);
        deleted_S3 = zeros(1,2);
        ori_S3 = S3;
        j = 1;
        while j < s1(1)+1 %SN�̈����f�[�^������
            if S3(j,2) < MAX3*th_ratio
                deleted_S3 = cat(1,deleted_S3,S3(j,:));
                S3(j,:) = [];
            else
                j = j+1;
            end
            s1 = size(S3);
        end
        % f_L3 = fit(ax_pixel(Min_L_L3:Max_L_L3),Y_L3,'gauss1');
        % plot(f_L3,ax_pixel(Min_L_L3:Max_L_L3),Y_L3)
        f_L3 = fit(S3(:,1),S3(:,2),'gauss1');
        plot(f_L3,'r-',ori_S3(:,1),ori_S3(:,2),'.w');
        hold on
        plot(S3(:,1),S3(:,2),'bo');
        hold on
        plot(deleted_S3(:,1),deleted_S3(:,2),'kx');
        hold on
        yline(MAX3*th_ratio,'g','LineWidth',3);
        title(['Detecting Lambda position' newline 'for the 2nd peak'])
        legend('off')
        xlabel('Pixel number in Lambda direction')
        ylabel('Intensity [cnt]')
        xlim([Min_L_L3 Max_L_L3])
        Coef_L3 = coeffvalues(f_L3);
        cal_result(3,2) = Coef_L3(1,2);%��2�s�[�N���擾
        cal_result(3,3) = Coef_L3(1,3);%��2�s�[�N��(���u�֐�)���擾
        cal_result(3,4) = Coef_L3(1,1)/1e5;%���Ί��x���擾
    end

    sgtitle(['CH',num2str(calib_CH)])

    if save_fit
        if not(exist(num2str(date),'dir'))
            mkdir(num2str(date));
        end
        saveas(gcf,[num2str(date),'/fit_CH',num2str(calib_CH),'.png'])
        hold off
        close
    end
    % cal_result
    Lambda = [lambda_1st lambda_2nd lambda_3rd];
    figure
    x_poly = cal_result(1:3,2);
    y_poly = Lambda(1:3);
    p = polyfit(x_poly,y_poly,2);
    x2 = linspace(1,1024,1024);
    y2 = polyval(p,x2);
    plot(x_poly,y_poly,'o',x2,y2)
    % calib_output = [calib_CH cal_result(1,1:5)]
    if cal_1st && save_cal
        if not(exist(num2str(date),'dir'))
            mkdir(num2str(date));
        end
        save([num2str(date),'/calibation',num2str(calib_CH),'.mat'],'calib_output')
    end
end