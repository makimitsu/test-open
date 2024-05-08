function [IDSPdata] = cal_ionflow(IDSP,pathname,show_offset,plot_fit,save_fit)
%IDSP�ϐ�/pathname/offset��\��/�K�E�X�t�B�b�e�B���O���v���b�g/������ۑ�
%------�t�B�b�e�B���O�����p�yinput�z-------
w_CH = 6;%�yinput�z�`�����l������(Y����)�������킹��[px]6
l_mov_flow = 15;%�yinput�z�g�������ړ����ϒ���[px]15
th_ratio_flow = 0.8;%�yinput�z�t�B�b�e�B���O��臒l[px]0.8
l_mov_temp = l_mov_flow ;%�yinput�z�g�������ړ����ϒ���[px]15
th_ratio_temp = th_ratio_flow;%�yinput�z�t�B�b�e�B���O��臒l[px]0.8
l_L1 = 61;%�yinput�z�g�����̐؂��蒷��[px]
d_L1 = 38;%�yinput�z�g�������ʒu���꒲��[px]
d_CH = -7;%�yinput�zCH�����ʒu���꒲��[px]
pre_offset = 0.086;%�yinput�z������g�������ʒu���꒲��[nm]

savename = [pathname.mat,'/ionflow/',num2str(IDSP.date),'_shot',num2str(IDSP.shot),'_',num2str(IDSP.delay),'us_w=',num2str(IDSP.width),'_gain=',num2str(IDSP.gain),'.mat'];
if exist(savename,"file")
    load(savename,'IDSPdata')
else
    %�����萔
    Vc = 299792.458;%����(km/s)
    Angle = 30;%���ˊp(�x)
    switch char(IDSP.line)
        case 'Ar'%�A���S���̎�
            A = 40;%���q��
            lambda0 = 480.602;%�g�p�X�y�N�g��(nm)
            lambda1 = 480.7019;%�Z�������v�X�y�N�g��(nm)
            lambda2 = 479.2619;%�Z�������v�X�y�N�g��(nm)
        case 'H'%���f�̎�
            A = 1;%���q��
            lambda0 = 486.135;%�g�p�X�y�N�g��(nm)
            warning('Sorry, not ready for H experiment.')%ICCD.line�̓��̓G���[
            return;
        otherwise
            warning('Input error in ICCD.line.')%ICCD.line�̓��̓G���[
            return;
    end
    center = load_IDSPcalib(IDSP);

    dir1 = [pathname.IDSP '/' num2str(IDSP.date)];%�f�B���N�g��1
    if IDSP.n_z == 1
        filename1 = [dir1 '/shot' num2str(IDSP.shot) '_' num2str(IDSP.delay) 'us_w=' num2str(IDSP.width) '_gain=' num2str(IDSP.gain) '.asc'];%ICCD�t�@�C����
        if not(exist(filename1,"file"))
            warning([filename1,' does not exist.']);
            return
        end
        %data1����X�y�N�g���Q���擾
        data1 = importdata(filename1);
    end

    %----------------------------------------------------------------------------
    % centerX = repmat(center(:,3),1,IDSP.n_z);%�`�����l���Ή����S����X���W
    switch IDSP.n_z
        case 1
            centerY = center(:,2);%�`�����l���Ή����SY���W
        case 2
            warning('Sorry, not ready for n_z = 2.')%ICCD.line�̓��̓G���[
            return
            centerY = repmat(center(:,2),1,IDSP.n_z);%�`�����l���Ή����SY���W
    end

    X1 = data1(:,1);%X1(�s�N�Z��)�����`
    [l_X1,~]=size(data1);%X1���̒������擾
    L1 = zeros(l_X1,IDSP.n_CH);%L1(�g��)�����`
    L1_shaped = zeros(l_L1,IDSP.n_CH);%�؂������g����
    px2nm = zeros(IDSP.n_CH,1);%nm/pixel
    switch char(IDSP.line)
        case 'Ar'%�A���S���̎�
            %     lambda = [lambda1 lambda2];
            for i = 1:IDSP.n_CH
                %         pixel = [center(i,3) center(i,4)];
                %         p = polyfit(pixel,lambda,1);
                %         px2nm(i,1) = p(1);
                %         L1(:,i) = polyval(p,X1);
                px2nm(i,1) = center(i,4);
                L1(:,i) = lambda1 - px2nm(i,1)*(X1(:,1)-center(i,3));
                L1_shaped(:,i) = L1(round(center(i,3))-(l_L1-1)/2+d_L1:round(center(i,3))+(l_L1-1)/2+d_L1,i)+pre_offset;
            end
        case 'H'%���f�̎�
            for i = 1:IDSP.n_CH
                px2nm(i,1) = 0.00536;
                L1(:,i) = px2nm(i,1)*(X1-center(i,3))+lambda0;
            end
    end
    spectrum1 = zeros(l_L1,IDSP.n_CH);%data1�̕������ʂ�����
    for i = 1:IDSP.n_CH
        spectrum1(:,i) = ...
            sum(data1(round(center(i,3))-(l_L1-1)/2+d_L1:round(center(i,3))+(l_L1-1)/2+d_L1,round(centerY(i,1)+d_CH-w_CH):round(centerY(i,1)+d_CH+w_CH)),2);
        spectrum1(:,i) = spectrum1(:,i)./center(i,6);
    end

    % %data2����X�y�N�g���Q���擾
    % if IDSP.n_z > 1
    %     data2 = importdata(filename2);
    %     X2 = data2(:,1);%X2(�s�N�Z��)�����`
    %     [LofX2,~]=size(data2);%X2���̒������擾
    %     spectrum2=zeros(LofX2,IDSP.n_CH);%data2�̕������ʂ�����
    %     for i = 1:IDSP.n_CH
    %         spectrum2(:,i) = ...
    %             sum(data2(:,round(centerY(i,2)-width):round(centerY(i,2)+width)),2);
    %     end
    % end

    %�v�Z���ʎ擾�p�̔z���p��
    amp = zeros(IDSP.n_CH,IDSP.n_z);%�U��1���data1
    shift = zeros(IDSP.n_CH,IDSP.n_z);%���S1���data1
    err_shift = zeros(IDSP.n_CH,IDSP.n_z);%���S1���data1
    sigma = zeros(IDSP.n_CH,IDSP.n_z);%�V�O�}1���data1
    err_sigma = zeros(IDSP.n_CH,IDSP.n_z);%�V�O�}1���data1
    IDSPdata.V_CH = zeros(IDSP.n_CH,IDSP.n_z);%1���data1�EV(km/s)
    IDSPdata.err_V_CH = zeros(IDSP.n_CH,IDSP.n_z);%1���data1�EV(km/s)
    IDSPdata.V_i = zeros(IDSP.n_r,IDSP.n_z*2);%1���data1�EVz(km/s)�A2���data1�EVr(km/s)
    IDSPdata.err_V_i = zeros(IDSP.n_r,IDSP.n_z*2);%1���data1�EVz(km/s)�A2���data1�EVr(km/s)
    IDSPdata.absV = zeros(IDSP.n_r,IDSP.n_z);%1���data1�EV(km/s)
    IDSPdata.T_CH = zeros(IDSP.n_CH,IDSP.n_z);%1���data1�E���x(eV)
    IDSPdata.err_T_CH = zeros(IDSP.n_CH,IDSP.n_z);%1���data1�E���x(eV)
    IDSPdata.T_i = zeros(IDSP.n_r,IDSP.n_z);%1���data1�E���x(eV)
    IDSPdata.err_T_i = zeros(IDSP.n_r,IDSP.n_z);%1���data1�E���x(eV)
    IDSPdata.trim_T_i = zeros(IDSP.n_r,IDSP.n_z);%1���data1�E���x(eV)
    IDSPdata.err_trim_T_i = zeros(IDSP.n_r,IDSP.n_z);%1���data1�E���x(eV)
    IDSPdata.offset = zeros(IDSP.n_r,IDSP.n_z);
    error_CH = zeros(IDSP.n_CH,1);

    %�g�������̈ړ����ς��痣�ꂽ�O��l���ړ����ςɕύX(�m�C�Y����)
    M1_flow = movmean(spectrum1,l_mov_flow);
    spectrum1_flow = spectrum1;
    for i = 1:IDSP.n_CH
        for j = 1:l_L1
            if spectrum1_flow(j,i) > M1_flow(j,i)*1.05
                spectrum1_flow(j,i) = M1_flow(j,i);
            elseif spectrum1_flow(j,i) < M1_flow(j,i)*0.95
                spectrum1_flow(j,i) = M1_flow(j,i);
            end
        end
    end

    for k = 1:IDSP.n_r
        %0�x�y�A�X�y�N�g������I�t�Z�b�g�����o
        for i = 1:2
            i_CH = (k-1)*4+i;%CH�ԍ�
            S1 = [L1_shaped(:,i_CH) spectrum1_flow(:,i_CH)]; %[�g��,���x]
            s1 = size(S1);%S1��[�s��,��]
            MAX1 = max(spectrum1_flow(:,i_CH)); %�X�y�N�g���̍ő�l
            j = 1;
            while j < s1(1)+1 %SN�̈����f�[�^������
                if S1(j,2) < MAX1*th_ratio_flow
                    S1(j,:) = [];
                else
                    j = j+1;
                end
                s1 = size(S1);
            end
            try
                f = fit(S1(:,1),S1(:,2),'gauss1');
            catch ME
                warning(['Fitting failed in CH ',num2str(i_CH),'.']);
                error_CH(i_CH,1) = 1;
                break
            end
            coef=coeffvalues(f);
            shift(i_CH,1) = coef(2)-lambda0;
        end
        IDSPdata.offset(k,1) = (shift((k-1)*4+1,1) + shift((k-1)*4+2,1))/2;%�Ό��������瓾��ꂽ�I�t�Z�b�g[nm]
    end

    % %�ُ�l������offset����}���Ēu������
    % ori_offset = IDSPdata.offset;%�ύX�O��ۊ�
    % buf = IDSPdata.offset(:,1);%offset���R�s�[
    % idx_0 = find(buf == 0);%0�����o
    % idx_non0 = find(buf ~= 0);%non0�����o
    % buf(idx_0) = [];%0���폜
    % buf = filloutliers(buf,"linear");%�O��l����^���}�Œu������
    % IDSPdata.offset(idx_non0,1) = buf;%offset���X�V
    % IDSPdata.offset(idx_0,1) = mean(IDSPdata.offset(idx_non0,1));%offset���X�V
    % for k = 1:IDSP.n_r
    %     if IDSPdata.offset(k,1) ~= ori_offset(k,1)
    %         disp(['Offset in Row ',num2str(k),' is replaced from ',num2str(ori_offset(k,1)),' to ',num2str(IDSPdata.offset(k,1))'.'])
    %     end
    % end


    % offset_mean = (IDSPdata.offset(2,1) + IDSPdata.offset(3,1))/2;
    for k = [1 5:7]
        IDSPdata.offset(k,1) = IDSPdata.offset(3,1);
    end
    % IDSPdata.offset(1,1) = IDSPdata.offset(2,1);%1�Ԃ̎���2���g�`����offset���M�p�ł��Ȃ����ߓ���(20230807�Z���̂Ƃ�)
    % IDSPdata.offset(6,1) = IDSPdata.offset(5,1);
    % IDSPdata.offset(5,1) = IDSPdata.offset(6,1);%5�Ԃ̎���1�Ǝ���2���A�ԂłȂ�offset���M�p�ł��Ȃ����ߓ���(20230807�Z���̂Ƃ�)

    %data1�̗����v�Z�p�K�E�X�t�B�b�e�B���O
    if plot_fit
        if save_fit
            figure('Position',[300 50 1000 1000],'visible','off')
        else
            figure('Position',[300 50 1000 1000],'visible','on')
        end
        sgtitle(['Fitting data for flow (Horizontal�FView Line, Vertical�FMeasured Position)',newline, ...
            'shot',num2str(IDSP.shot),'-',num2str(IDSP.delay),'us-w=',num2str(IDSP.width),'-gain=',num2str(IDSP.gain),'.asc'])
    end
    for k = 1:IDSP.n_r
        % if (error_CH((k-1)*4+1,1) == 0) && (error_CH((k-1)*4+2,1) == 0)%�΍R�����̃t�B�b�e�B���O�ɐ���
        if show_offset
            disp(['Offset = ',num2str(IDSPdata.offset(k,1)),' in Row ',num2str(k),'.'])
        end
        %�I�t�Z�b�g�������ăK�E�X�t�B�b�e�B���O
        for i = 1:4
            i_CH = (k-1)*4+i;%CH�ԍ�
            S1 = [L1_shaped(:,i_CH)-IDSPdata.offset(k,1) spectrum1_flow(:,i_CH)]; %[�g��,���x]
            s1 = size(S1);%S1��[�s��,��]
            ori_S1 = S1;
            deleted_S1 = zeros(1,2);
            MAX1 = max(spectrum1_flow(:,i_CH)); %�X�y�N�g���̍ő�l
            j = 1;
            while j < s1(1)+1 %SN�̈����f�[�^������
                if S1(j,2) < MAX1*th_ratio_flow
                    deleted_S1 = cat(1,deleted_S1,S1(j,:));
                    S1(j,:) = [];
                else
                    j = j+1;
                end
                s1 = size(S1);
            end
            try
                f = fit(S1(:,1),S1(:,2),'gauss1');
                coef=coeffvalues(f);
                confi = confint(f);%�t�B�b�e�B���O�W����95%�M�����(1�s�ډ���, 2�s�ڏ��)
                amp(i_CH,1) = coef(1);
                shift(i_CH,1) = coef(2)-lambda0;
                err_shift(i_CH,1) = confi(2,2)-coef(2);
                IDSPdata.V_CH(i_CH,1) = shift(i_CH,1)/lambda0*Vc;
                IDSPdata.err_V_CH(i_CH,1) = err_shift(i_CH,1)/lambda0*Vc;
                if (l_mov_flow == l_mov_temp) && (th_ratio_flow == th_ratio_temp)
                    sigma(i_CH,1) = sqrt(coef(3)^2-(center(i_CH,5)*px2nm(i_CH,1))^2);
                    err_sigma(i_CH,1) = sqrt(confi(2,3)^2-(center(i_CH,5)*px2nm(i_CH,1))^2);
                    IDSPdata.T_CH(i_CH,1) = 1.69e8*A*(2*sigma(i_CH,1)*sqrt(log(2))/lambda0)^2;
                    IDSPdata.err_T_CH(i_CH,1) = 1.69e8*A*(2*err_sigma(i_CH,1)*sqrt(log(2))/lambda0)^2 - IDSPdata.T_CH(i_CH,1);
                end
                if plot_fit
                    subplot(IDSP.n_r,4,i_CH);
                    plot(f,'r-',ori_S1(:,1),ori_S1(:,2),'.w');
                    hold on
                    plot(S1(:,1),S1(:,2),'bo');
                    hold on
                    plot(deleted_S1(:,1),deleted_S1(:,2),'kx');
                    hold on
                    xline(lambda0);
                    yline(MAX1*th_ratio_flow,'g','LineWidth',3);
                    title(sprintf('%2.1f �} %2.1f km/s',IDSPdata.V_CH(i_CH,1),IDSPdata.err_V_CH(i_CH,1)))
                    legend('off')
                    xlabel('Wavelength [nm]')
                    ylabel('Intensity [cnt]')
                    xlim([lambda0-0.1,lambda0+0.1])
                    ylim([0,MAX1*1.1])
                    legend('off')
                end
            catch ME
                warning(['Fitting failed in CH ',num2str(i_CH),'.']);
                break
            end
        end
        % else
        %     warning(['Offset failed in Row ',num2str(k),'.']);
        % end
    end
    if plot_fit
        if save_fit
            time = round(IDSP.delay+IDSP.width/2);%�v������
            if (l_mov_flow == l_mov_temp) && (th_ratio_flow == th_ratio_temp)
                saveas(gcf,[pathname.fig,'/fit/',num2str(IDSP.date),'_shot', num2str(IDSP.shot),'_',num2str(time),'us_fit_for_both.png'])
            else
                saveas(gcf,[pathname.fig,'/fit/',num2str(IDSP.date),'_shot', num2str(IDSP.shot),'_',num2str(time),'us_fit_for_flow.png'])
            end
            hold off
            close
        else
            hold off
        end
    end

    if not((l_mov_flow == l_mov_temp) && (th_ratio_flow == th_ratio_temp))
        %�g�������̈ړ����ς��痣�ꂽ�O��l���ړ����ςɕύX(�m�C�Y����)
        M1_temp = movmean(spectrum1,l_mov_temp);
        spectrum1_temp = spectrum1;
        for i = 1:IDSP.n_CH
            for j = 1:l_L1
                if spectrum1_temp(j,i) > M1_temp(j,i)*1.05
                    spectrum1_temp(j,i) = M1_temp(j,i);
                elseif spectrum1_temp(j,i) < M1_temp(j,i)*0.95
                    spectrum1_temp(j,i) = M1_temp(j,i);
                end
            end
        end

        %data1�̉��x�v�Z�p�K�E�X�t�B�b�e�B���O
        if plot_fit
            if save_fit
                figure('Position',[300 50 1000 1000],'visible','off')
            else
                figure('Position',[300 50 1000 1000],'visible','on')
            end
            sgtitle(['Fitting data for temp. (Horizontal�FView Line, Vertical�FMeasured Position)',newline, ...
                'shot',num2str(IDSP.shot),'-',num2str(IDSP.delay),'us-w=',num2str(IDSP.width),'-gain=',num2str(IDSP.gain),'.asc'])
        end
        for k = 1:IDSP.n_r
            for i = 1:4
                i_CH = (k-1)*4+i;%CH�ԍ�
                S1 = [L1_shaped(:,i_CH)-IDSPdata.offset(k,1) spectrum1_temp(:,i_CH)]; %[�g��,���x]
                s1 = size(S1);%S1��[�s��,��]
                ori_S1 = S1;
                deleted_S1 = zeros(1,2);
                MAX1 = max(spectrum1_temp(:,i_CH)); %�X�y�N�g���̍ő�l
                j = 1;
                while j < s1(1)+1 %SN�̈����f�[�^������
                    if S1(j,2) < MAX1*th_ratio_temp
                        deleted_S1 = cat(1,deleted_S1,S1(j,:));
                        S1(j,:) = [];
                    else
                        j = j+1;
                    end
                    s1 = size(S1);
                end
                try
                    f = fit(S1(:,1),S1(:,2),'gauss1');
                    coef=coeffvalues(f);
                    confi = confint(f);%�t�B�b�e�B���O�W����95%�M�����(1�s�ډ���, 2�s�ڏ��)
                    sigma(i_CH,1) = sqrt(coef(3)^2-(center(i_CH,5)*px2nm(i_CH,1))^2);
                    err_sigma(i_CH,1) = sqrt(confi(2,3)^2-(center(i_CH,5)*px2nm(i_CH,1))^2);
                    IDSPdata.T_CH(i_CH,1) = 1.69e8*A*(2*sigma(i_CH,1)*sqrt(log(2))/lambda0)^2;
                    IDSPdata.err_T_CH(i_CH,1) = 1.69e8*A*(2*err_sigma(i_CH,1)*sqrt(log(2))/lambda0)^2 - IDSPdata.T_CH(i_CH,1);
                    if plot_fit
                        subplot(IDSP.n_r,4,i_CH);
                        plot(f,'r-',ori_S1(:,1),ori_S1(:,2),'.w');
                        hold on
                        plot(S1(:,1),S1(:,2),'bo');
                        hold on
                        plot(deleted_S1(:,1),deleted_S1(:,2),'kx');
                        hold on
                        xline(lambda0);
                        yline(MAX1*th_ratio_temp,'g','LineWidth',3);
                        title(sprintf('%3.0f �} %3.0f eV',IDSPdata.T_CH(i_CH,1),IDSPdata.err_T_CH(i_CH,1)))
                        legend('off')
                        xlabel('Wavelength [nm]')
                        ylabel('Intensity [cnt]')
                        xlim([lambda0-0.1,lambda0+0.1])
                        ylim([0,MAX1*1.1])
                        legend('off')
                    end
                catch ME
                    warning(['Fitting failed in CH ',num2str(i_CH),'.']);
                    break
                end
            end
            % else
            %     warning(['Offset failed in Row ',num2str(k),'.']);
            % end
        end
        if plot_fit
            if save_fit
                % time = round(IDSP.delay+IDSP.width/2);%�v������
                % saveas(gcf,[pathname.fig,'/fit/',num2str(IDSP.date),'_shot', num2str(IDSP.shot),'_',num2str(time),'us_fit_for_temp.png'])
                hold off
                close
            else
                hold off
            end
        end
    end
    %data1�̗����A���x���v�Z
    for i = 1:IDSP.n_r
        set = (i-1)*4;
        Va = -shift(set+4,1)/lambda0*Vc;
        Vb = -shift(set+3,1)/lambda0*Vc;
        Vz2 = -shift(set+2,1)/lambda0*Vc;
        err_Va = err_shift(set+4,1)/lambda0*Vc;
        err_Vb = err_shift(set+3,1)/lambda0*Vc;
        err_Vz2 = err_shift(set+2,1)/lambda0*Vc;
        %����3�Ǝ���4�Ōv�Z
        IDSPdata.V_i(i,1) = (Va-Vb)/(2*cos(Angle*pi/180));%Vz
        IDSPdata.V_i(i,2) = -(Va+Vb)/(2*sin(Angle*pi/180));%Vr
        IDSPdata.err_V_i(i,1) = (err_Va+err_Vb)/(2*cos(Angle*pi/180));%err_Vz
        IDSPdata.err_V_i(i,2) = (err_Va+err_Vb)/(2*sin(Angle*pi/180));%err_Vr
        IDSPdata.absV(i,1) = sqrt(IDSPdata.V_i(i,1)^2 + IDSPdata.V_i(i,2)^2);
        % if i == 5
        %     % ����2�Ǝ���4�Ōv�Z
        %     IDSPdata.V_i(i,1) = Vz2;%Vz
        %     IDSPdata.V_i(i,2) = Vz2/tan(Angle*pi/180)-Va/sin(Angle*pi/180);%Vr
        %     IDSPdata.err_V_i(i,1) = err_Vz2;%err_Vz
        %     IDSPdata.err_V_i(i,2) = err_Vz2/tan(Angle*pi/180)+err_Va/sin(Angle*pi/180);%err_Vr
        %     IDSPdata.absV(i,1) = sqrt(IDSPdata.V_i(i,1)^2 + IDSPdata.V_i(i,2)^2);
        % end
        % if i == 6
        %     %����2�Ǝ���3�Ōv�Z
        %     IDSPdata.V_i(i,1) = Vz2;%Vz
        %     IDSPdata.V_i(i,2) = -Vz2/tan(Angle*pi/180)-Vb/sin(Angle*pi/180);%Vr
        %     IDSPdata.err_V_i(i,1) = err_Vz2;%err_Vz
        %     IDSPdata.err_V_i(i,2) = err_Vz2/tan(Angle*pi/180)+err_Vb/sin(Angle*pi/180);%err_Vr
        %     IDSPdata.absV(i,1) = sqrt(IDSPdata.V_i(i,1)^2 + IDSPdata.V_i(i,2)^2);
        % end
        % IDSPdata.T_i(i,1) = mean(IDSPdata.T_CH(set+1:set+4,1));
        % IDSPdata.err_T_i(i,1) = std(IDSPdata.T_CH(set+1:set+4,1));
        IDSPdata.trim_T_i(i,1) = trimmean(IDSPdata.T_CH(set+1:set+4,1),50);
        sum_weight = 0;
        for j = 1:4
            weight = 1./IDSPdata.err_T_CH(set+j,1).^2;%�d��1/��^2
            IDSPdata.T_i(i,1) = IDSPdata.T_i(i,1) + IDSPdata.T_CH(set+j,1)*weight;
            sum_weight = sum_weight + weight;
        end
        IDSPdata.T_i(i,1) = IDSPdata.T_i(i,1)/sum_weight;
        IDSPdata.err_T_i(i,1) = 1./sqrt(sum_weight);
    end

    IDSPdata.z = IDSP.z;
    IDSPdata.r = IDSP.r;
    IDSPdata.delay = IDSP.delay;
    IDSPdata.width = IDSP.width;
    IDSPdata.gain = IDSP.gain;
    IDSPdata.time = IDSP.time;
    save(savename,'IDSPdata')
end

