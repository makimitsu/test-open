function [] = main(t,ndata,gas,nz,plot_fitting,cal_flow,save_flow,plot_flow,factor)
%����/�f�[�^�ԍ�/�K�X��(Ar:1,H:2)/z�����f�[�^��(���l)/�t�B�b�e�B���O��\��(TorF)/�������v�Z(TorF)/������ۑ�(TorF)/������\��(TorF)/���T�C�Y(���l)

%�p�����[�^���`
run define/parameter.m
r_measured = zeros(NofCH/4,nz);%�x�N�g���v���b�gr���W1���data1�A2���data2
z_measured = zeros(NofCH/4,nz);%�x�N�g���v���b�gz���W1���data1�A2���data2

%------------------------�t�@�C�����ƑΉ�����v���_�����--------------------------
run data_221220_TFscan.m
% z_measured(:,1) = [3 3 3 3 3 3 3 3 0 0 0 0 0 0 0 0 -3 -3 -3 -3 -3 -3 -3 -3];%data1��z���W
% r_measured(:,1) = [10 10 15 15 20 20 25 25 10 10 15 15 20 20 25 25 10 10 15 15 20 20 25 25];%data1��r���W
z_measured(:,1) = [3 3 3 3 3 3 3 3 0 0 0 0 0 0 0 0];%data1��z���W
r_measured(:,1) = [10 10 15 15 20 20 25 25 10 10 15 15 20 20 25 25];%data1��r���W
% z_measured(:,1) = 0;%data1��z���W
% r_measured(:,1) = [10 10 15 15 20 20 25 25];%data1��r���W
% r_measured(:,1) = [17.5 17.5 22.5 22.5 27.5 27.5 32.5 32.5];%data1��r���W
% filename1 = '/Users/itsuki/Downloads/shot22_470us_w=2_gain=3400.asc';%data1�t�@�C����
% filename1 = '/Users/itsuki/Documents/����/��������/Doppler_probe_data/221220/shot33_474us_w=2_gain=4095_30kV.asc';%data1�t�@�C����
filename1 = '/Volumes/experiment/results/Doppler/Andor/IDSP/230110/shot9-16.asc';%data1�t�@�C����
% filename1 = '/Volumes/experiment/results/Doppler/Andor/IDSP/230110/shot20_480us_w=3_gain=4095.asc';%data1�t�@�C����
% filename1 = '/Users/itsuki/Documents/����/��������/Doppler_probe_data/221221/sum-shot1-3.asc';%data1�t�@�C����
% z_measured(:,1) = -2.1;%data1��z���W
% r_measured(:,1) = [12.5, 15, 17.5, 20, 22.5, 25, 27.5];%data1��r���W
% filename2 = '/Users/Ryo/Doppler/211110/ascii/shot50_464us_w=2_gain=3800.asc';%data2�t�@�C����
% filename2 = '/Volumes/experiment/results/Doppler/Andor/IDSP/221101/shot31_474_w=2_3100.asc';%data2�t�@�C����
% z_measured(:,2) = 2.1;%data2��z���W
% r_measured(:,1) = [22.5, 25, 27.5, 30, 32.5, 35, 37.5];%data2��r���W
%----------------------------------------------------------------------------

%�r���f�[�^����v���[�u�e�`�����l���Ή��ӏ����擾
if gas == 1
    %�A���S���̎�
        center_file = 'calibation_96ch_230110_shot9.txt';%���S�f�[�^�t�@�C����%96ch�p
%     if NofCH == 32
%         center_file = 'calibation_32ch_221220.txt';%���S�f�[�^�t�@�C����
%     elseif NofCH == 96
%         center_file = 'calibation_96ch.txt';%���S�f�[�^�t�@�C����
%     end
elseif gas == 2
    %���f�̎�
    center_file = 'Hbeta_calibration.txt';%���S�f�[�^�t�@�C����
end

center = importdata(center_file);%���S���W���擾
% centerX = repmat(center(:,3),1,nz);%�`�����l���Ή����S����X���W
if nz == 1
    centerY = center(:,2);%�`�����l���Ή����SY���W
elseif nz == 2
    centerY = repmat(center(:,2),1,nz);%�`�����l���Ή����SY���W
end

%data1����X�y�N�g���Q���擾
data1 = importdata(filename1);
X1 = data1(:,1);%X1(�s�N�Z��)�����`
[LofX1,n1]=size(data1);%X1���̒������擾
haba = 30;
slide = 0;
SN = 0.75;
LofX1 = haba*2+1;%128ch�p
L1 = zeros(LofX1,NofCH);%L1(�g��)�����`
px2nm = zeros(NofCH,1);%nm/pixel
% load('magflow/mat/221209_Xe4834_calibation.mat','CH_number','CH_position','lambda1_position','nm/pixel','instrument')
if gas == 1
    %�A���S���̎�
%     lambda = [lambda1 lambda2];%96ch�p
    for i = 1:NofCH
%96ch������
%         pixel = [center(i,3) center(i,4)];
%         p = polyfit(pixel,lambda,1);
%         px2nm(i,1) = p(1);
%         L1(:,i) = polyval(p,X1);
%128ch������
        px2nm(i,1) = center(i,4);
        L1(:,i) = -px2nm(i,1)*(X1(round(center(i,3))-haba+slide:round(center(i,3))+haba+slide,1)-center(i,3))+lambda1;
        
    end
elseif gas == 2
    %���f�̎�
    for i = 1:NofCH
        px2nm(i,1) = 0.00536;
        L1(:,i) = px2nm(i,1)*(X1-center(i,3))+lambda0;
    end
end
spectrum1 = zeros(LofX1,NofCH);%data1�̕������ʂ�����
for i = 1:NofCH
    spectrum1(:,i) = ...
        sum(data1(round(center(i,3))-haba+slide:round(center(i,3))+haba+slide,round(centerY(i,1)-W):round(centerY(i,1)+W)),2);
end


%data2����X�y�N�g���Q���擾
if nz == 2
    data2 = importdata(filename2);
    X2 = data2(:,1);%X2(�s�N�Z��)�����`
    [LofX2,n2]=size(data2);%X2���̒������擾
    spectrum2=zeros(LofX2,NofCH);%data2�̕������ʂ�����
    for i = 1:NofCH
        spectrum2(:,i) = ...
            sum(data2(:,round(centerY(i,2)-W):round(centerY(i,2)+W)),2);
    end
end

%�g�������̈ړ����ς��痣�ꂽ�O��l���ړ����ςɕύX(�m�C�Y����)
M1 = movmean(spectrum1,l);
% M1 = (M1*l - spectrum1) / (l-1);
if nz == 2
    M2 = movmean(spectrum2,l);
    % M2 = (M2*l - spectrum2) / (l-1);
end
for i = 1:NofCH
    for j = 1:LofX1
        if spectrum1(j,i) > M1(j,i)*1.1
            spectrum1(j,i) = M1(j,i);
        elseif spectrum1(j,i) < M1(j,i)*0.9
            spectrum1(j,i) = M1(j,i);
        end
        if nz == 2
            if spectrum2(j,i) > M2(j,i)*1.1
                spectrum2(j,i) = M2(j,i);
            elseif spectrum2(j,i) < M2(j,i)*0.9
                spectrum2(j,i) = M2(j,i);
            end
        end
    end
end


% %�o�b�N�O���E���h�����ĂȂ��Ƃ��p�̃m�C�Y����
% for j = 1:NofCH
%     goukei = 0;
%     for i = 1:200
%         goukei = goukei + spectrum1(i,j);
%     end
%     heikin = goukei/200;
%     for i = 1:1024
%         spectrum1(i,j) = spectrum1(i,j) - heikin;
%     end
% end
% % spectrum1(1024,8) 


%�v�Z���ʎ擾�p�̔z���p��
amp = zeros(NofCH,nz);%�U��1���data1�A2���data2
shift = zeros(NofCH,nz);%���S1���data1�A2���data2
sigma = zeros(NofCH,nz);%�V�O�}1���data1�A2���data2
V = zeros(NofCH/4,nz);%1���data1�EVz(km/s)�A2���data1�EVr(km/s)�A3���data2�EVz(km/s)�A4���data2�EVr(km/s)
absV = zeros(NofCH/4,nz);%1���data1�EV(km/s)�A2���data2�EV(km/s)
T = zeros(NofCH,nz);%1���data1�E���x(eV)�A2���data1�E���x(eV)
Ti = zeros(NofCH/4,nz);%1���data1�E���x(eV)�A2���data1�E���ω��x(eV)
offset = zeros(NofCH/4,nz);

temp = zeros(NofCH,1);

% data1�̃K�E�X�t�B�b�e�B���O
if plot_fitting
    figure('Position',[300 50 1400 1000])
end
for k = 1:NofCH/4
    if plot_fitting
%         if mod(k,8) == 1%96�����̎�
%             figure('Position',[300 50 1400 1000])
%         end
    end
    %0�x�y�A�X�y�N�g������I�t�Z�b�g�����o
    for i = 1:2
        S1 = [L1(:,(k-1)*4+i) spectrum1(:,(k-1)*4+i)] %(k-1)*4+i�Ԗڂ̃`�����l����[�g��,���x]
        s1 = size(S1);%S1��[�s��,��]
        MAX1 = max(spectrum1(:,(k-1)*4+i)); %�X�y�N�g���̍ő�l
        j = 1;
        while j < s1(1)+1 %SN�̈����f�[�^������
            if S1(j,2) < MAX1*SN
                S1(j,:) = [];
            else
                j = j+1;
            end
            s1 = size(S1);
        end
        f = fit(S1(:,1),S1(:,2),'gauss1');
         
        %���[�����g�Œ��S�o��
        INT_spect_lamda = 0;
        INT_spect = 0;
        CENTER = 0;
        for j = 1:s1(1,1)
            INT_spect_lamda = INT_spect_lamda + S1(j,1)*S1(j,2);
            INT_spect = INT_spect + S1(j,2);
        end
%         INT_spect_lamda
%         INT_spect
        CENTER = INT_spect_lamda ./ INT_spect;%�Z���^�[�̈ʒu���o�Ă�
        
        
        coef=coeffvalues(f);
        amp((k-1)*4+i,1) = coef(1);
        shift((k-1)*4+i,1) = coef(2)-lambda0;%�ʏ펞�̃V�t�g
%         shift((k-1)*4+i,1) = CENTER-lambda0;%���[�����g�̎��̃V�t�g
        sigma((k-1)*4+i,1) = sqrt(coef(3)^2-(center((k-1)*4+i,5)*px2nm((k-1)*4+i,1))^2);
        T((k-1)*4+i,1) = 1.69e8*A*(2*sigma((k-1)*4+i,1)*sqrt(2*log(2))/lambda0)^2;
    end
    offset(k,1) = (shift((k-1)*4+1,1) + shift((k-1)*4+2,1))/2; %�Ό��������瓾��ꂽ�I�t�Z�b�g[nm]
%     offset(k,1) = 0; 
    %�I�t�Z�b�g�������ăK�E�X�t�B�b�e�B���O
    for i = 1:4
        S1 = [L1(:,(k-1)*4+i)-offset(k,1) spectrum1(:,(k-1)*4+i)]; %(k-1)*4+i�Ԗڂ̃`�����l����[�g��,���x]
        s1 = size(S1);%S1��[�s��,��]
        MAX1 = max(spectrum1(:,(k-1)*4+i)); %�X�y�N�g���̍ő�l
        j = 1;
        while j < s1(1)+1 %SN�̈����f�[�^������
            if S1(j,2) < MAX1*SN
                S1(j,:) = [];
            else
                j = j+1;
            end
            s1 = size(S1);
        end
%         S1(:,1)
%         S1(:,2)
        f = fit(S1(:,1),S1(:,2),'gauss1');
        INT_spect_lamda = 0;
        INT_spect = 0;
        CENTER = 0;
        for j = 1:s1(1,1)
            INT_spect_lamda = INT_spect_lamda + S1(j,1)*S1(j,2);
            INT_spect = INT_spect + S1(j,2);
        end
        CENTER = INT_spect_lamda ./ INT_spect;%�Z���^�[�̈ʒu���o�Ă�
        
        coef=coeffvalues(f);
        amp((k-1)*4+i,1) = coef(1);
        shift((k-1)*4+i,1) = coef(2)-lambda0;%�ʏ펞�̃V�t�g
%         shift((k-1)*4+i,1) = CENTER-lambda0;%���[�����g�̎��̃V�t�g
        sigma((k-1)*4+i,1) = sqrt(coef(3)^2-(center((k-1)*4+i,5)*px2nm((k-1)*4+i,1))^2);
        T((k-1)*4+i,1) = 1.69e8*A*(2*sigma((k-1)*4+i,1)*sqrt(2*log(2))/lambda0)^2;
        temp((k-1)*4+i) = T((k-1)*4+i,1);%���ׂẴ`�����l���̉��x��ۑ�
                    
        if plot_fitting
%             subplot(NofCH/8/3,8,mod((k-1),8)*4+i);%96�̎�
            subplot(NofCH/8,8,(k-1)*4+i);
            plot(f,S1(:,1),S1(:,2));
            xline(lambda0);
            title([num2str(T((k-1)*4+i,1)),' eV'])
            legend('off')
            xlabel('Wavelength [nm]')
            ylabel('Intensity [cnt]')
        end
    end
end
if plot_fitting
    sgtitle('Fitting data1 (Horizontal�FChannel Vertical�FPosition)')
end

%���x�ω��̊m�F�p�t�@�C������
save(['magflow/mat/temp_',num2str(t),'us_',num2str(ndata),'.mat'],'temp')


%data1�̗����A���x���v�Z
if cal_flow
    for i = 1:NofCH/4
        set = (i-1)*4;
        Va = -shift(set+4,1)/lambda0*Vc;
        Vb = -shift(set+3,1)/lambda0*Vc;
%         Va = -SHIFT(set+4,1)/lambda0*Vc;
%         Vb = -SHIFT(set+3,1)/lambda0*Vc;
        V(i,1) = (Va-Vb)/(2*cos(Angle*pi/180));%Vz
        V(i,2) = -(Va+Vb)/(2*sin(Angle*pi/180));%Vr
        absV(i,1) = sqrt(V(i,1)^2 + V(i,2)^2);
        %        fprintf('(z, r) = (%.2f, %.2f) cm\n (Vz, Vr) = (%.2f, %.2f) km/s\n'...
        %             ,z(i,1),r(i,1),V(i,1),V(i,2));
        %             fprintf('Ti: %.2f, %.2f, %.2f, %.2f eV\n'...
        %                     ,T(set+1,1),T(set+2,1),T(set+3,1),T(set+4,1));
        %         Ti(i,1) = (T(set+1,1)+T(set+2,1)+T(set+3,1)+T(set+4,1))/4;
        Ti(i,1) = trimmean(T(set+1:set+4,1),50);
%         fprintf('%.2f\n',Ti(i,1));
    end
end

%data1�̌v�Z���ʂ�ۑ�
if save_flow
    save(['magflow/mat/save_',num2str(t),'us_',num2str(ndata),'.mat'],'V','Ti')
end

%data2�̃K�E�X�t�B�b�e�B���O
if nz == 2
    if plot_fitting
        figure('Position',[300 50 1000 1000])
    end
    for k = 1:NofCH/4
        %0�x�y�A�X�y�N�g������I�t�Z�b�g�����o
        for i = 1:2
            S2 = [L1(:,(k-1)*4+i) spectrum2(:,(k-1)*4+i)]; %(k-1)*4+i�Ԗڂ̃`�����l����[�g��,���x]
            s2 = size(S1);%S1��[�s��,��]
            MAX2 = max(spectrum2(:,(k-1)*4+i)); %�X�y�N�g���̍ő�l
            j = 1;
            while j < s2(1)+1 %SN�̈����f�[�^������
                if S2(j,2) < MAX2*0.5
                    S2(j,:) = [];
                else
                    j = j+1;
                end
                s2 = size(S2);
            end
            f = fit(S2(:,1),S2(:,2),'gauss1');
            coef=coeffvalues(f);
            amp((k-1)*4+i,2) = coef(1);
            shift((k-1)*4+i,2) = coef(2)-lambda0;
            sigma((k-1)*4+i,2) = sqrt(coef(3)^2-(center((k-1)*4+i,5)*px2nm((k-1)*4+i,1))^2);
            T((k-1)*4+i,2) = 1.69e8*A*(2*sigma((k-1)*4+i,2)*sqrt(2*log(2))/lambda0)^2;
        end
        offset(k,2) = (shift((k-1)*4+1,2) + shift((k-1)*4+2,2))/2;%�Ό��������瓾��ꂽ�I�t�Z�b�g[nm]
        %�I�t�Z�b�g�������ăK�E�X�t�B�b�e�B���O
        for i = 1:4
            S2 = [L1(:,(k-1)*4+i) spectrum2(:,(k-1)*4+i)]; %(k-1)*4+i�Ԗڂ̃`�����l����[�g��,���x]
            s2 = size(S1);%S1��[�s��,��]
            MAX2 = max(spectrum2(:,(k-1)*4+i)); %�X�y�N�g���̍ő�l
            j = 1;
            while j < s2(1)+1 %SN�̈����f�[�^������
                if S2(j,2) < MAX2*0.5
                    S2(j,:) = [];
                else
                    j = j+1;
                end
                s2 = size(S2);
            end
            f = fit(S2(:,1),S2(:,2),'gauss1');
            coef=coeffvalues(f);
            amp((k-1)*4+i,2) = coef(1);
            shift((k-1)*4+i,2) = coef(2)-lambda0;
            sigma((k-1)*4+i,2) = sqrt(coef(3)^2-(center((k-1)*4+i,5)*px2nm((k-1)*4+i,1))^2);
            T((k-1)*4+i,2) = 1.69e8*A*(2*sigma((k-1)*4+i,2)*sqrt(2*log(2))/lambda0)^2;
            if plot_fitting
                subplot(NofCH/4,4,(k-1)*4+i);
                plot(f,S2(:,1),S2(:,2));
                xline(lambda0);
                title([num2str(T((k-1)*4+i,2)),' eV'])
                legend('off')
                xlabel('Wavelength [nm]')
                ylabel('Intensity [cnt]')
            end
        end
    end
    if plot_fitting
        sgtitle('Fitting data1 (Horizontal�FChannel Vertical�FPosition)')
    end
    
    %data2�̗����A���x���v�Z
    if cal_flow
        for i = 1:NofCH/4
            set = (i-1)*4;
            Va = -shift(set+4,2)/lambda0*Vc;
            Vb = -shift(set+3,2)/lambda0*Vc;
            V(i,3) = (Va-Vb)/(2*cos(Angle*pi/180));%Vz
            V(i,4) = -(Va+Vb)/(2*sin(Angle*pi/180));%Vr
            absV(i,2) = sqrt(V(i,3)^2 + V(i,4)^2);
            %        fprintf('(z, r) = (%.2f, %.2f) cm\n (Vz, Vr) = (%.2f, %.2f) km/s\n'...
            %             ,z(i,2),r(i,2),V(i,3),V(i,4));
            %             fprintf('Ti: %.2f, %.2f, %.2f, %.2f eV\n'...
            %                     ,T(set+1,2),T(set+2,2),T(set+3,2),T(set+4,2));
            Ti1=(T(set+1,2)+T(set+2,2)+T(set+3,2)+T(set+4,2))/4;
            fprintf('%.2f\n',Ti1);
        end
    end
end

%�����A���x���v���b�g
if cal_flow
    if plot_flow
        figure('Position',[600 150 300 600])
        Ticon = repmat(Ti(:,1),1,5);%�������}�p���σC�I�����x
        zcon = z_measured(1,1);
        s = pcolor([zcon-2 zcon-1 zcon zcon+1 zcon+2],r_measured,Ticon);
        s.FaceColor = 'interp';
        s.EdgeAlpha = 0;
        colormap('jet')
        colorbar
        hold on
        plot(z_measured,r_measured,'xr');
        hold on
        if nz == 1
            q = quiver(z_measured,r_measured,V(:,1)*factor,V(:,2)*factor);
        end
        if nz == 2
            q = quiver(z_measured,r_measured,[V(:,1),V(:,3)]*factor,[V(:,2),V(:,4)]*factor);
        end
        q.LineWidth = 2;
        q.MaxHeadSize = 10;
        q.AutoScale = 'off';
        q.Color = 'k';
        xlim([min(z_measured,[],'all')-2.5 max(z_measured,[],'all')+2.5])
        ylim([min(r_measured,[],'all')-2.5 max(r_measured,[],'all')+2.5])
        title('Ion Flow [km/s]','Color','black','FontWeight','bold')
        absV = round(absV,1);
        for j = 1:nz
            for i = 1:NofCH/4
                txt = text(z_measured(i,j)+0.5,r_measured(i,j)-0.2,num2str(absV(i,j)));
                txt.FontSize = 18;
                txt.Color = 'k';
                txt.FontWeight = 'bold';
            end
        end
        xlabel('Z [cm]')
        ylabel('R [cm]')
        ax = gca;
        ax.FontSize = 14;
        grid on
        daspect([1 1 1])
        hold off
    end
end
