function [] = main_multi_temp(inversion_method,draw_spectra,draw_result,draw_type,draw_analisis,draw_compare)
%�t�ϊ�(1~5),�X�y�N�g���`��,���ʕ`��,������(1)/�O����(2),���͕`��,��r�`��

nz = 1;%�v���_Z�����ԍ�(nz = 1:-2.1cm, nz = 2:+2.1cm)
NT = 6;%�v���g���K���Ԑ�(1~6)
NR = 7;%�v���_��(1~7)
NData = 5;%�f�[�^��(1~5)
Tz_t = zeros(NT,NR);%����t�̕���Tz
Tr_t = zeros(NT,NR);%����t�̕���Tr
Sz_t = zeros(NT,NR);%����t��Tz�̕W���΍�
Sr_t = zeros(NT,NR);%����t��Tr�̕W���΍�

% Vmean_z = zeros(NData,NR);
% Vmean_r = zeros(NData,NR);
% Tz = zeros(NData,NR);
% Tr = zeros(NData,NR);

% for t = 1:NT
%     for j = 1:NR
%         for k = 1:NData
%             [Vmean_z(k,j),Vmean_r(k,j),Tz(k,j),Tr(k,j)] = main(463+2*t,nz,j,k,inversion_method,draw_spectra,draw_result,draw_type,draw_analisis,draw_compare);
%         end
%     end
%     %�ۑ�
%     save(['mat/save_fwhm_plus_',num2str(463+2*t),'us.mat'],'Vmean_z','Vmean_r','Tz','Tr')
%     Tz_t(t,:) = mean(Tz);
%     Tr_t(t,:) = mean(Tr);
%     Sz_t(t,:) = std(Tz);
%     Sr_t(t,:) = std(Tr);
% end

% for t = 1:NT
%     %�ǂݍ���
%     load(['mat/save_fwhm_minus_',num2str(463+2*t),'us.mat'],'Tz','Tr')
%     Tz_t(t,:) = mean(Tz);
%     Tr_t(t,:) = mean(Tr);
%     Sz_t(t,:) = std(Tz);
%     Sr_t(t,:) = std(Tr);
% end
% %-------------------�[�����x�̈ꎟ�����z���v���b�g-------------------
% plot_temp(NT,NR,Tz_t,Tr_t,Sz_t,Sr_t)

for t = 1:NT
    %�ǂݍ���
    load(['mat/save_fwhm_minus_',num2str(463+2*t),'us.mat'],'Vmean_z','Vmean_r')
    Tz_t(t,:) = mean(Vmean_z);
    Tr_t(t,:) = mean(Vmean_r);
    Sz_t(t,:) = std(Vmean_z);
    Sr_t(t,:) = std(Vmean_r);
end
%-------------------���x�����̈ꎟ�����z���v���b�g-------------------
plot_velociy(NT,NR,Tz_t,Tr_t,Sz_t,Sr_t)


end
