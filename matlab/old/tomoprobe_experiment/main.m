function [Vmean_z,Vmean_r,Ti_z,Ti_r] = main(t,nz,nr,ndata,inversion_method,draw_spectra,draw_result,draw_type,draw_analisis,draw_compare)
%����,�v���_Z�����ԍ�,�v���_R�����ԍ�,�f�[�^�ԍ�,�t�ϊ�(1~5),�X�y�N�g���`��,���ʕ`��,������(1)/�O����(2),���͕`��,��r�`��
tic
time = t;
r_measured = (100 + 25*nr);
%�p�����[�^���`
run define/parameter.m

%�X�y�N�g���f�[�^P����荞��
run data/dataP.m

%�ϑ��X�y�N�g����`��
if draw_spectra
    run draw/spectra.m
end

%W(F->P)���쐬
run make/funcW.m

%---------------Invertion theory----------------

%pinv���g����W^(-1)�����߂�B
if inversion_method == 1
    run method/pseudo.m
    run method/goto0.m
end

%Tikhonov 0th
if inversion_method == 2
    run method/Tikhonov0.m
    run method/goto0.m
end

%Tikhonov 1st
if inversion_method == 3
    run method/Tikhonov1.m
    run method/goto0.m
end

%Tikhonov 2nd
if inversion_method == 4
    run method/Tikhonov2.m
    run method/goto0.m
end

%�񕉐���SIRT
if inversion_method == 5
    run method/SIRT.m
end

%minimum Fischer information
if inversion_method == 6
    run method/MFI.m
    run method/goto0.m
end

%���ʂ̕\��
if inversion_method ~= false
    drReF = reshape(ReF, [numVx,numVy]);%�`��p�ɐ��`
    drReF = flipud(drReF);
    [Summity,Summitx] = find(drReF == max(drReF(:)));
    
    %�č\�����ʂ�`��
    if draw_result
        run draw/result.m
    end
    
    %�č\�����ʂ��瓾����X�y�N�g�����v�Z�A�X�y�N�g���f�[�^�Ɣ�r
    if draw_compare
        run draw/compare.m
    end
    %�������v�Z
    Fsum = 0;
    Vmean_z = 0;
    Vmean_r = 0;
    for i = 1:numVy
        for j = 1:numVx
            if drReF(i,j)> drReF(Summity,Summitx)*0.5
                Vmean_z = Vmean_z + drReF(i,j)*Vx(j);
                Vmean_r = Vmean_r + drReF(i,j)*Vy(i);
                Fsum = Fsum + drReF(i,j);
            end
        end
    end
    Vmean_z = Vmean_z/Fsum;
    Vmean_r = Vmean_r/Fsum;
    Ti_z = 0;
    Ti_r = 0;
    %     fprintf('(Vz = %.1f, Vr = %.1f)[km/s].\n',Vmean_z,Vmean_r);
    %     %�[�����x���v�Z(���U����v�Z)
    %     Vdisp = zeros(2,1);
    %     for i = 1:numVy
    %         for j = 1:numVx
    %             if drReF(i,j)> drReF(Summity,Summitx)*0
    %                 Vdisp(1,1) = Vdisp(1,1) + drReF(i,j)*(Vx(j)-Vmean_z)^2;
    %                 Vdisp(2,1) = Vdisp(2,1) + drReF(i,j)*(Vy(i)-Vmean_r)^2;
    %             end
    %         end
    %     end
    %
    %     Vdisp(1,1) = Vdisp(1,1)/Fsum - (mean(inst(1+4*(nr-1):4+4*(nr-1)),'all')*deltaLambda/lambda0*Vc)^2;
    %     Vdisp(2,1) = Vdisp(2,1)/Fsum - (mean(inst(1+4*(nr-1):4+4*(nr-1)),'all')*deltaLambda/lambda0*Vc)^2;
    % %     Vdisp(1,1) = Vdisp(1,1)/Fsum;%���u�֐��Ȃ�
    % %     Vdisp(2,1) = Vdisp(2,1)/Fsum;%���u�֐��Ȃ�
    %     Ti_z = Vdisp(1,1)*10^6*A*mp/(2*kB);
    %     Ti_r = Vdisp(2,1)*10^6*A*mp/(2*kB);
    %     fprintf('(Tz = %.1f, Tr = %.1f)[eV].\n',Ti_z,Ti_r);
    
%     %�[�����x���v�Z(���l������v�Z)
%     W_0 = zeros(numLambda,numVx*numVy);%0�x�X�y�N�g�����v�Z
%     for i = 1:numLambda
%         Unit_0 = [cos(0) sin(0)];
%         for j = 1:numVx*numVy
%             %W(:,j)�̑Ή����x�ԍ�x,y���擾�F�Ή����x��Vx(x),Vy(y)
%             x = 1 + idivide(j-1, int8(numVy), 'floor');
%             y = int8(mod(j,numVy));
%             if y == 0
%                 y = numLambda;
%             end
%             V = [Vx(x), Vy(y)];
%             D = dot(V, Unit_0);%�����������x
%             %      W(i,j) = 1 - abs(Lambda(r)/lambda0*Vc-D)*(1/(deltaLambda/lambda0*Vc));%�ꎟ�֐��t�B���^
%             W_0(i,j) = 1 - ((ShiftLambda(i,1)/lambda0*Vc-D)*(1/(deltaLambda/lambda0*Vc)))^2;%�񎟊֐��t�B���^
%             if W_0(i,j) < 0
%                 W_0(i,j) = 0;
%             end
%         end
%     end
%     P_0 = W_0*ReF;
%     SummitP_0 = find(P_0 == max(P_0));
%     %�v���b�g
% %     figure('Position',[300 550 600 400])
% %     plot(Vx,P_0/max(P_0), 'b','LineWidth',2)
% %     hold on
% %     yline(0.5, '--r','LineWidth',2);
% %     xlabel('V_{Z} [km/s]')
% %     ylabel('Relative intensity')
% %     xlim([-60 60])
% %     xticks(-60:20:60)
% %     yticks(0:0.1:1)
% %     ax = gca;
% %     ax.FontSize = 18;
% %     hold off
%     for i = 1:numLambda - SummitP_0
%         if P_0(i+SummitP_0) < max(P_0)*0.5
%             I_xmax = i+SummitP_0;
%             break
%         end
%     end
%     fwhm_xmax = Vx(I_xmax);%���lx�E�[
%     for i = 1:SummitP_0
%         if P_0(SummitP_0-i+1) < max(P_0)*0.5
%             I_xmin = SummitP_0-i+1;
%             break
%         end
%     end
%     fwhm_xmin = Vx(I_xmin);%���lx���[
%     fwhm_x = fwhm_xmax - fwhm_xmin;
%     sigma_x = fwhm_x/2.35;
% %     S_x = sigma_x^2 - (mean(inst(1+4*(nr-1):4+4*(nr-1)),'all')*deltaLambda/lambda0*Vc)^2;
%     S_x = sigma_x^2;%���u�֐��Ȃ�
%     Ti_z = S_x*10^6*A*mp/(2*kB);
%     
%     W_90 = zeros(numLambda,numVx*numVy);%0�x�X�y�N�g�����v�Z
%     for i = 1:numLambda
%         Unit_90 = [cos(pi/2) sin(pi/2)];
%         for j = 1:numVx*numVy
%             %W(:,j)�̑Ή����x�ԍ�x,y���擾�F�Ή����x��Vx(x),Vy(y)
%             x = 1 + idivide(j-1, int8(numVy), 'floor');
%             y = int8(mod(j,numVy));
%             if y == 0
%                 y = numLambda;
%             end
%             V = [Vx(x), Vy(y)];
%             D = dot(V, Unit_90);%�����������x
%             %      W(i,j) = 1 - abs(Lambda(r)/lambda0*Vc-D)*(1/(deltaLambda/lambda0*Vc));%�ꎟ�֐��t�B���^
%             W_90(i,j) = 1 - ((ShiftLambda(i,1)/lambda0*Vc-D)*(1/(deltaLambda/lambda0*Vc)))^2;%�񎟊֐��t�B���^
%             if W_90(i,j) < 0
%                 W_90(i,j) = 0;
%             end
%         end
%     end
%     P_90 = W_90*ReF;
%     SummitP_90 = find(P_90 == max(P_90));
%     %�v���b�g
% %     figure('Position',[300 550 600 400])
% %     plot(Vx,P_90/max(P_90), 'b','LineWidth',2)
% %     hold on
% %     yline(0.5, '--r','LineWidth',2);
% %     xlabel('V_{R} [km/s]')
% %     ylabel('Relative intensity')
% %     xlim([-60 60])
% %     xticks(-60:20:60)
% %     yticks(0:0.1:1)
% %     ax = gca;
% %     ax.FontSize = 18;
% %     hold off
%     for i = 1:numLambda - SummitP_90
%         if P_90(i+SummitP_90) < max(P_90)*0.5
%             I_ymax = i+SummitP_90;
%             break
%         end
%     end
%     fwhm_ymax = Vy(I_ymax);%���ly�E�[
%     for i = 1:SummitP_90
%         if P_0(SummitP_90-i+1) < max(P_90)*0.5
%             I_ymin = SummitP_90-i+1;
%             break
%         end
%     end
%     fwhm_ymin = Vy(I_ymin);%���lx���[
%     fwhm_y = fwhm_ymax - fwhm_ymin;
%     sigma_y = fwhm_y/2.35;
% %     S_y = sigma_y^2 - (mean(inst(1+4*(nr-1):4+4*(nr-1)),'all')*deltaLambda/lambda0*Vc)^2;
%     S_y = sigma_y^2;%���u�֐��Ȃ�
%     Ti_r = S_y*10^6*A*mp/(2*kB);
%         fprintf('Reconstructd: (Tz = %.1f, Tr = %.1f)[eV].\n',Ti_z,Ti_r);
end
toc
end
