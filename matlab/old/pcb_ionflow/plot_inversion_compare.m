function [] = plot_inversion_compare(F,W,P,L1_trans,mpoints,Angle,ICCD)
n_Theta = numel(Angle);
[n_L1,~] = size(L1_trans);
% figure('Position',[500 150 400 700])
figure('Position',[200 150 1000 900])
for k = 1:mpoints.n_r
    ReP = W(:,:,k)*F(:,k);
    draw_P = reshape(P(:,k), [n_L1,n_Theta]);%描画用に整形
    draw_ReP = reshape(ReP, [n_L1,n_Theta]);%描画用に整形
    for i = 1:n_Theta
        % subplot(n_Theta,1,i);
        subplot(mpoints.n_r,n_Theta,n_Theta*(k-1)+i);
        switch i
            case 1
                plot(L1_trans(:,4*(k-1)+2), draw_P(:,i),'r', L1_trans(:,4*(k-1)+2), draw_ReP(:,i),'b')
            case 2
                plot(L1_trans(:,4*(k-1)+4), draw_P(:,i),'r', L1_trans(:,4*(k-1)+4), draw_ReP(:,i),'b')
            case 3
                plot(L1_trans(:,4*(k-1)+3), draw_P(:,i),'r', L1_trans(:,4*(k-1)+3), draw_ReP(:,i),'b')
        end
        legend({'measured','recon.'},'Location', 'northwest')
        title(['r = ',num2str(mpoints.r(k)),'[cm], θ = ',num2str(Angle(i)), '°'])
        xlabel('Shift wavelength [nm]');
        ylabel('Intensity [a.u.]');
        % xlim([-0.04 0.04])
    end
    sgt = sgtitle(['Spectra comparison (Horizontal：View Line, Vertical：Measured Position)',newline, ...
        'shot',num2str(ICCD.shot),'-',num2str(ICCD.trg),'us-w=',num2str(ICCD.exp_w),'-gain=',num2str(ICCD.gain),'.asc']);
    sgt.FontSize = 16;
end
end

