function [] = plot_inversion_compare(IDSPdata,IDSP)
n_Theta = numel(IDSPdata.Angle);
[n_L1,~] = size(IDSPdata.Lambda);
% figure('Position',[500 150 400 700])
figure('Position',[200 150 1000 900])
for k = 1:IDSP.n_r
    ReP = IDSPdata.W(:,:,k)*IDSPdata.F(:,k);
    draw_P = reshape(IDSPdata.P(:,k), [n_L1,n_Theta]);%描画用に整形
    draw_ReP = reshape(ReP, [n_L1,n_Theta]);%描画用に整形
    for i = 1:n_Theta
        % subplot(n_Theta,1,i);
        subplot(IDSP.n_r,n_Theta,n_Theta*(k-1)+i);
        switch i
            case 1
                plot(IDSPdata.Lambda(:,4*(k-1)+2), draw_P(:,i),'r', IDSPdata.Lambda(:,4*(k-1)+2), draw_ReP(:,i),'b')
            case 2
                plot(IDSPdata.Lambda(:,4*(k-1)+4), draw_P(:,i),'r', IDSPdata.Lambda(:,4*(k-1)+4), draw_ReP(:,i),'b')
            case 3
                plot(IDSPdata.Lambda(:,4*(k-1)+3), draw_P(:,i),'r', IDSPdata.Lambda(:,4*(k-1)+3), draw_ReP(:,i),'b')
        end
        legend({'measured','recon.'},'Location', 'northwest')
        title(['r = ',num2str(IDSP.r(k)),'[cm], θ = ',num2str(IDSPdata.Angle(i)), '°'])
        xlabel('Shift wavelength [nm]');
        ylabel('Intensity [a.u.]');
        % xlim([-0.04 0.04])
    end
    sgt = sgtitle(['Spectra comparison (Horizontal：View Line, Vertical：Measured Position)',newline, ...
        'shot',num2str(IDSP.shot),'-',num2str(IDSP.delay),'us-w=',num2str(IDSP.width),'-gain=',num2str(IDSP.gain),'.asc']);
    sgt.FontSize = 16;
end
end

