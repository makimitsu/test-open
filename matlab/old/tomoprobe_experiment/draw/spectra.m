figure('Position',[300 550 1000 600])
for i = 1:numCh
    S3 = [ShiftLambda(:,i) ArrangedSpectrum(:,i)]; %i�Ԗڂ̃`�����l����[X,Y]
    subplot(numCh/2,2,i);
    plot(S3(:,1),S3(:,2));
    xline(0);
    legend('off')
    title(['CH',num2str(i)])
    xlabel('Shift wavelength [nm]')
    ylabel('Intensity [a.u.]')
end
sgt = sgtitle('Measured spectra');
sgt.FontSize = 20;
