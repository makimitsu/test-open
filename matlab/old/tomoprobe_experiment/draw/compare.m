ReP = W*ReF;
drP = reshape(P, [numLambda,numTheta]);%描画用に整形
drReP = reshape(ReP, [numLambda,numTheta]);%描画用に整形
figure('Position',[500 150 400 700])
for k = 1:numTheta
    subplot(numTheta,1,k);
    plot(ShiftLambda(:,k), drP(:,k),'r', ShiftLambda(:,k), drReP(:,k),'b')
    legend({'measured','recon.'},'Location', 'northwest')
    title(['θ = ',num2str(Angle(k)), '°'])
    xlabel('Shift wavelength [nm]');
    ylabel('Intensity [a.u.]');
end
sgt = sgtitle('Spectra comparison');
sgt.FontSize = 20;
