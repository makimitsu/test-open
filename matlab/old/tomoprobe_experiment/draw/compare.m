ReP = W*ReF;
drP = reshape(P, [numLambda,numTheta]);%�`��p�ɐ��`
drReP = reshape(ReP, [numLambda,numTheta]);%�`��p�ɐ��`
figure('Position',[500 150 400 700])
for k = 1:numTheta
    subplot(numTheta,1,k);
    plot(ShiftLambda(:,k), drP(:,k),'r', ShiftLambda(:,k), drReP(:,k),'b')
    legend({'measured','recon.'},'Location', 'northwest')
    title(['�� = ',num2str(Angle(k)), '��'])
    xlabel('Shift wavelength [nm]');
    ylabel('Intensity [a.u.]');
end
sgt = sgtitle('Spectra comparison');
sgt.FontSize = 20;
