function [coeffvals,ci,spectra_deconv] = deconvolution(lambda,I,sigma_instru)
stepsize = abs(lambda(1)-lambda(2));
x = linspace(-abs(lambda(1)-lambda(end))/2,abs(lambda(1)-lambda(end))/2,length(I));
h = exp(-x.^2/(2*sigma_instru^2))/((2*pi)^0.5*sigma_instru);
[v,~] = deconv(I,h,"same",Method="least-squares",RegularizationFactor=0.1);
f = fit(lambda,v'/stepsize,'gauss1');
coeffvals= coeffvalues(f)/sqrt(2);
ci = confint(f,0.68)/sqrt(2);
spectra_deconv = f(lambda);

% figure(1);
% subplot(2,1,1)
% hold on;
% plot(lambda,v/stepsize);
% plot(lambda,I);
% hold off;
% legend(["deconv","raw"])
% subplot(2,1,2)
% plot(x'+lambda(i_max),h);
% legend("transfer")
end