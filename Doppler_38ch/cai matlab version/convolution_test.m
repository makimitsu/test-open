clear all
stepsize = 0.001;
x_signal = 480:stepsize:481;
x_trans = -0.5:stepsize:0.5;
% the input function
a1 = 4000; b1=480.5; c1=.05;
func_real = a1*exp(-((x_signal-b1).^2)/(2*c1^2))/((2*pi)^0.5*c1);

% the transfer function
a2 = 1; b2=0; c2=.1;
func_trans = a2*exp(-(x_trans-b2).^2/(2*c2^2))/((2*pi)^0.5*c2);

% calculate the output function using conv
func_conv = conv(func_trans,func_real,"same")*stepsize;

% deconvolution

stepsize_interp = 0.001;
x_signal_interp = 480:stepsize_interp:481;
func_conv_interp = interp1(x_signal,func_conv,x_signal_interp);
x_trans_interp = -0.5:stepsize_interp:0.5;
func_trans_interp = a2*exp(-(x_trans_interp-b2).^2/(2*c2^2))/((2*pi)^0.5*c2);

[v,~] = deconv(func_conv_interp,func_trans_interp,"same",Method="least-squares",RegularizationFactor=0);
figure(1),clf(1)
subplot(2,1,1),hold on
yyaxis left
plot(x_signal, func_real,'b','LineStyle','-');
plot(x_signal_interp, func_conv_interp,'r','LineStyle','-');
ylim([-5*a1 10*a1]),xlim([480 481]);
yyaxis right
plot(x_trans_interp+480.5, func_trans_interp,'k','LineStyle','-');
legend('input func', 'convoluted func', 'transfer func');
ylim([0 max(func_trans)]),xlim([480 481]),hold off
subplot(2,1,2)
plot(x_signal_interp,v/stepsize_interp,'b'),legend('deconvoluted func'),ylim([-5*a1 10*a1]),xlim([480 481]);

% alpha = 0.1:0.1:100;
% error_deconv = zeros(size(alpha));
% for i = 1:length(alpha)
%     [v,~] = deconv(y,h,"same",Method="least-squares",RegularizationFactor=alpha(i));
%     error_deconv(i) = sum((v/stepsize-u).^2);
% end
% figure,loglog(alpha,error_deconv);