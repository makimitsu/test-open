date = 230826;%230523;
shot = 3;%8;
ch = [7 8];
calib(ch) = [63.9568,223.2319];

rgw2txt(date,shot)
rogowski_plot(date,shot,ch,calib)