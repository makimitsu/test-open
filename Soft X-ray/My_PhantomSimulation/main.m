setenv('SXR_IMAGE_DIR','G:\My Drive\X-ray\Data\TIF');
setenv('SXR_RECONSTRUCTED_DIR','G:\My Drive\X-ray\Data\SXROUT');
setenv('SXR_MATRIX_DIR','G:\My Drive\X-ray\Data\SXROUT');
% --
% load('parameters.mat');
load('G:\My Drive\X-ray\Data\SXROUT\LF_LR\230920\shot15\6.mat');
EE = EE3;
% --

prompt = {'doSave:','doFilter:','doNLR:'};
definput = {'','',''};
if exist('doSave','var')
    definput{1} = num2str(doSave);
end
if exist('doFilter','var')
    definput{2} = num2str(doFilter);
end
if exist('doNLR','var')
    definput{3} = num2str(doNLR);
end
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims,definput);
doSave = logical(str2double(cell2mat(answer(1))));
doFilter = logical(str2double(cell2mat(answer(2))));
doNLR = logical(str2double(cell2mat(answer(3))));
plot_sxr_multi(doSave,doFilter,doNLR,M,K,gm2d2,U2,s2,v2,EE,range);