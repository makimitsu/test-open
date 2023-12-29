function [] = ProcessTIFImageWithNLMFilter(TIFfolder,Shotlist)
% input1: String 'G:\My Drive\X-ray\Data\TIF\231213'
% input2: Array (2:10)

if exist(TIFfolder,'dir') == 0
    disp('No directory found');
    return
elseif exist(strcat(TIFfolder,'/filtered'),'dir') == 0
    mkdir(strcat(TIFfolder,'/filtered'));
end

for i = Shotlist
    rawTIFPath = strcat(TIFfolder,'/shot',num2str(i,'%03i'),'.tif');
    filteredTIFPath = strcat(TIFfolder,'/filtered/shot',num2str(i,'%03i'));
    rawTIF = imread(rawTIFPath);
    patch = rawTIF(1:200,1:200);
    patchSq = patch.^2;
    edist = sqrt(sum(patchSq,3));
    patchSigma = sqrt(var(edist(:)));
    [filteredTIF,~] = imnlmfilt(rawTIF,'SearchWindowSize',91,'ComparisonWindowSize',13,'DegreeOfSmoothing',patchSigma*0.625);
    save(filteredTIFPath,"filteredTIF");disp(strcat('Shot ',num2str(i),' completed.'));
end

end