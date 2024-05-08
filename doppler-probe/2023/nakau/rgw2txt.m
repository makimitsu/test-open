function [] = rgw2txt(date)

% Get all rgw files in the current folder
current_folder = strcat('/Users/itsuki/Documents/東大/研究資料/Doppler_probe_data/');

files = dir(strcat(current_folder,num2str(date),'/','*.rgw'));
% Loop through each
for id = 1:length(files)
    % Get the file name (minus the extension)
    [path,name,~] = fileparts(files(id).name);
    rename = strcat(current_folder,name,'.txt');
    old_name = strcat(current_folder,files(id).name);
    copyfile(old_name,rename);
end

end