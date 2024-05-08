date = 230310;
if exist(['fig/',num2str(date)],'dir') == 0
    mkdir(sprintf("fig/%s", num2str(date)));
end
