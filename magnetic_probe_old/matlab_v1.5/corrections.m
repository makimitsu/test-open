function data = corrections(data,date)
% This function assigns values to dead channels
% input: 
%   raw data
%   date of experiment(yymmdd)
% output:
%   data corrected

if date == 191209
data(39,:) = (data(38,:) + data(40,:))/2;
data(93,:) = (data(92,:) + data(94,:))/2;
data(120,:) = (data(119,:) + data(121,:))/2;
data(185,:) = (data(176,:) + data(186,:))/2;
data(205,:) = (data(204,:) + data(206,:))/2;
data(239,:) = 2*data(110,:) - data(109,:);
data(240,:) = 2*data(239,:) - data(110,:);

data(225:236,:) = data(99:110,:);
data(237:238,:) = data(239:240,:);

elseif date >= 191210 && date < 191218 
data(39,:) = (data(38,:) + data(40,:))/2;
data(93,:) = (data(92,:) + data(94,:))/2;
data(120,:) = (data(119,:) + data(121,:))/2;
data(185,:) = (data(176,:) + data(186,:))/2;
data(205,:) = (data(204,:) + data(206,:))/2;
data(239,:) = 2*data(110,:) - data(109,:);
data(240,:) = 2*data(239,:) - data(110,:);

data(185:188,:) = data(62:65,:);
data(225:236,:) = data(99:110,:);
data(237:238,:) = data(239:240,:);

elseif date >= 191218 && date < 200121
data(39,:) = (data(38,:) + data(40,:))/2;
data(42,:) = (data(28,:) + data(98,:))/2;
data(93,:) = (data(92,:) + data(94,:))/2;
data(120,:) = (data(119,:) + data(121,:))/2;
data(205,:) = (data(204,:) + data(206,:))/2;
data(232,:) = (data(84,:) + data(98,:))/2;
data(234,:) = (data(84,:) + data(56,:))/2;
data(236,:) = 2*data(56,:) - data(234,:);
data(70,:) = 2*data(236,:) - data(56,:);
data(239,:) = 2*data(110,:) - data(109,:);

elseif date >= 200121 && date < 200128
data(93,:) = (data(92,:) + data(94,:))/2;
data(219,:) = (data(110,:) + data(244,:))/2;
data(243,:) = 2*data(242,:) - data(259,:);
data(129,:) = (data(130,:) + data(128,:))/2;

elseif date >= 200128 && date < 200121
data(93,:) = (data(92,:) + data(94,:))/2;
data(219,:) = (data(110,:) + data(244,:))/2;
data(243,:) = 2*data(242,:) - data(259,:);

elseif date >= 200130 && date < 200414
data(93,:) = (data(92,:) + data(94,:))/2;
data(219,:) = (data(110,:) + data(244,:))/2;
data(243,:) = 2*data(242,:) - data(259,:);
data(173,:) = (data(174,:) + data(172,:))/2;
data(266,:) = (data(265,:) + data(267,:))/2;

elseif date >= 200414
data(58,:) = (data(57,:) + data(59,:))/2;
data(93,:) = (data(92,:) + data(94,:))/2;
data(219,:) = (data(110,:) + data(244,:))/2;
data(243,:) = 2*data(242,:) - data(259,:);
data(259,:) = (data(242,:) + data(240,:))/2;
data(174,:) = (data(173,:) + data(175,:))/2;
data(266,:) = (data(265,:) + data(267,:))/2;

elseif date >= 200416
data(58,:) = (data(57,:) + data(59,:))/2;
data(93,:) = (data(92,:) + data(94,:))/2;
data(219,:) = (data(110,:) + data(244,:))/2;
data(243,:) = 2*data(242,:) - data(259,:);
data(259,:) = (data(242,:) + data(240,:))/2;
data(174,:) = (data(173,:) + data(175,:))/2;
data(266,:) = (data(265,:) + data(267,:))/2;

else
    disp('Wrong date!!!');
    return
end