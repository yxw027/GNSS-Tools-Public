function [epochs, clk, clk_sig] = Read_GPS_05sec_CLK(filename, flag32)
% Read_GPS_05sec_CLK
% Reads high rate (5 second) CODE time data from RINEX
if nargin < 2
    flag32 = 0;
end
temp = importdata(filename);
textData = temp.textdata;

% Pull PRN List
inds = [strfind(textData(:,1),'PRN LIST')];
inds = find(~cellfun(@isempty,inds));
prns = [];
for i = 1:length(inds)
    linetext = textData{inds(i)};
    gInds = strfind(linetext,'G');
    prns = [prns; arrayfun(@(x) str2double(linetext(x+(1:2))), gInds)'];
end

% Ensure there are PRNs available
if isempty(prns)
    fprintf(2, 'Error no GPS PRNs identified in file: %d\n', length(idx));
    epochs = [];
    clk = [];
    clk_sig = [];
    return
end

% Initialize output variables
epochs = (0:5:86395)';
clk = NaN(max(prns),length(epochs));
if flag32
    clk = NaN(32,length(epochs));
end
clk_sig = clk;

indAS = find(strcmp(textData(:,1),'AS'));
indOffset = size(temp.textdata,1)-size(temp.data,1);
for i = 1:length(prns)
    indG = intersect(find(strcmp(textData(:,2),['G' num2str(prns(i), '%02d')])),indAS)-indOffset;
    
    year = temp.data(indG,1);
    month = temp.data(indG,2);
    day = temp.data(indG,3);
    hours = temp.data(indG,4);
    minutes = temp.data(indG,5);
    seconds = temp.data(indG,6);
    
    jdx = round(hours*720 + minutes*12 + seconds/5) + 1;
    clk(prns(i),jdx) = temp.data(indG,8);
    k = find(mod(minutes,5) == 0 & mod(seconds, 60) ==0);
    if size(temp.data,2) >= 9
        clk_sig(prns(i),jdx(k)) = temp.data(indG(k),9);
    end
end


end


% Below is the old version of the above code:


% inds = find([strcmp([textData(:,1)],'AS')]);

% filetext = fileread(filename);
% idx = strfind(filetext, 'PRN LIST');
% prns = [];
% for i = 1:length(idx)
%     j=idx(i)-60;
%     while j < idx(i) && filetext(j) == 'G'
%         prns = [prns; str2double(filetext((j+1):(j+2)))];
%         j = j+4;
%     end
% end
% if isempty(prns)
%     fprintf(2, 'Error no GPS PRNs identified in file: %d\n', length(idx));
%     epochs = [];
%     clk = [];
%     clk_sig = [];
% 	return
% end
% 
% epochs = (0:5:86395)';
% clk = NaN(max(prns),length(epochs));
% clk_sig = clk;
% 
% for i = 1:length(prns)
%     text_str = ['AS G' num2str(prns(i), '%02d')];
%     idx = strfind(filetext, text_str);
%     if isempty(idx)
% %         fprintf(2, 'Error no data for GPS PRN: %02d\n', prns(i));
% %         epochs = [];
% %         clk = [];
% %         clk_sig = [];
% %         return
%         continue
%     end
% %     tic
%     year = arrayfun(@(x) str2double(filetext(x+(8:11))), idx);
%     if length(unique(year)) > 1
%         fprintf(2, 'Error more than one year in the daily file\n');
%         epochs = [];
%         clk = [];
%         clk_sig = [];
%         return
%     end
%     month = arrayfun(@(x) str2double(filetext(x+(13:14))), idx);
%     if length(unique(month)) > 1
%         fprintf(2, 'Error more than one month in the daily file\n');
%         epochs = [];
%         clk = [];
%         clk_sig = [];
%         return
%     end    
%     day = arrayfun(@(x) str2double(filetext(x+(16:17))), idx);
%     if length(unique(day)) > 1
%         fprintf(2, 'Error more than one day in the daily  file\n');
%         epochs = [];
%         clk = [];
%         clk_sig = [];
%         return
%     end
%     hours = arrayfun(@(x) str2double(filetext(x+(19:20))), idx);
%     minutes = arrayfun(@(x) str2double(filetext(x+(22:23))), idx);
%     seconds = arrayfun(@(x) str2double(filetext(x+(25:33))), idx);
%     jdx = round(hours*720 + minutes*12 + seconds/5) + 1;
%     
%     clk(prns(i),jdx) = arrayfun(@(x) str2double(filetext(x+(40:59))), idx);
%     k = find(mod(minutes,5) == 0 & mod(seconds, 60) ==0);
%     clk_sig(prns(i),jdx(k)) = arrayfun(@(x) str2double(filetext(x+(61:78))), idx(k));
%     
%     toc
% end