function [YearChange,dayChange] = ftp_file_general(YearList,dayList,ftpStruc,varargin)


ftpSite      = ['' ftpStruc.ftpSite];
sourceFormat = ftpStruc.sourceFormat;
destDir      = ftpStruc.destDir;
destFormat   = ftpStruc.destFormat;
fileFormat   = ftpStruc.fileFormat;
unzipFlag    = ftpStruc.unzipFlag;

% Years and days when files have been updated
YearChange = [];
dayChange  = [];

if nargin >= 4
    val1 = varargin{1};
end
if nargin >= 5
    val2 = varargin{2};
end
%% Download available files that do not already exist in local directory
% from ftp
mw = ftp(ftpSite);

for ddx = 1:length(dayList)
    dayNum = dayList(ddx);
    Year = YearList(ddx);
    
    jdi = doy2jd(Year,dayNum);
    [yri,mni,dyi] = jd2cal(jdi);
    [gpsWeek,gpsTow] = jd2gps(jdi);
    gpsDow = floor(gpsTow/86400);
    
    % Initial week of year
    [gpsWeek0,tow0] = jd2gps(cal2jd(Year,1,1));
%     if tow0 ~= 0
%         gpsWeek0 = gpsWeek0+1;
%     end
    woy = gpsWeek-gpsWeek0+1;
    
    target_dir = [destDir eval(destFormat)];
    
    ftpDir = eval(sourceFormat);
    
    % If we're not currently in the desired folder, go there.
    if ~strcmp(cd(mw),ftpDir)
        cd(mw,ftpDir);
    end
    
    if ~exist(target_dir,'dir') && strcmp(fileFormat{1},'[''*'']')
        mget(mw, '*', target_dir);
        change = 1;
    else
        %check to see what we already have
        llist = dir(target_dir);
        
        lname = {llist(:).name};
        ldate = [llist(:).datenum];
        
        change = 0;
        rlist = dir(mw);
        for i = 1:length(rlist)
            %only download if we do not have it or the remote version is newer
            
            serverName = rlist(i).name;
            if isempty(lname)
                have = 0;
            else
                [have, idx] = ismember(serverName,lname);
            end
            
            % Check if this is one of the files we want
            desired = 0;
            
            for fdx = 1:length(fileFormat)
                desiredName = regexptranslate('wildcard',eval(fileFormat{fdx}));
                
                desired = desired || ~isempty(regexp(serverName, desiredName,'once'));
                
                if desired
                    break
                end
            end
            
            if (~have || rlist(i).datenum > ldate(idx)) && desired
                mget(mw, serverName, target_dir);
                change = 1;
                
                if unzipFlag
                    unzipFile([target_dir '\' serverName]);
                end
            end
          
        end
        
    end
    if change
        disp(['File(s) updated on day ' int2str(dayNum)]);
        YearChange = [YearChange; Year];
        dayChange  = [dayChange; dayNum];
    end
    
end
close(mw);


end