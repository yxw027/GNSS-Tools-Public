function [ionoparams, clockParams] = readRinexIonoParams(filename)
% Pulls ionospheric parameters from RINEX navigation file if available. 
% If unavailable, it crashes. 
ionoparams = [];
clockParams = nan(6,4);

% File needs to exist!
if ~exist(filename,'file')
    return;
end

fid = fopen(filename);

if fid == -1
    return;
end

while 1
    tline = fgetl(fid);
    
    % Look for Iono params
    if contains(tline,'ION ALPHA')
        % pull out iono alpha params
        
        ionoparams(1:4) = cell2mat(textscan(tline,'%f %f %f %f'));
%         ionoparams(1) = str2double(tline(4:14));
%         ionoparams(2) = str2double(tline(16:26));
%         ionoparams(3) = str2double(tline(28:38));
%         ionoparams(4) = str2double(tline(40:50));        
    end
    if contains(tline,'ION BETA')
        % pull out beta alpha params
                ionoparams(5:8) = cell2mat(textscan(tline,'%f %f %f %f'));
%         ionoparams(5) = str2double(tline(4:14));
%         ionoparams(6) = str2double(tline(16:26));
%         ionoparams(7) = str2double(tline(28:38));
%         ionoparams(8) = str2double(tline(40:50));
    end
    if contains(tline,'GAUT')
        % pull out GAL 2 UTC
        clockParams(1,1) = str2double(tline(6:22));
        clockParams(1,2) = str2double(tline(23:38));
        clockParams(1,3) = str2double(tline(40:45));
        clockParams(1,4) = str2double(tline(47:50));
    end
    if contains(tline,'GPUT')
        % pull out GPS 2 UTC
        clockParams(2,1) = str2double(tline(6:22));
        clockParams(2,2) = str2double(tline(23:38));
        clockParams(2,3) = str2double(tline(40:45));
        clockParams(2,4) = str2double(tline(47:50));
    end
    
    if contains(tline,'SBUT')
        % pull out SBAS 2 UTC
        clockParams(3,1) = str2double(tline(6:22));
        clockParams(3,2) = str2double(tline(23:38));
        clockParams(3,3) = str2double(tline(40:45));
        clockParams(3,4) = str2double(tline(47:50));
    end
    
    if contains(tline,'GLUT')
        % pull out GLONASS 2 UTC
        clockParams(4,1) = str2double(tline(6:22));
        clockParams(4,2) = str2double(tline(23:38));
        clockParams(4,3) = str2double(tline(40:45));
        clockParams(4,4) = str2double(tline(47:50));
    end
    if contains(tline,'GPGA')
        % pull out GPS 2 GAL
        clockParams(5,1) = str2double(tline(6:22));
        clockParams(5,2) = str2double(tline(23:38));
        clockParams(5,3) = str2double(tline(40:45));
        clockParams(5,4) = str2double(tline(47:50));
    end
    if contains(tline,'GLGP')
        % pull out GLONASS 2 GPS
        clockParams(6,1) = str2double(tline(6:22));
        clockParams(6,2) = str2double(tline(23:38));
        clockParams(6,3) = str2double(tline(40:45));
        clockParams(6,4) = str2double(tline(47:50));
    end
    
    if contains(tline,'END OF HEADER')
        break
    end
    
    
end
fclose(fid);
end