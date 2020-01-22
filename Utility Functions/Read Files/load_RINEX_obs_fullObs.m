function [obsStruc, time_ref, time, week, date, pos, interval, antoff, antmod,...
    rxmod,obsColumns,rinexVer] = ...
    load_RINEX_obs_fullObs(filename, constellations, wait_dlg,rnx3Format,varargin)

% SYNTAX:
%   [pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2, ...
%    time_ref, time, week, date, pos, interval, antoff, antmod] = ...
%    load_RINEX_obs(filename, constellations, wait_dlg);
%
% INPUT:
%   filename = RINEX observation file(s)
%   constellations = struct with multi-constellation settings
%                   (see 'multi_constellation_settings.m' - empty if not available)
%   wait_dlg = optional handler to waitbar figure (optional)
%
% OUTPUT:
%   pr1 = code observation (L1 carrier)
%   ph1 = phase observation (L1 carrier)
%   pr2 = code observation (L2 carrier)
%   ph2 = phase observation (L2 carrier)
%   dop1 = Doppler observation (L1 carrier)
%   dop2 = Doppler observation (L2 carrier)
%   snr1 = signal-to-noise ratio (L1 carrier)
%   snr2 = signal-to-noise ratio (L2 carrier)
%   time = receiver seconds-of-week
%   week = GPS week
%   date = date (year,month,day,hour,minute,second)
%   pos = rover approximate position
%   interval = observation time interval [s]
%   antoff = antenna offset [m]
%   antmod = antenna model [string]
%
% DESCRIPTION:
%   Parses RINEX observation files.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2013 Mirko Reguzzoni,Eugenio Realini
% Portions of code contributed by Damiano Triglione and Stefano Caldera
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

global report

% Check the input arguments
if (nargin < 3)
    wait_dlg_PresenceFlag = false;
else
    wait_dlg_PresenceFlag = wait_dlg;
end

if nargin < 4
    rnx3Format  = false;
end

if (nargin < 2 || isempty(constellations)) %then use only GPS as default
    [constellations] = initConstellation(1, 0, 0, 0, 0, 0);
end

p = inputParser;
validNumberFn = @(x) (x > 0) && isnumeric(x);
p.addParameter('headerOnly', false);  

% parse the results
parse(p, varargin{:});
res = p.Results;
headerOnly = res.headerOnly;

%number of satellite slots for enabled constellations
nSatTot = constellations.nEnabledSat;

%number of RINEX files to be read
if (iscell(filename))
    nFiles = size(filename,1);
else
    nFiles = 1;
end

%variable initialization
nEpochsAdd = 3000;
nEpochs = 3000;
time = NaN(nEpochs,1,nFiles);
tow = NaN(nEpochs,1,nFiles);
week = NaN(nEpochs,1,nFiles);
date = NaN(nEpochs,6,nFiles);
pos = zeros(3,1,nFiles);
interval = zeros(1,1,nFiles);
antoff = zeros(3,1,nFiles);
antmod = cell(1,1,nFiles);
rxmod = cell(1,1,nFiles);

for f = 1 : nFiles
    
    if (iscell(filename))
        current_file = filename{f,1};
    else
        current_file = filename;
    end
        
    %open RINEX observation file
    fid = fopen(current_file,'r');
    
%     try
    
    if (wait_dlg_PresenceFlag)
        waitbar(0.5,wait_dlg,['RINEX file ' current_file ': parsing header...'])
    end
    
    %parse RINEX header
    [obs_type, pos(:,1,f), basic_info, interval(1,1,f), sysId, antoff(:,1,f), ...
        antmod{1,1,f}, ~,rxmod{1,1,f},rinexVer] = RINEX_parse_hdr(fid);
    
    %check the availability of basic data to parse the RINEX file
    if (basic_info == 0)
        error(['RINEX file ' current_file ': basic data is missing in the file header'])
    end
    
    % Pull each signal name out
    [obsColumns,nObsTypes,obsTypes ] = obs_type_find_any(obs_type,sysId);
    
    % Convert to RINEX 3 Codes if necessary and desired
    if isempty(sysId) && rnx3Format
        constTypes = {'G','R','E','C','J','S'};
        
        obsTypes3 = {};
        for cdx = 1:length(constTypes)
            obsColumnsi = convertRinex3ObsCodes(obsColumns,constTypes{cdx});
            obsColumns3.(constTypes{cdx}) = obsColumnsi;
            
            obsTypes3 = unique([obsTypes3; obsColumnsi(~cellfun(@isempty,obsColumnsi))]);
        end
        
        obsColumns = obsColumns3;
        obsTypes = obsTypes3;
        sysId = 'CONVERTED_3';        
    end
    
    if f == 1
       % initialize storage variables
       for odx = 1:length(obsTypes)
          obsStruc.(obsTypes{odx}) = NaN(nSatTot,nEpochs,nFiles);
          
       end
       obsOut = nan(nSatTot,nEpochs,nFiles,length(obsTypes));
    end

    if isempty(sysId) % RINEX v2.xx
        obsColumnsMat = repmat(1:length(obsTypes),6,1);
        
    else % RINEX v3.xx
        constTypes = {'G','R','E','C','J','S'};
        obsColumnsMat = zeros(length(constTypes),length(obsTypes));
        for cdx = 1:length(constTypes)
            for odx = 1:length(obsTypes)
                if isfield(obsColumns,constTypes{cdx})
                    %                 coli = obsColumns.(constTypes{cdx}).(obsTypes{odx});
                    coli = find(~cellfun(@isempty,strfind(obsColumns.(constTypes{cdx}) , obsTypes{odx})));
                    
                    if ~isempty(coli)
                        obsColumnsMat(cdx,odx) = coli(1);
                    end
                end
            end
        end        
    end
    
    if (wait_dlg_PresenceFlag)
        waitbar(1,wait_dlg)
    end
    
    %-------------------------------------------------------------------------------
    
    if (wait_dlg_PresenceFlag)
        waitbar(0.5,wait_dlg,['RINEX file ' current_file ': reading observations...'])
    end
    
    if headerOnly
        % if we only want the header, ignore the rest of the file
        fclose(fid);
        continue;
    end
    
    k = 1;
    while (~feof(fid))
        %read data for the current epoch (ROVER)
        try
            [time(k,1,f), date(k,:,f), num_sat, sat, sat_types, tow(k,1,f)] = RINEX_get_epoch(fid);
        catch
            % Occasionally, files are busted. 
             'fdafas';
            continue
           
        end    
           
        if (k > nEpochs)
            obsOut2 = nan(size(obsOut,1),nEpochs+nEpochsAdd,size(obsOut,3),size(obsOut,4));
            obsOut2(:,1:size(obsOut,2),:,:) = obsOut;
            obsOut = obsOut2;
            obsOut2 = [];
            
            date2 = nan(nEpochs+nEpochsAdd,size(date,2),size(date,3));
            date2(1:size(date,1),:,:) = date;
            date = date2;
            date2 = [];
            
            tow2 = nan(nEpochs+nEpochsAdd,size(tow,2),size(tow,3));
            tow2(1:size(tow,1),:) = tow;
            tow = tow2;
            tow2 = [];

            time2 = nan(nEpochs+nEpochsAdd,size(time,2),size(time,3));
            time2(1:size(time,1),:) = time;
            time = time2;
            time2 = [];
            
            week2 = zeros(nEpochs+nEpochsAdd,size(week,2),size(week,3));
            week2(1:size(week,1),:) = week;
            week = week2;
            week2 = [];
            
            nEpochs = nEpochs  + nEpochsAdd;
        end
        
        %read ROVER observations
        [~,obsMati] = RINEX_get_obs_fullObs(fid, num_sat, sat, sat_types, obsColumns, ...
            nObsTypes, constellations,obsColumnsMat,obsTypes);
        
        obsOut(:,k,f,:) = obsMati;
        
        k = k + 1;
    end
    
    if (wait_dlg_PresenceFlag)
        waitbar(1,wait_dlg)
    end
    
    %GPS week number
    week(:,1,f) = date2gps(date(:,:,f));
    
    %observation rate
    if (interval(:,1,f) == 0)
        interval(:,1,f) = round((median(time(2:k-1,1,f) - time(1:k-2,1,f)))*1000)/1000;
    end
    
%     -------------------------------------------------------------------------------
%     catch ME
%          disp(['Error- filename = ' filename])
%          disp(ME.message)
%          
%     end
    %close RINEX files
    fclose(fid);
    
end

for odx = 1:length(obsTypes)
   obsStruc.(obsTypes{odx}) = squeeze(obsOut(:,:,:,odx)); 
end

[time_ref, time, week, date, obsStruc, interval] = ...
    sync_obs_fullObs(time, week, date, obsStruc, interval);

end

