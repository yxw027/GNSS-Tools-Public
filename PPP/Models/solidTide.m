function [posOffset,rangeOffset] = solidTide(epoch,pos,varargin)

% Parse optional inputs
p = inputParser;
p.addParameter('sunPos',  []); % where the sun at (ECEF, meters)
p.addParameter('moonPos', []); % where the moon at (ECEF, meters)
p.addParameter('svPos',   []); % where the satellites at  (ECEF, meters)

parse(p, varargin{:});
res = p.Results;
sunPos  = res.sunPos;
moonPos = res.moonPos;
svPos   = res.svPos;


%% Initialize outputs
posOffset   = nan(3,1); % User position offset
rangeOffset = nan(3,1); % Offset per range

%% If sun and moon position are not available, need to find them
if isempty(sunPos) || isempty(moonPos)
    jd = epochs2jd(epoch);
    et = cspice_str2et(['jd ' num2str(jd)]);
end

if isempty(sunPos)
    sunPos   = cspice_spkezr( 'sun',et , 'itrf93', 'none', 'earth')*1000;
end

if isempty(moonPos)
    moonPos  = cspice_spkezr( 'moon',et , 'itrf93', 'none', 'earth')*1000;
end

%% need time in year, month, day, hour of day
[yr,mn,dy,hour,min,sec] = epochs2cal(epoch);
hour = hour+min/60+sec/3600;

% Solid tides
posOffset = dehanttideinel(pos,yr,mn,dy,hour,sunPos,moonPos)';

%% If satellite positions are available and the output is requested, also 
% include the offsets per range
if nargout > 1 && ~isempty(svPos)
    los = pos-svPos';
    los = los./(sqrt(sum(los.^2,1)));
    
    rangeOffset = sum(los.*posOffset,1)';    
elseif nargout > 1 && isempty(svPos)
    error('Satellite positions required to output solid tide range offset')
end

end




















