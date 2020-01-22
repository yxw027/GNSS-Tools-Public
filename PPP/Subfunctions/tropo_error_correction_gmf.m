function [trop0,m,tropDataSave] = tropo_error_correction_gmf(el,h,lat,lon,mjd,epochi)

% Compute Dhz
pres = zeros(size(lat));
temp = zeros(size(lat));
undu = zeros(size(lat));
for idx= 1:length(lat)
    [pres(idx), temp(idx), undu(idx)] = gpt(mjd, lat(idx)*pi/180, lon(idx)*pi/180, h(idx));
end

% Equation 9.4 in IERS TN 36
fs = 1-0.00266.*cos(2*lat)-0.00000028.*h;

Dhz = 0.0022768*pres./fs;

Dwz = zeros(size(Dhz));


% FILLING IN HERE LOL
% load('C:\Users\kgunning\Box Sync\Novatel Data\IGS Data\STFU60\ztdNrc.mat')
% indsFill = find(round(epochi) == epochs);
% ztdNrci = ztdNrc(indsFill);
% ztdNrci(isnan(ztdNrci)) = ztdNrc(2);
% Dwz = (ztdNrci-Dhz(1))*ones(size(Dhz));

% Mapping functions
gmfh = zeros(size(lat));
gmfw = zeros(size(lat));
for idx = 1:length(lat)
    [gmfh(idx),gmfw(idx)] = gmf_f_hu (mjd,lat(idx)*pi/180,lon(idx)*pi/180,h(idx),pi/2-el(idx)*pi/180);
end


trop0 = Dhz.*gmfh+Dwz.*gmfw;
m = gmfw;

tropDataSave.trototSave = Dhz;
tropDataSave.gmfwSave   = gmfw;
tropDataSave.tzd        = Dhz(1)+Dwz(1);

end
