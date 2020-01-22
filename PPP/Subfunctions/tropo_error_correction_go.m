function [tropo_corr, Mw,tropDataExtra] = tropo_error_correction_go(el, h, lat,doy)

a = 2.3;
b = 0.116e-3;
d_dry = a.*exp(-b.*h);

[Mw, Md] = mappingOfNiell(el, h, lat,doy);

Mw = Mw;
Md = Md;

% Mw = m(:,1);
% Md = m(:,2);
d_wet = 0.1;
tropo_corr = Md.*d_dry+Mw.*d_wet;
tropo_corr = tropo_corr;

Mw = Mw;

tropDataExtra.trototSave = d_dry;
tropDataExtra.gmfwSave   = Mw;
tropDataExtra.tzd        = d_dry(1)+d_wet(1);


end























