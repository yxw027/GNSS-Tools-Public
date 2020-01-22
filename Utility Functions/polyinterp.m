function [yy, yy_dot, chi2, p] = polyinterp(x, y, m_order, xx, flag, var,polyIn)


if nargin < 6
    var = ones(size(y));
elseif var == 0 | isempty(var)
    var = ones(size(y));
end

if nargin > 4
    %remove flagged data and NaN
    idx  = ~flag & ~isnan(x) & ~isnan(y);
else
    idx  = ~isnan(x) & ~isnan(y);    
end
if nargin < 7
   polyIn = []; 
end

x = x(idx);
y = y(idx);
var = var(idx);

% make sure enough points remain to fit
if length(y) < m_order + 1
%     warning('Warning bad fit in polyinterp\n');
    yy     = NaN(size(xx));
    yy_dot = NaN(size(xx));
    chi2   = NaN;    
    p = NaN(m_order+1,1);
    return
end



xscale = (max(x) - min(x))/2;
xbias  = mean(x);
xp = (x-xbias)/xscale;

yscale = (max(y) - min(y))/2;
ybias  = mean(y);
if yscale ~= 0
    yp = (y-ybias)/yscale;
else
    yp = y*0;
end

n = size(xx,1);
xx = [xx; x];
if isempty(polyIn)
    p = polyfit2(xp, yp, m_order);
else
    p = polyIn;
end
z = yscale*polyval(p, (xx - xbias)/xscale) + ybias;
yy = z(1:n);

if nargout > 1
    pdot = p(1:end-1).*(m_order:-1:1);
    yy_dot = yscale*polyval(pdot, (xx(1:n) - xbias)/xscale)/xscale;
    if nargout > 2
        chi2 = (y - z((n+1):end))'*((y - z((n+1):end))./var);
    end
end




















