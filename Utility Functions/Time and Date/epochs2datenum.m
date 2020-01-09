function tDatenum = epochs2datenum(epochs)
% Convert from GPS epoch to MATLAB datenum

tDatenum = datenum(epochs2cal(epochs,1));


end