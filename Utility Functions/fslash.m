function strOut = fslash(strIn)
% Change all slashes to forward slashes

strOut = strIn;

backInds = strfind(strIn,'\');

strOut(backInds) = '/';


end