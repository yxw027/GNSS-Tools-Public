function initIonoData(obj,year,doy,settings,varargin)


[~,~,~,IFileNameFull] = loadIonex(year,doy,settings,1);
dlFlag = 0;
for idx = 1:length(IFileNameFull)
    if ~exist(IFileNameFull{idx},'file')
        dlFlag = 1;
    end
end
if dlFlag
    % Download necessary orbit products
    ftpHelper(15,year,doy,settings,settings.clkCenter{idx});
end
[ionoData,~,~,filenameIono] = loadIonex(year,doy,settings,0);


obj.iono = ionoData;

end