function ionoData = loadCodeKlobIono(Year,dayNum,settings)


destDir      = settings.dcbMgexDir;
dirPath   = [int2str(Year) '\' num2str(dayNum,'%03d') '/'];
filenamei   =  ['CGIM' num2str(dayNum,'%03d') '0.' num2str(mod(Year,100),'%02d') 'N'] ;

ionoData = readRinexIonoParams([destDir dirPath filenamei]);


end