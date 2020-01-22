function outData = nanunique(inData)

outData = unique(inData);
outData(isnan(outData)) = [];

end