function biasOut = findDcbElement(PRNs,consts,obs1,obs2,epochs,dcbData,statCode,noNan)


if nargin < 7
    statCode = [];
end

if nargin < 8
    noNan = 0;
end

if noNan
    biasOut = zeros(size(PRNs));
else
    biasOut = nan(size(PRNs));
end

for idx = 1:length(PRNs)
    prni = PRNs(idx);
    consti = consts(idx);
    
    obs1i = obs1{idx};
    obs2i = obs2{idx};
    
    epochi = epochs(idx);
    if isempty(statCode)
        inds = find(dcbData.PRNs == prni & dcbData.constInd == consti & strFindCell(dcbData.obs1,obs1i,1,1) & ...
            strFindCell(dcbData.obs2,obs2i,1,1) & epochi >= dcbData.startEpoch & epochi <= dcbData.endEpoch);
    else
        inds = find(dcbData.PRNs == prni & dcbData.constInd == consti & strFindCell(dcbData.obs1,obs1i,1,1) & ...
            strFindCell(dcbData.obs2,obs2i,1,1) & epochi >= dcbData.startEpoch & epochi <= dcbData.endEpoch & ...
            strFindCell(dcbData.sites,statCode,1,1));
    end
    
    
    if ~isempty(inds)
        inds = inds(1);
        
        biasOut(idx) = dcbData.bias(inds)*1e-9;
    end
    
end


end