function manageSubsets(ss,obs,PARAMS)

%% Check what satellites are in view- this is the basis for the subset list
indsRange = find(obs.range.obs ~= 0);
indsDoppler = find(obs.doppler.obs ~= 0);
% Full list of satellites that require a subset
prnConstInds  = unique([obs.range.PRN(indsRange) obs.range.constInds(indsRange);
    obs.doppler.PRN(indsDoppler) obs.doppler.constInds(indsDoppler)],'rows');

nSats = size(prnConstInds,1);
for idx = 1:nSats
    if ~ismember(prnConstInds(idx,:),ss.ssGroupRows,'rows')
        % need to add this one satellite
        % check if there's space to add it in a group
        
         indAdd = min(find(isnan(squeeze(ss.ssGroups(:,:,1)))));
        if isempty(indAdd)
            groupAdd = nan(1,PARAMS.solSep.subsetGroupSize,size(ss.ssGroups,3));
            groupAdd(1,1,:) = prnConstInds(idx,:);
            
            ss.ssGroups = cat(1,ss.ssGroups,groupAdd);
            
        else
           % This can be added to an existing group 
           
           % In previous versions, this is added to a big exclusion matrix
           % here as well. 
           [rowAdd,colAdd] = ind2sub([size(ss.ssGroups,1) size(ss.ssGroups,2)],indAdd);
           
           ss.ssGroups(rowAdd,colAdd,:) = prnConstInds(idx,:);
           
           % should add this one to existing subsets to avoid adding an
           % entirely new subset ( unless the group is also new)
%            error('fix me');
        end
    end    
end

nPerGroup = sum(~isnan(squeeze(ss.ssGroups(:,:,1))),2);
pSatGroups = 1-(1-PARAMS.solSep.Psat).^(nPerGroup);

%% this doesnt really matter right now, but here is where you can add other types of measurements
posMeasAvail = 0;
if posMeasAvail
    posMeasInd = nan(1,PARAMS.solSep.subsetGroupSize,size(ss.ssGroups,3));
    posMeasInd(1) = 0;
    posMeasPsat = 0;
else
    posMeasInd = [];
    posMeasPsat = [];
end

ssGroupsFull = cat(1,posMeasInd, ss.ssGroups);
pSatFull     = cat(1,posMeasPsat,pSatGroups);

nGroupsOut = PARAMS.solSep.nOutSubset+PARAMS.solSep.faultResetSubsets;

nMeasOutMax = nGroupsOut.*PARAMS.solSep.subsetGroupSize;
ssNeeded     = [];
pSatNeeded   = [];
ssTypeNeeded = [];
%% Build the list of subsets to create- each has a certain number out as well as 
for idx = 1:nGroupsOut
    indsi  = nchoosek(1:size(ssGroupsFull,1),idx);

    ssGroupsi = [];
    pSatSsi = ones(size(indsi,1),1);
    for jdx = 1:size(indsi,2)
        ssGroupsi = cat(2,ssGroupsi, ssGroupsFull(indsi(:,jdx),:,:));
        pSatSsi = pSatSsi.*pSatFull(indsi(:,jdx));
    end
    
    % Put these in a nan-padded full size matrix to make concatenation
    % easier
    ssGroupsPad = nan(size(ssGroupsi,1),nMeasOutMax,size(ssGroupsi,3));
    ssGroupsPad(:,1:size(ssGroupsi,2),:) = ssGroupsi;

    if idx == nGroupsOut && PARAMS.solSep.faultResetSubsets
        typei = 2; % on-deck subset
    else
        typei = 1;
    end
    
    ssNeeded = cat(1,ssNeeded,ssGroupsPad);
    pSatNeeded = cat(1,pSatNeeded,pSatSsi);
    ssTypeNeeded = cat(1,ssTypeNeeded,typei*ones(size(ssGroupsi,1),1));
end

%% Add something here to remove old subsets
ssAlreadyExist = permute([ss.ssList.measOut],[2 1 3]);
ssAlreadyExist = reshape(ssAlreadyExist,size(ssAlreadyExist,1),size(ssAlreadyExist,2)*size(ssAlreadyExist,3),1);
ssNeededRows = reshape(ssNeeded,size(ssNeeded,1),size(ssNeeded,2)*size(ssNeeded,3),1);
if ~isempty(ssAlreadyExist)
    indsAdd = find(~ismember(ssNeededRows,ssAlreadyExist,'rows'));
else
    indsAdd = (1:size(ssNeededRows,1))';
end

%% Add necessary subsets
% only add enough so that we don't exceed the maximum subset count
nSubsetCurr = length(ss.ssList);
nSubsetMax = PARAMS.solSep.nMaxSubset;

indsAdd = indsAdd(1:min([nSubsetMax-nSubsetCurr length(indsAdd)]));

for idx = 1:length(indsAdd)
    indNewSs =length(ss.ssList)+1; 
    ss.ssList(indNewSs).subsetType = ssTypeNeeded(indsAdd(idx));
    ss.ssList(indNewSs).ekf  = copy(ss.ssList(1).ekf);
    ss.ssList(indNewSs).aivFlag = 0;
    ss.ssList(indNewSs).measOut = permute(ssNeeded(indsAdd(idx),:,:),[2 1 3]);
    ss.ssList(indNewSs).pSat    = pSatNeeded(indsAdd(idx));
    ss.ssList(indNewSs).ID      = ss.lastID+1;
    
    ss.lastID = ss.lastID+1;
end

end