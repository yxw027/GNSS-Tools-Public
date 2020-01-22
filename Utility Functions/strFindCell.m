function indsFind = strFindCell(str,pattern,noFind,exact)

if nargin < 3
    noFind = 0;
end

if nargin < 4
    exact = 0;
end

if ~exact
    if noFind
        indsFind = ~cellfun(@isempty, strfind(str,pattern));
    else
        indsFind = find(~cellfun(@isempty, strfind(str,pattern)));
    end
else
    if noFind
        indsFind = strcmp(str,pattern);
    else
        indsFind = find(strcmp(str,pattern));
    end
end


end