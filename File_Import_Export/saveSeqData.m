%saveSeqData will save VDJdata
%
%  saveSeqData(FullSaveName,VDJdata,VDJheader)
%  saveSeqData(FullSaveName,VDJdata,VDJheader,'Delimiter','\t')
%
%  Delimiter can be ',' or ';' or '\t'

function saveSeqData(FullSaveName,VDJdata,VDJheader,varargin)
P = inputParser;
addParameter(P,'Delimiter','\t',@(x) ismember(x,{';' ',' '\t'}));
parse(P,varargin{:});
Delimiter = P.Results.Delimiter;

H = getHeaderVar(VDJheader);

[~,~,FileExt] = parseFileName(FullSaveName);

%Before saving to xlsx, convert columns with matrix values into char
for q = 1:size(VDJdata,1)
    for w = 1:3
        VDJdata{q,H.FamNumLoc(w)} = mat2str(VDJdata{q,H.FamNumLoc(w)});
    end
end

%Save to excel or csv file, depending on OS
if ~isempty(regexpi(FileExt,'.xls','once'))
    xlswrite(FullSaveName,cat(1,VDJheader,VDJdataSave));
elseif ~isempty(regexpi(FileExt,'.csv','once'))
    writeDlmFile([VDJheader;VDJdata],FullSaveName,Delimiter);
end
