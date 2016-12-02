%extractFrequentSeq will look through files for frequent sequences either
%by highest number of unique sequence per group, largest number of template
%count per group, or maximum template count per group.
%
% GrpTempMat is GrpNum GrpSize SumTemplateCtOfGrp MaxTemplateCtInGrp
function [NewVDJdata,GrpTempMat] = extractFrequentSeq(varargin)
Option = [];
TopPicks = 1; %Select the first number of top picsk to extract

if ~isempty(varargin)
    CheckedInput = ones(size(varargin))==1;
    
    %Determine how you want the frequency to be collected
    for q = 1:length(varargin)
        if ischar(varargin{q})
            if strcmpi(varargin{q},'Option')
                Option = varargin{q+1};
                CheckedInput(q:q+1) = 0;
                break
            end
        end
    end

    %Determine how you want the frequency to be collected
    for q = 1:length(varargin)
        if ischar(varargin{q})
            if strcmpi(varargin{q},'TopPicks')
                TopPicks = varargin{q+1};
                CheckedInput(q:q+1) = 0;
                break
            end
        end
    end
    varargin = varargin(CheckedInput);
end
if isempty(Option)
    disp('Option 1 : most sequence diversity')
    disp('Option 2 : most template count per group')
    disp('Option 3 : highest maxmium template count per group')
    Option = input('Choose the option');
end

if isempty(varargin)
    [FileNames, FilePath] = uigetfile('*.xlsx','Select Files to look for shared seq','multiselect','on');
    if ~iscell(FileNames)
        FileNames = {FileNames};
    end
    SkipLoad = 0;
    Fnum = length(FileNames);
else
    VDJdata = varargin{1};
    NewHeader = varargin{2};
    SkipLoad = 1;
    Fnum = 1;
end

for f = 1:Fnum
    
    if SkipLoad == 0
        FileName = FileNames{f};
        [VDJdata,NewHeader,~,~] = openSeqData([FilePath FileName]);
    end
    
    getHeaderVar;
    
    GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
    GrpNumUnq = unique(GrpNum);
    GrpTempMat = zeros(length(GrpNumUnq),4);
    
    for y = 1:length(GrpNumUnq)
        GrpIdx = find(GrpNum == GrpNumUnq(y));
        TempCt = cell2mat(VDJdata(GrpIdx,TemplateLoc));
        GrpTempMat(y,1) = GrpNumUnq(y);
        GrpTempMat(y,2) = length(GrpIdx);
        GrpTempMat(y,3) = sum(TempCt);
        GrpTempMat(y,4) = max(TempCt);
    end
    
    %Want top picks
    if size(GrpTempMat,1) < TopPicks
        TopPicks = size(GrpTempMat,1);
    end
    GrpTempMat = sortrows(GrpTempMat,-(Option+1));
    SelectGrpNum = GrpTempMat(1:TopPicks,1);
    NumOfElem = sum(GrpTempMat(1:TopPicks,2));
    
    %Extract sequence with the same group and CDR3
    NewVDJdata = cell(NumOfElem,length(NewHeader));
    r = 1;
    for q = 1:length(SelectGrpNum)
        IdxLoc = find(GrpNum == SelectGrpNum(q));
        NewVDJdata(r:r+length(IdxLoc)-1,:) = VDJdata(IdxLoc,:);
        r = r+length(IdxLoc);
    end
    
    %Before saving to xlsx, convert columns with matrix values into char
    NewVDJdata = reformatAlignment(NewVDJdata,1);
    for q = 1:size(NewVDJdata,1)
        for w = 1:3
            NewVDJdata{q,FamNumLoc(w)} = mat2str(NewVDJdata{q,FamNumLoc(w)});
        end
    end
    
    if SkipLoad == 0
        %Saving just the frequentSeq and the group
        SaveName = ['FreqSeqOPT' num2str(Option) '_' FileName];
        xlswrite(SaveName,[NewHeader;NewVDJdata]);
    end
end
        