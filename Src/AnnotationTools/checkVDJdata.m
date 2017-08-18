%checkVDJdata is more of a debugging script made to check entries of
%VDJdata to ensure everything makes sense. If it encounters a bad set, will
%save the bad entries into a [DebugFileName].mat file.
%
%  [BadLoc, BadReason] = checkVDJdata(VDJdata,VDJheader,DebugFileName,DebugModeOn)
%
%  INPUT
%    DebugFileName: prefix name for the debug file to generate when an
%      erroneous entry is encountered
%    DebugModeOn: either 1 or 0 to turn on or off debugging.
%
%  OUTPUT
%    BadLoc: a Mx1 logical matrix showing where an error was encounter
%    BadReason: a Mx1 cell matrix of text of the type of error encountered.
%    
%  NOTE
%    Outputs automatically generated in the current working directory of
%    only VDJdata entries where there is an error.
%
%    You could use this code to filter out bad entries, by using the
%    following commands:
%      >> [BadLoc, BadReason] = checkVDJdata(VDJdata,VDJheader,'debug1');
%      >> VDJdata(BadLoc,:) = [];
function [BadLoc, BadReason] = checkVDJdata(VDJdata,VDJheader,DebugFileName,DebugModeOn)
%Initialize outputs
BadLoc = zeros(size(VDJdata,1),1,'logical');
BadReason = cell(size(VDJdata,1),1);

%Ensure debug mode is on
if DebugModeOn == 0
    return
end

%Determine chain and extract key locations
H = getHeavyHeaderVar(VDJheader);
L = getLightHeaderVar(VDJheader);
if H.SeqLoc > 0 && L.SeqLoc > 0
    Chain = 'HL';
elseif H.SeqLoc > 0
    Chain = 'H';
elseif L.SeqLoc > 0
    Chain = 'L';
end

for k = 1:length(Chain)
    if Chain(k) == 'H'
        B = H;
    else
        B = L;
    end
    
    %Check each entry for validity
    for j = 1:size(VDJdata,1)
        %Check if Seq and RefSeq have same lengths
        Seq = VDJdata{j,B.SeqLoc};
        RefSeq = VDJdata{j,B.RefSeqLoc};

        if ~isempty(RefSeq) && length(RefSeq) ~= length(Seq)
            BadLoc(j) = 1;
            BadReason{j} = 'Length of RefSeq and Seq not the same';
            continue;
        end

        %Check if segment lengths make sense
        VMDNJ = cell2mat(VDJdata(j,B.LengthLoc));
        if sum(VMDNJ) ~= length(Seq)
            BadLoc(j) = 1;
            BadReason{j} = 'Length of Seq differs from the sum of VMDNJ';
            continue;
        end
        if min(VMDNJ) < 0 
            BadLoc(j) = 1;
            BadReason{j} = 'Negative value in VMDNJ';
            continue;
        end
        if VMDNJ(1) <= 0
            BadLoc(j) = 1;
            BadReason{j} = 'Vlen is 0';
            continue;
        end
        if VMDNJ(end) <= 0
            BadLoc(j) = 1;
            BadReason{j} = 'Jlen is 0';
            continue;
        end
        if length(VMDNJ) > 3 && VMDNJ(3) <= 0 %For heavy chain only
            BadLoc(j) = 1;
            BadReason{j} = 'Dlen is 0';
            continue;
        end

        %Check if deletion counts make sense
        VDDJdel = cell2mat(VDJdata(j,B.DelLoc));
        if min(VDDJdel) < 0
            BadLoc(j) = 1;
            BadReason{j} = 'Negative value in gene deletion count';
            continue;
        end
    end
end

%Will save the bad VDJdata only, but will not modify VDJdata
if max(BadLoc > 0) && ~isempty(DebugFileName)
    %Setup the save file name
    [FilePath,FileName,FileExt] = parseFileName(DebugFileName);
    if isempty(FileExt); 
        FileExt = '.csv';
        FileName = sprintf('%s%s',FileName,FileExt);
    end
    FullSaveName = [FilePath FileName];
    
    %Prep VDJdata and VDJheader for saving with error message
    VDJdata = cat(2,VDJdata(BadLoc,:),BadReason(BadLoc));
    VDJheader = cat(2,VDJheader,{'ErrMsg'});
    saveSeqData(FullSaveName,VDJdata,VDJheader,'delimiter',';');
end
