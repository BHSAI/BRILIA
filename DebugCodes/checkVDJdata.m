%checkVDJdata is more of a debugging script made to check entries of
%VDJdata to ensure everything makes sense. If it encounters a bad set, will
%save the bad entries into a [DebugFileName].mat file.
%
%  [BadLoc, BadReason] = checkVDJdata(VDJdata,NewHeader,DebugFileName,DebugModeOn)
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
%      only VDJdata entries where there is an error.
%
%    You could use this code to filter out bad entries, by using the
%      following commands:
%      >> [BadLoc, BadReason] = checkVDJdata(VDJdata,NewHeader,'debug1');
%      >> VDJdata(BadLoc,:) = [];
function [BadLoc, BadReason] = checkVDJdata(VDJdata,NewHeader,DebugFileName,DebugModeOn)
BadLoc = zeros(size(VDJdata,1),1,'logical');
BadReason = cell(size(VDJdata,1),1);

if DebugModeOn == 0
    return
end

getHeaderVar;
for j = 1:size(VDJdata,1)
    %Check if Seq and RefSeq have same lengths, or RefSeq is empty
    Seq = VDJdata{j,SeqLoc};
    RefSeq =VDJdata{j,RefSeqLoc};
    
    if ~isempty(RefSeq)
        if length(RefSeq) ~= length(Seq)
            BadLoc(j) = 1;
            BadReason{j} = 'Length of RefSeq and Seq not the same';
            continue;
        end
    end
    
    %Check if segment lengths make sense
    VMDNJ = cell2mat(VDJdata(j,LengthLoc));
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
    if VMDNJ(5) <= 0
        BadLoc(j) = 1;
        BadReason{j} = 'Jlen is 0';
        continue;
    end
    if VMDNJ(3) <= 0
        BadLoc(j) = 1;
        BadReason{j} = 'Dlen is 0';
        continue;
    end
        
    %Check if deletion counts make sense
    VDDJdel = cell2mat(VDJdata(j,DelLoc));
    if min(VDDJdel) < 0
        BadLoc(j) = 1;
        BadReason{j} = 'Negative value in V3 D5 D3 J5 deletion count';
        continue;
    end
end

%Will save the bad VDJdata only, but will not modify VDJdata
if max(BadLoc > 0)
    VDJdata = VDJdata(BadLoc,:);
    VDJdata(:,MiscLoc) = BadReason(BadLoc);
    save([DebugFileName '.mat'],'VDJdata','NewHeader');
end