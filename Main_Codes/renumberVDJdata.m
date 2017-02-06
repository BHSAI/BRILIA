%renumberVDJdata will renumber AND reorder the sequence numbers in VDJdata
%to go chronologically in order per each grp num.
%
%  VDJdata = renumberVDJdata(VDJdata,VDJheader,Option1,Option2)
%
%  INPUT
%    Option1 ['seq']: Renumber sequence numbers
%    Option2 ['grp']: Renumber group numbers
%
%  NOTE
%    IF no options are specified, then will only sort groups together.
%
%    You can sepecify both 'seq' and 'grp' in any order if you want both to
%    be sorted.
%
%    It will ALWAYS resort by group number, but NEVER by sequence number,
%    because the 1st Seq in group is assumed most ancestral sequence.
%  
%  EXAMPLES
%    VDJdata = num2cell([[1 3 4 7 2 8 13]' [1 1 1 3 1 3 3]']);
%    VDJheader = {'SeqNum','GroupNum'};
%
%    VDJdata = renumberVDJdata(VDJdata,VDJheader,'seq')
%    SeqNum    GrpNum  -->   SeqNum    GrpNum
%    1         1             1         1
%    3         1             2         1
%    4         1             3         1
%    7         3             4         1
%    2         1             5         3
%    8         3             6         3
%    13        3             7         3
%
%    VDJdata = renumberVDJdata(VDJdata,VDJheader,'grp')
%    SeqNum    GrpNum  -->   SeqNum    GrpNum
%    1         1             1         1
%    3         1             3         1
%    4         1             4         1
%    7         3             2         1
%    2         1             7         2
%    8         3             8         2
%    13        3             13        2
%
%    VDJdata = renumberVDJdata(VDJdata,VDJheader,'seq','grp')
%    SeqNum    GrpNum  -->   SeqNum    GrpNum
%    1         1             1         1
%    3         1             2         1
%    4         1             3         1
%    7         3             4         1
%    2         1             5         2
%    8         3             6         2
%    13        3             7         2

function VDJdata = renumberVDJdata(VDJdata,VDJheader,varargin)
H = getHeaderVar(VDJheader);

%Determine which options were given
SeqOption = 0; %Default, renumber sequences. Logical value.
GrpOption = 0; %Default, renumber group numbers. Logical value.
if max(findCell(varargin,'seq','MatchCase','any')) > 0
    SeqOption = 1;
end
if max(findCell(varargin,'grp','MatchCase','any')) > 0
    GrpOption = 1;
end

%Begin grp sorting, preserving relative sequence number
GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
GrpMap = [GrpNum [1:length(GrpNum)]'];
[GrpMap,GrpSortIdx] = sortrows(GrpMap,[1 2]); %Sort by group, then temp seq order number
GrpNum = GrpMap(:,1);
VDJdata = VDJdata(GrpSortIdx,:);

%Renumber the group numbers
if GrpOption == 1
    GrpCt = 1;
    CurGrp = GrpNum(1);
    for j = 1:length(GrpNum)
        if GrpNum(j) ~= CurGrp
            GrpCt = GrpCt + 1;
            CurGrp = GrpNum(j);
        end
        GrpNum(j) = GrpCt;
    end   
    VDJdata(:,H.GrpNumLoc) = num2cell(GrpNum);
end

%Renubmer to sequence numbers
if SeqOption == 1 
    VDJdata(:,H.SeqNumLoc) = num2cell(1:size(VDJdata,1));
end
