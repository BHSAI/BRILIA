%sortNTseq takes a trimmed, aligned group of sequences from an Excel file
%and rearranges it based on how similar the sequences are from each other.
%
%  Output = sortNTseq()   This will ask user to select a file to open,
%  assuming Column 1 = Seq. The outputfile is an Excel file with these
%  data, plus a last column with the group number. Similarity is determined
%  as 6 NTs different between members within the same family.
%
%  Output = sortNTseq(6,NTseq) will process the input of sequence, grouping
%  Hamming distances of 6 or less into same group.
%
%  Output = sortNTseq(6,NTseq,'Option',[ham or shmham] will calculate
%  pairwise distance using shmham
%  
%  NOTE! Consecutive mismatches causes a ^2 hamming distance increases.
%    EX:    001100 is a hamming distance of 4, NOT 2.

function NTmaps = sortNTseq(varargin)

%Extract the inputs
    %Case 1) Maximum deviation is given, but no file name
    %Case 2) Maximum deviation + file name given
    %Case 3) Maximum deviation + NTseq given

%Looking for Maximum deviation info, if none, ask the user. 
Option = 'shmham';
MaxDev = 1;
if nargin == 0
    MaxDev = input('What is the Hamming distance to cluster?');
    if isempty(MaxDev)
        MaxDev = 1;
    else
        MaxDev = round(abs(MaxDev));
        fprintf('Hamming Distance is: %d \n', MaxDev')
    end
else
    DelThis = zeros(1,length(varargin)) > 1;
    for v = 1:length(varargin)
        if isnumeric(varargin{v})
            MaxDev = round(abs(varargin{v}));
            DelThis(v) = 1; %Remove from varargin
        elseif ischar(varargin{v})
            if strcmpi(varargin{v},'option')
                Option = varargin{v+1};
                DelThis(v:v+1) = 1; %Remove from varargin
            end
        end
    end
    varargin(DelThis) = [];
end

%Looking for sequence information as cells
if isempty(varargin) %Ask the user for the file Name
    [VDJdata, NewHeader, FileName, FilePath] = openSeqData();
    getHeaderVar;
    NTseq = VDJdata(:,SeqLoc);
else
    if iscell(varargin{1}) %Probably sequences
        NTseq = varargin{1};
    else
        [VDJdata, NewHeader, FileName, FilePath] = openSeqData(varargin{v});
        getHeaderVar;
        NTseq = VDJdata(:,SeqLoc);
    end
end

%Preallocate scoring matrix, using uint8,16,32 depending on number of seq
Bits = [8 16 32];
BitLoc = sum(size(NTseq,2) >= 2.^Bits)+1;

%Calculate matrix of NT differences
SeqNum = size(NTseq,1);
switch lower(Option)
    case 'ham'
        ScoreMat = zeros(size(NTseq,1),['uint' num2str(Bits(BitLoc))]);
        if iscell(NTseq)
            NTseq = char(NTseq);
        end
        SeqLen = size(NTseq,2);
        for x = 1:size(NTseq,1)   
            CompareMat = repmat(NTseq(x,:),SeqNum-x+1,1) == NTseq(x:end,:);
            ScoreMat(x:end,x) = SeqLen - sum(CompareMat,2);    
            ScoreMat(x,x:end) = ScoreMat(x:end,x);
        end
    case 'shmham'
        ScoreMat = calcPairDist(NTseq,'shmham'); %,zeros(4)%!!!THIS CAN BE ALTERED, to true SHMHAM score by removing "zeros(4)"
        ScoreMat(eye(size(ScoreMat))>0) = 0;
        MaxDev = MaxDev*2; %Remember, we had to double shmham scores, to fit to uint formats!
end

%Perform grouping
NTgrp = zeros(size(NTseq,1),1); %Index of group #
FmapIdx = 1:SeqNum; %Index of entries to accelerate find command.
FmapLog = ones(1,SeqNum)==1; %Binary location of unassigned entries.
Group = 1;
while 1
    CurLoc = FmapIdx(FmapLog);
    if isempty(CurLoc); break; end
    CurLoc = CurLoc(1);
    
    %Look for initial group member besides self
    AllMatch = ((ScoreMat(CurLoc,:) <= MaxDev) & FmapLog);
    
    %Mark the initial group size
    AllPath = sum(AllMatch);

    %Look at each hit, and look for another hit close to those hits
    while 1
        FindMatch = FmapIdx(AllMatch);
        for k = 1:length(FindMatch)
            AllMatch = (AllMatch | (ScoreMat(FindMatch(k),:) <= MaxDev)) & FmapLog;
        end
        if sum(AllMatch) == AllPath %No more paths available
            break
        else
            AllPath = sum(AllMatch);%Redo as new paths were found
        end
    end
    
    %Save the group numbers and update FmapLog.
    NTgrp(AllMatch) = Group;
    Group = Group + 1;
    FmapLog(AllMatch) = 0;
end

%Resort everything according to ascending GroupNum
NTmaps = sortrows([NTgrp [1:size(NTseq,1)]']);

if exist('FileName','var')
    VDJdata = VDJdata(NTmaps(:,2),:);

    %Add the SeqNum and GroupNum into the file.
    SeqNumLoc = findHeader(NewHeader,'SeqNum');
    if SeqNumLoc == 0    
        NewHeader = cat(2,NewHeader,{'SeqNum'});
        VDJdata = cat(2,VDJdata,num2cell(NTmaps(:,2)));
    else
        %IF there was a seqnum before, you need to map another map, so you can
        %always track back to the original seq number.
        if ~isempty(VDJdata{1,SeqNumLoc})
            PrevMap = cell2mat(VDJdata(:,SeqNumLoc));
            VDJdata(:,SeqNumLoc) = num2cell(PrevMap(NTmaps(:,2)));
        else
            VDJdata(:,SeqNumLoc) = num2cell(NTmaps(:,2));
        end
    end

    %Add the Group Num information
    GroupNumLoc = findHeader(NewHeader,{'GroupNum'});
    if GroupNumLoc == 0
        NewHeader = cat(2,NewHeader,{'GroupNum'});
        VDJdata = cat(2,VDJdata,num2cell(NTmaps(:,1)));
    else
        VDJdata(:,GroupNumLoc) = num2cell(NTmaps(:,1));
    end
    
    if strcmpi(Option','shmham') %Divide by 2
        MaxDev = MaxDev/2;
    end

    %Save the sorted file
    DotLoc = find(FileName == '.');
    OutputFilePre = [FileName(1:DotLoc(end)) 'ham' num2str(MaxDev)]; %Add hamming distance info
    if ispc %Write to Excel if it's a pc
        xlswrite([FilePath OutputFilePre '.xlsx'],cat(1,NewHeader,VDJdata));
    else %Write to csv, tab-delimited file
        writeDlmFile(cat(1,NewHeader,VDJdata),[FilePath OutputFilePre '.csv'],'\t');
    end
end