%sortAAseq takes a trimmed, aligned group of sequences from an Excel file
%and rearranges it based on how similar the AA sequences are from each
%other. Will group first based on AA length, and then hamming distances. 
%
%  Output = sortAAseq()   This will ask user to select a file to open,
%  assuming Column 1 = Seq. The outputfile is an Excel file with these
%  data, plus a last column with the group number. Similarity is determined
%  as 6 NTs different between members within the same family.
%
%  Output = sortAAseq(6,AAseq) will process the input of sequence, grouping
%  Hamming distances of 6 or less into same group.
%
%  Output is [GrpNum SeqNum], sorted by GrpNum.
%  
%  NOTE! Consecutive mismatches causes a M^2 hamming distance increases.
%    EX:    001100 is a hamming distance of 4, NOT 2.

function AAmaps = sortAAseq(varargin)

%Extract the inputs
    %Case 1) Maximum deviation is given, but no file name
    %Case 2) Maximum deviation + file name given
    %Case 3) Maximum deviation + AAseq given

%Looking for Maximum deviation info, if none, ask the user. 
MaxDevPerc = 10; %Percentage deviation. It can be high, as we are looking for point AA mutations.
if nargin == 0
    MaxDevPerc = input('What is the Hamming distance % to cluster?');
    if isempty(MaxDevPerc)
        MaxDevPerc = 1;
    else
        MaxDevPerc = round(abs(MaxDevPerc));
        fprintf('Hamming Distance is: %d \n', MaxDevPerc')
    end
else
    DelThis = zeros(1,length(varargin)) > 1;
    for v = 1:length(varargin)
        if isnumeric(varargin{v})
            MaxDevPerc = round(abs(varargin{v}));
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
    AAseq = VDJdata(:,CDR3Loc(1));
else
    if iscell(varargin{1}) %Probably sequences
        AAseq = varargin{1};
    else
        [VDJdata, NewHeader, FileName, FilePath] = openSeqData(varargin{v});
        getHeaderVar;
        AAseq = VDJdata(:,CDR3Loc(1));
    end
end

%Determine the CDR3 lengths
SortMat = zeros(size(AAseq,1),3); %[OrigSeqNum, CDR3Length, UnqClustNum]
for j = 1:size(SortMat,1)
    SortMat(j,1) = j;
    SortMat(j,2) = length(AAseq{j});
end
SortMat = sortrows(SortMat,2);

%Now divide each cluster, and subcluster
UnqAALen = unique(SortMat(:,2));
UnqAALen(UnqAALen == 0) = [];
GrpCt = 0;
for y = 1:length(UnqAALen)
    IdxLoc = SortMat(:,2) == UnqAALen(y);
    IdxLocAbs = SortMat(IdxLoc,1);
    AAseqT = AAseq(IdxLocAbs);
    
    %CONSIDER UNIQUE AAseqT, then do this. Reduces comp time. Need to remap
    %the unq seq.
    
    %Calculate matrix of NT differences
    ScoreMat = calcPairDist(AAseqT,'hampen'); %Hamming distance with M^2 penalty on consec mismatches, M.
    
    %Determine the clusters
    MaxDev = ceil(MaxDevPerc/100*UnqAALen(y));
    AncMap = calcAncMap(ScoreMat);
    AncMap(AncMap(:,3) > MaxDev, 2) = 0;
    AAmapT = findTreeClust(AncMap);

    %Update SortMat with unique grp numbers
    SortMat(IdxLoc,3) = AAmapT(:,2) + GrpCt;
    GrpCt = max(SortMat(IdxLoc,3));    
end
    
%Save this clustering, which will be used to extract sequences individual
%files and identify it's sharedness.
AAmaps = SortMat(:,[3 1]); %Want GrpNum | SeqNum