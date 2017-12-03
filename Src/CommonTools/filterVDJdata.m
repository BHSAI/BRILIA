%filterVDJdata will search through VDJdata for certain information,
%returning a filtered VDJdata cell matrix of all clusters that have at
%least one member that matched the filter criteria. This function will also
%return a binary index of size(VDJdata,1)x1 that can be used for more
%complicated AND + OR filtering of VDJdata.
%
%  FiltVDJdata = filterVDJdata(HeaderName,Value,...)
%
%  FiltVDJdata = filterVDJdata(HeaderLoc,Value,...)
%
%  FiltVDJdata = filterVDJdata(VDJdata,VDJheader,...)
%
%  FiltVDJdata = filterVDJdata(FileName,...)
%
%  FiltVDJdata = filterVDJdata(...,Logic)
%
%  [FiltVDJdata, GrpBinIdx] = filterVDJdata(...)
%
%  [FiltVDJdata, GrpBinIdx, SeqBinIdx] = filterVDJdata(...)
%
%  INPUT
%    HeaderName: is the name of the VDJdata data column of VDJdata. Can
%      also use column index instead (see next)
%    HeaderLoc: is the nth column of VDJdata to search. It's related to
%      HeaderName by HeaderName = VDJheader{HeaderLoc}. Can also use name
%      instead.
%    Value: is either a string, list {n1,n2}, range [min max], or number N
%      that is being sought in a data column of VDJdata.
%    Logic ['AND','OR']: uses either the AND or OR logic for all the filter
%      parameters. Can be specified anywhere in the input, though try to
%      make it before the HeaderName-Value pairs.
%
%    Below is the list of HeaderName that can be used. Value is either
%      numeric or character that is being sought. Certain rules apply to
%      different value types as listed below.
%        alphanum: can match any case, including partial matches
%        num: list is inputted as a cell like {N1,N2,N3}
%             range is inputted as a matrix like [min,max]
%             single number is just an integer value like N
%        seq: does a same-length string comparison, but X is a wildcard
%        prop: similar to seq, but uses AA property codes which are
%          different from AA letter. AA with similar properties gets a
%          single letter code. See the AminoAcidProp.csv file. Unrecognized
%          property code letters are wildcard.
%
%      HeaderName         Value type    Special Rules
%      ---------------    ------------  -----------------------------------
%      SeqName            alphanum      -
%      SeqNum             num           -
%      GroupNum           num           -
%      TemplateCount      num           -
%      TreeChildCount     num           -
%
%   *  H-CDR3_Property    prop          will look for CDR3 data, translate
%                                       them to AA property codes, and then
%                                       search with query property code.
%      H-CDR3_AminoAcid   seq           -
%      H-CDR3_Length      num           -
%      H-Seq              seq           -
%      H-RefSeq           seq           -
%      H-Functional       alphanum      can be 'y','n','m'
%      H-V_GeneName       alphanum      -
%      H-V_MapNum         num           -
%      H-V_Deletion3      num           -
%      H-D_GeneName       alphanum      -
%      H-D_MapNum         num           -
%      H-D_Deletion5      num           -
%      H-D_Deletion3      num           -
%      H-J_GeneName       alphanum      -
%      H-J_MapNum         num           -
%      H-J_Deletion5      num           -
%      H-Length_V         num           -
%      H-Length_Nvd       num           -
%      H-Length_D         num           -
%      H-Length_Ndj       num           -
%      H-Length_J         num           -
%      H-SHM_V            num           -
%      H-SHM_Nvd          num           -
%      H-SHM_D            num           -
%      H-SHM_Ndj          num           -
%      H-SHM_J            num           -
%
%   *  L-CDR3_Property    prop          will look for CDR3 data, translate
%                                       them to AA property codes, and then
%                                       search with query property code.
%      L-CDR3_AminoAcid   seq           -
%      L-CDR3_Length      num           -
%      L-Seq              seq           -
%      L-RefSeq           seq           -
%      L-Functional       alphanum      can be 'y','n','m'
%      L-V_GeneName       alphanum      -
%      L-V_MapNum         num           -
%      L-V_Deletion3      num           -
%      L-J_GeneName       alphanum      -
%      L-J_MapNum         num           -
%      L-J_Deletion5      num           -
%      L-Length_V         num           -
%      L-Length_N         num           -
%      L-Length_J         num           -
%      L-SHM_V            num           -
%      L-SHM_N            num           -
%      L-SHM_J            num           -
%
%  OUTPUT
%    FiltVDJdata: filtered VDJdata containing groups that have matched the
%      Header-Value pair by at least one group member.
%    GrpBinIdx: Binary index of size size(VDJdata,1)x1 of VDJdata groups
%      that were extracted.
%    SeqBinIdx: Binary index of size size(VDJdata,1)x1 of VDJdata sequences
%      that triggered match.
%
%  NOTE
%    If nothing matches, will return an empty cell and all BinIdx outputs
%    will be zeros(size(VDJdata,1),1,'logical');
%
%    This uses the AND logic by default. To combine AND and OR, must use
%    filterVDJdata multiple times to get the desired results, and then
%    properly use the binary index outputs to get the right results.
%
%    Range values are always inclusive, and placed in brackets. Example: To
%    find a range within 0 to 100, set filter value to [0,100]. To get
%    greater than 100, set filter value to [100.01,inf].
%
%    List values are always in a cells. Example: To get item number 1,2,4,
%    you set the filter Value to {1,2,3}.
%
%  EXAMPLE
%    Finding groups with
%    TemplateCount >= 10 AND CDR3seq = 'CARXSYW'
%
%      [VDJdata,VDJheader] = openSeqData; %Opens the BRILIA seq file.
%      H = getHeaderVar(VDJheader); %Gets structure of main column
%      [filtVDJdata,GrpBinIdx,SeqBinIdx] =
%        filterVDJdata(VDJdata,VDJheader,'AND',VDJheader{H.TemplateLoc},[10
%        Inf],'CDR3_AminoAcid','CARXSYW')
%
%    Finding Sequences with:
%    TemplateCount >= 10 AND (CDR3seq = 'CARXSYW' OR CDR3seq = 'CSRFXYW')
%
%      [~,~,SeqBinIdx1] = 
%        filterVDJdata(VDJdata,VDJheader,'TemplateCount',[10,Inf]);
%      [~,~,SeqBinIdx2] = 
%        filterVDJdata(VDJdata,VDJheader,'CDR3_AminoAcid','CARXSYW');
%      [~,~,SeqBinIdx3] = 
%        filterVDJdata(VDJdata,VDJheader,'CDR3_AminoAcid','CSRFXYW');
%      SeqBinIdx = SeqBinIdx1 & (SeqBinIdx2 | SeqBinIdx3);
%      filterVDJdata = VDJdata(SeqBinIdx,:);

function varargout = filterVDJdata(varargin)
%See if user gave VDJdata, VDJheader  OR  FileName
if nargin < 3
    error('Not enough inputs');
else
    if iscell(varargin{1}) && iscell(varargin{2}) %VDJdata and VDJheader were given
        VDJdata = varargin{1};
        VDJheader = varargin{2};
    elseif ischar(varargin{1}) %has to be char, but can't tell if it's a file name or header name
        if ~isempty(varargin{1} == '.') %Filenames have periods
            FileName = varargin{1};
            try
                [VDJdata,VDJheader] = openSeqData(FileName);
            catch
                [VDJdata,VDJheader] = openSeqData;
            end
            varargin(1) = [];
        else %Need to select that file
            [VDJdata,VDJheader] = openSeqData;
            if isempty(varargin{1})
                varargin(1) = [];
            end
        end
    else
        error('Check inputs');
    end
end

%Look now for presence of logic input.
Logic = 'AND';
LogicLoc = findCell(varargin,{'and','or'},'MatchCase','any');
if ~isempty(LogicLoc) && LogicLoc(1) > 0
    Logic = upper(varargin{LogicLoc(1)});
    varargin(LogicLoc) = [];
end

%Whatever is left of varargin are Param-Value pairs.
if mod(length(varargin),2) ~= 0
    error('%s: Mismatched number of param-value pairs.',mfilename);
end

H = getHeavyHeaderVar(VDJheader);

%Just reorganize inputs into Nx3 cell of {HeaderName HeaderLoc FilterValue DataType}
InputCell = cell(length(varargin)/2,4);
q = 1;
for j = 1:2:length(varargin)
    if isnumeric(varargin{j}) %Provided HeaderLoc
        InputCell{q,2} = varargin{j};
        InputCell{q,1} = VDJheader{varargin{j}};
    elseif ischar(varargin{j}) %Provided HeaderName
        InputCell{q,1} = varargin{j};
        InputCell{q,2} = findCell(VDJheader,varargin{j});
    end
    InputCell{q,3} = varargin{j+1};
    InputCell{q,4} = getDataType(InputCell{q,1});

    %Special case for CDR3_Property, set the column location to CDR3Loc(1)
    if strcmpi(InputCell{q,1},'CDR3_Property')
        InputCell{q,2} = H.CDR3Loc(1);
    end

    q = q + 1;
end


%Now go through parameters and extract the selected groups
if strcmp(Logic,'AND')
    SeqBinIdx = ones(size(VDJdata,1),1,'logical'); %Start with ones
else
    SeqBinIdx = zeros(size(VDJdata,1),1,'logical'); %Start with zeros
end

%For each critiria, determine valid sequence location in VDJdata
for j = 1:size(InputCell,1)
    %Check if there is a data column to search
    DataColNum = InputCell{j,2};
    if isempty(DataColNum) || DataColNum < 1 || DataColNum > size(VDJdata,2) %Invalid number
        continue
    end
    
    %Prepare data type and data value for searching
    DataValue = InputCell{j,3};
    DataType = InputCell{j,4};    
    switch DataType
        case 'seq'
            BinIdxT = searchSeq(VDJdata(:,DataColNum),DataValue);
        case 'num'
            BinIdxT = searchNumeric(VDJdata(:,DataColNum),DataValue);
        case 'alphanum'
            BinIdxT = searchAlphaNum(VDJdata(:,DataColNum),DataValue);
        case 'prop'
            BinIdxT = searchProp(VDJdata(:,DataColNum),DataValue);
    end
    
    %Save the cumulative results
    if strcmp(Logic,'AND')
        SeqBinIdx = SeqBinIdx & BinIdxT;
    else strcmp(Logic,'OR')
        SeqBinIdx = SeqBinIdx | BinIdxT;
    end    
end

%Now go and fill in the group binary idx, which differs from seq binary idx
if H.GrpNumLoc > 0
    GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
    UnqGrpNum = unique(GrpNum(SeqBinIdx));
    GrpBinIdx = zeros(size(VDJdata,1),1,'logical');
    for y = 1:length(UnqGrpNum)
        GrpBinIdx = GrpBinIdx | (UnqGrpNum(y) == GrpNum);
    end
else
    GrpBinIdx = SeqBinIdx;
end

if nargout >= 1
    varargout{1} = VDJdata(GrpBinIdx,:);
    if nargout >= 2
        varargout{2} = GrpBinIdx;
        if nargout >= 3
            varargout{3} = SeqBinIdx;
        end
    end
end

function BinIdxT = searchSeq(M,QuerySeq)
%Prepare query
QuerySeq = upper(QuerySeq);
XlocQ = (QuerySeq == 'X');

%Search for match to query
BinIdxT = zeros(size(M),'logical');
for j = 1:length(M(:))
    if ~ischar(M{j}); continue; end
    if length(M{j}) == length(QuerySeq)
        SeqT = upper(M{j});
        XlocT = SeqT == 'X';
        MatchLoc = (SeqT == QuerySeq) | XlocQ | XlocT;
        if min(MatchLoc) == 1
            BinIdxT(j) = 1;
        end
    end
end

function BinIdxT = searchNumeric(M,QueryNum)
%Prepare query
if iscell(QueryNum) %Have a list
    while iscell(QueryNum) %Unwrap a cell, until it is no longer one
        QueryNumT = QueryNum{1};
        if ~iscell(QueryNumT) 
            break;
        else
            QueryNum = QueryNumT;
        end
    end
    QueryNum = cell2mat(QueryNum);
    NumType = 1; %List of single values
else
    QueryNum = sort(QueryNum);   
    if length(QueryNum) > 1 && QueryNum(1) < QueryNum(2)
        NumType = 2; %Range
    else
        NumType = 1; %Single value
    end
end

%Convert M to matrix, but only after ensuring no empty or characters
ValidIdx = ones(size(M),'logical'); %Locations where there are number, and not empties
for j = 1:length(M(:))
    if ~isnumeric(M{j});
        M{j} = 0;
        ValidIdx(j) = 0;
    elseif isempty(M{j})
        M{j} = 0;
        ValidIdx(j) = 0;
    end
end
M = cell2mat(M);

%Search for match to query
if NumType == 1 %Single value(s)
    BinIdxT = zeros(size(M),'logical');
    for k = 1:length(QueryNum)
        BinIdxT = BinIdxT | (M == QueryNum(k));
    end
elseif NumType == 2 %Range value
    BinIdxT = (QueryNum(1) <= M) & (M <= QueryNum(2));
end

%Remember to account for empty M entries
BinIdxT = BinIdxT & ValidIdx;

function BinIdxT = searchAlphaNum(M,QueryAlphaNum)
%Search for match to query
BinIdxT = zeros(size(M),'logical');
FindMatchIdx = findCell(M,QueryAlphaNum,'MatchCase','any','MatchWord','partial');
if ~isempty(FindMatchIdx) && FindMatchIdx(1) ~= 0
    BinIdxT(FindMatchIdx) = 1;
end

function BinIdxT = searchProp(M,QueryProp)
%Prepare Query
M = getPropertyCode(M); %Assumes M is already a cell of CDR3 seq

%Search for match to query
BinIdxT = zeros(size(M),'logical');
for j = 1:length(M(:))
    if length(M{j}) == length(QueryProp)
        SeqT = upper(M{j});
        MatchLoc = (SeqT == QueryProp);
        if min(MatchLoc) == 1
            BinIdxT(j) = 1;
        end
    end
end
