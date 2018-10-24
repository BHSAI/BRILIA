%getVDJmapper returns a structure where the fields refer to column
%number(s) corresponding to a data type in VDJdata cell array. It uses the
%VDJheader information.
%
%  [Map, NumIdx, CharIdx] = getVDJmappper(VDJheader)
%
%  INPUT
%    VDJheader: 1xN cell of header string for VDJdata
%
%  OUTPUT
%    Map: structure of standard commonf field names corresponding to
%      VDJheader index or indices (Mx1 matrix). Values are 0 for
%      non-existing fields in VDJheader.
%    NumIdx: indices of VDJdata column that are numeric
%    CharIdx: indices of VDJdata column that are strings
%
function [Map, NumIdx, CharIdx] = getVDJmapper(VDJheader)
persistent C CType H HType L LType %No need to recreate static variables 

if isstruct(VDJheader) %Keep this here for backward compatibility issues
    Map = VDJheader;
    return
end

%Setup the static variables once, if not initialized
[C, CType] = COMMON_CHAIN_VARNAME(C, CType); 
[H, HType] = HEAVY_CHAIN_VARNAME(H, HType);
[L, LType] = LIGHT_CHAIN_VARNAME(L, LType);

%Create the merged structure (workaround to Matlab's inability to do this)
Header = ['Chain' fieldnames(C)',  fieldnames(H)',  fieldnames(L)'];
Values = ['XX'    struct2cell(C)', struct2cell(H)', struct2cell(L)'];
Map = struct;
for j = 1:length(Header)
    Map.(Header{j}) = Values{j};
end

%Look for the indices in VDJheader for each column name
Map = replaceStrWithIdx(Map, VDJheader);

%Fill in the Map.Chain information, which is used by other functions
Map.Chain = '';
if Map.hSeq > 0; Map.Chain(end+1) = 'H'; end
if Map.lSeq > 0; Map.Chain(end+1) = 'L'; end 

%Determine indices of numeric and string columns of VDJdata, based on VDJheader
if nargout >= 2
    %Determine the data types, 1 = numeric, 0 = string
    Types = [struct2cell(CType); struct2cell(HType); struct2cell(LType)];
    Types = vertcat(Types{:});

    %Extract the current position of non-zero indicies of all of Map
    MapIdx = struct2cell(Map);
    ChainLoc = strcmpi(fieldnames(Map), 'chain');
    MapIdx{ChainLoc} = []; %This is for Map.Chain
    MapIdx = vertcat(MapIdx{:});
    
    %Determine numeric and/or character indices
    NumLoc = repelem(false, 1, numel(VDJheader));
    NumLoc(MapIdx(MapIdx > 0)) = Types(MapIdx > 0);
    NumIdx  = find( NumLoc);
    if nargout >= 3
        CharIdx = find(~NumLoc);
    end
end

%EDIT_NOTE: Words that ends with "$" are string values, while all others
%are numbers. Use the "|" if the data header names have changed (due to
%version changes. Ex: "ParentNum|ParNum" or "hSeq|hSequence$").

%Store variables common to both heavy and light chain
%EDIT_NOTE: Edit this for column naming changes
function [C, CType] = COMMON_CHAIN_VARNAME(C, CType)
if ~isempty(C); return; end
C.SeqName    = 'SeqName$';
C.SeqNum     = 'SeqNum';
C.GrpNum     = 'GroupNum';
C.ParNum     = 'ParentNum|ParNum';
C.Template   = 'TemplateCount';
C.ChildCount = 'TreeChildCount';
[C, CType] = formatChainVarname(C, '');

%Store variables specific to heavy chain
%EDIT_NOTE: Edit this for column naming changes
function [H, HType] = HEAVY_CHAIN_VARNAME(H, HType)
if ~isempty(H); return; end
H.hSeq       = 'Seq$';
H.hRefSeq    = 'RefSeq$';
H.hOverSeq5  = 'OverhangSeq5$';
H.hOverSeq3  = 'OverhangSeq3$';
H.hFunct     = 'Functional$';
H.hCDR1      = 'CDR1AminoAcid$  CDR1Length  CDR1Start  CDR1End';
H.hCDR2      = 'CDR2AminoAcid$  CDR2Length  CDR2Start  CDR2End';
H.hCDR3      = 'CDR3AminoAcid$  CDR3Length  CDR3Start  CDR3End';
H.hGeneNum   = 'VMapNum  DMapNum  JMapNum';
H.hGeneName  = 'VGeneName$  DGeneName$  JGeneName$';
H.hLength    = 'LengthV  LengthNvd  LengthD  LengthNdj  LengthJ';
H.hDel       = 'VDeletion3  DDeletion5  DDeletion3  JDeletion5';
H.hVmut      = 'SHMV';
H.hMmut      = 'SHMNvd';
H.hDmut      = 'SHMD';
H.hNmut      = 'SHMNdj';
H.hJmut      = 'SHMJ';
H.hValign    = 'VAlign$';
H.hDalign    = 'DAlign$';
H.hJalign    = 'JAlign$';
[H, HType] = formatChainVarname(H, 'h');

%Store variables specific to light chain
%EDIT_NOTE: Edit this for column naming changes
function [L, LType] = LIGHT_CHAIN_VARNAME(L, LType)
if ~isempty(L); return; end
L.lSeq       = 'Seq$';
L.lRefSeq    = 'RefSeq$';
L.lOverSeq5  = 'OverhangSeq5$';
L.lOverSeq3  = 'OverhangSeq3$';
L.lFunct     = 'Functional$';
L.lCDR1      = 'CDR1AminoAcid$  CDR1Length  CDR1Start  CDR1End';
L.lCDR2      = 'CDR2AminoAcid$  CDR2Length  CDR2Start  CDR2End';
L.lCDR3      = 'CDR3AminoAcid$  CDR3Length  CDR3Start  CDR3End';
L.lGeneNum   = 'VMapNum  JMapNum';
L.lGeneName  = 'VGeneName$ JGeneName$';
L.lLength    = 'LengthV  LengthN  LengthJ';
L.lDel       = 'VDeletion3  JDeletion5';
L.lVmut      = 'SHMV';
L.lNmut      = 'SHMN';
L.lJmut      = 'SHMJ';
L.lValign    = 'VAlign$';
L.lJalign    = 'JAlign$';
[L, LType] = formatChainVarname(L, 'l');

%Used to format the chains above by append "h" or "l" before the string,
%remove the internal string marker "$", and save all string as a Mx1
%cell, where each cell will be replaced by the index in VDJheader
%EDIT_NOTE: This was used to make it easier to edit the above static string
%names, and to keep track of string vs number.
function [X, T] = formatChainVarname(X, Append)
FieldNames = fieldnames(X);
T = X;
for f = 1:numel(FieldNames)
    Str = strcat(Append, strsplit(X.(FieldNames{f})))'; %Make sure strsplit occurs BEFORE strcat, to get hCDR1AminoAcid, hCDR1Length, ...
    X.(FieldNames{f}) = strrep(Str, '$', ''); %Get rid of the string markers
    T.(FieldNames{f}) = ~endsWith(Str, '$');  %Mark where the string markers are
end

%Format string to be same for matching purposes. 
%EDIT_NOTE: Edit this code if the string matching criteria changes
function Str = formatStrSame(Str)
Str = lower(regexprep(Str, '[^a-zA-Z0-9\|]', ''));

%Replaces the string with matrix of index
%EDIT_NOTE: Edit this for search and replace optimization
function X = replaceStrWithIdx(X, Header)
Header = formatStrSame(Header);
Fields = fieldnames(X);
for f = 1:numel(Fields)
    Str = formatStrSame(X.(Fields{f}));
    if ~iscell(Str); continue; end %probably just a Map.Chain;
    X.(Fields{f}) = zeros(numel(Str), 1);
    for s = 1:numel(Str)
        MatchIdx = find(ismember(Header, strsplit(Str{s}, '|')), 1); %note the strsplit to enable degenerate name searches
        if ~isempty(MatchIdx) 
            X.(Fields{f})(s) = MatchIdx;
        end
    end
end