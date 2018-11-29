%formatSeqMut will convert the format of the getMutInfo data structure into
%either the scalar or nonscalar version.
%
%  MutInfo = formatMutInfo(MutInfo, Format)
%
%  INPUT
%    MutInfo: output from getMutInfo
%    Format ['cell', 'struct']: desired format of the output. If empty,
%      will switch between the two structure (cell <-> struct)
%
%  OUTPUT ('stuct' option)
%    MutInfo: nonscalar structure of each MUTATED nt
%      .Idx    - the nth position of nt mutation
%      .NT     - char like 'x>y' indicating nucleotide x to y mutation
%      .AA     - char like 'X>Y' indicating amino acid X to Y mutation
%      .Motif  - 5-letter motif around mutant (X = unknown/filler) nnXnn
%      .IsSyn  - 1 or 0 for synonymous mutations (requires Frame input)
%
%  OUTPUT ('cell' option)
%    MutInfo: scalar structure
%      .Header - 1xN cell of data header strings
%      .Data   - MxN cell of data
%
%  NOTE
%    If appending new fields to MutInfo 'struct' format, make sure it's a
%    single-cell value and not an array, such that you can convert between
%    the two formats.


function Out = formatMutInfo(MutInfo, Format)
CurFormat = getCurFormat(MutInfo);
if nargin == 1 %Need to determine Format as the other one
    switch CurFormat
        case 'cell'; Format = 'struct';
        case 'struct'; Format = 'cell';
    end
end

if ~any(startsWith({'cell', 'struct'}, Format, 'ignorecase', true))
    error('%s: The Format input is invalid. Use ''cell'' or ''struct''.', mfilename);
end

if startsWith('cell', Format, 'ignorecase', true)
    if strcmpi(CurFormat, 'cell')
        Out = MutInfo;
    else
        Out.Header = fieldnames(MutInfo)';
        Out.Data = squeeze(struct2cell(MutInfo))';
    end
elseif startsWith('struct', Format, 'ignorecase', true)
    if strcmpi(CurFormat, 'struct')
        Out = MutInfo;
    else
        Out = cell2struct(MutInfo.Data, MutInfo.Header, 2);
    end
end

function CurFormat = getCurFormat(MutInfo)
if ~isstruct(MutInfo)
    error('%s: MutInfo input must be a structure', mfilename);
end
if all(isfield(MutInfo, {'Header', 'Data'}))
    CurFormat = 'cell';
elseif all(isfield(MutInfo, {'Idx', 'NT', 'Motif'}))
    CurFormat = 'struct';
else
    error('%s: MutInfo does not have essential fields ''Idx'', ''NT'', ''Motif''.', mfilename);
end