%readSettingFile will read a text file containing setting required for
%BRILIA to operate. 
%
%  function P = readSettingFile(FileName,P)
%
%  INPUT
%    FileName: File name of the settings text file (see Example below).
%    P: Default settings. If empty, will use it's own defaults as shown.
%      P.Species = '';
%      P.Strain = '';
%      P.Ddirection = 'all';
%      P.Vfunction = 'all';
%      P.DevPerc = 3;
%      P.FileType = '';
%      P.Delimiter = ';';
%      P.CheckSeqDir = 'n';
%
%  OUTPUT
%    P: structure file containing all the settings required of BRILIA.
%
%  EXAMPLE SETTING FILE
%     Start of [FileName].txt----------------------------
%       Date = '27-Jan-2017';
%       Species = 'mouse';
%       Strain = 'all';
%       Ddirection = 'all';
%       Vfunction = 'all';
%       DevPerc = '5';
%       FileType = 'delimited';
%       Delimiter = ';';
%       CheckSeqDir = 'n';
%     End of [FileName].txt------------------------------
%
%  See also makeSettingFile

function P = readSettingFile(FileName,varargin)
%Read the setting file 
FID = fopen(FileName, 'r');
assert(FID > 0, '%s: Could not open setting file "%s".', mfilename, FileName);
T = textscan(FID, '%s', 'delimiter', '\n');
fclose(FID);

%Get current parsed input
if ~isempty(varargin)
    P0 = varargin{1};
else
    P0 = BRILIA('getinput');
end

%Clean setting file to get setting string and value
S = cellfun(@(x) regexpi(x, '\s*=\s*', 'split'), T{1}, 'UniformOutput', false);
Sstr = cellfun(@(x) x{1}, S, 'UniformOutput', false);
Sval = cellfun(@(x) strrep(strrep(x{2}, '''', ''), ';', ''), S, 'UniformOutput', false);

%Backward compatibility for version before v4.0.0
Sstr = regexprep(Sstr, 'Vfunction', 'Vgene', 'ignorecase');
Sstr = regexprep(Sstr, 'Ddirection', 'Dgene', 'ignorecase');
DevPercLoc = ~startsWith(Sstr, 'DevPerc', 'ignorecase', true);
Sstr = Sstr(DevPercLoc);
Sval = Sval(DevPercLoc);

SettingNames = fieldnames(BRILIA('getinput'));
[~, Idx, ~] = intersect(lower(Sstr), lower(SettingNames));

%Save new values to P structured input
P = P0;
for k = 1:length(Idx)
    try
        NumVal = convStr2Num(Sval{Idx(k)});
        if ~isempty(NumVal)
            P.(Sstr{Idx(k)}) = NumVal;
            continue
        end
    catch
    end
    P.(Sstr{Idx(k)}) = Sval{Idx(k)};
end