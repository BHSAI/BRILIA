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
%Setting up default structure
if ~isempty(varargin) && isstruct(varargin{1}) 
    P = varargin{1};
else
    P = BRILIA('getinput');
end

%Extract the setting that exists
FID = fopen(FileName,'r');
if FID < 0
    error('Could not open setting file')
end

SettingNames = fieldnames(P);
while feof(FID) == 0
    S = fgetl(FID);
    for j = 1:length(SettingNames)
        SearchPat = [SettingNames{j} '[\s*]='];
        MatchStart = regexpi(S,SearchPat,'end')+1;
        if ~isempty(MatchStart)
            MatchEnd = regexpi(S,';')-1;
            MatchEnd = MatchEnd(end);
        else
            continue
        end

        ExtractStr = S(MatchStart:MatchEnd);
        ExtractStr = strrep(ExtractStr,' ','');
        ExtractStr = strrep(ExtractStr,'''','');
        
        P.(SettingNames{j}) = ExtractStr;
        break
    end
end
fclose(FID);

%Convert field value of numbers to numbers
for j = 1:length(SettingNames)
    Value = P.(SettingNames{j});
    if ischar(Value) && ~isempty(Value)
        %Look for non-integer numbers
        NotNumLoc = regexpi(Value,'[^\d]');
        if isempty(NotNumLoc)
            P.(SettingNames{j}) = eval(Value);
        end
    end
end