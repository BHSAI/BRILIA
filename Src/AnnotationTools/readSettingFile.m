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
%Extract the setting that exists
FID = fopen(FileName,'r');
if FID < 0
    error('Could not open setting file')
end

%Get the parsed P input, which stores the current file names
if ~isempty(varargin)
    P0 = varargin{1};
else
    P0 = BRILIA('getinput');
end

%Remove the field name search for file names, as you want the inpuuted one.
SettingNames = fieldnames(P0);
DelLoc = findCell(SettingNames,'FullFileNames'); %Do no mess with file names from setting files
if DelLoc(1) > 0
    SettingNames(DelLoc) = [];
end

%Read the setting file and fill in P based on setting file data
P = P0; %The final setting structure
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

        Value = S(MatchStart:MatchEnd);
        Value = strrep(Value,'''',''); %Replace ' from strings       
        try  %Try converting to numerical value
            Value = eval(Value);
        catch %If you can't, just remove the spaces from here
            Value = strrep(Value,' ','');
        end        
        P.(SettingNames{j}) = Value;
        
        break
    end
end
fclose(FID);

%Convert field value of numbers to numbers
for j = 1:length(SettingNames)
    Value = P.(SettingNames{j});
    try
        P.(SettingNames{j}) = eval(Value);
    catch
    end
end
