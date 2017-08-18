%readSettingFile will read a text file containing setting required for
%BRILIA to operate. 
%
%  function readSettingFile(FileName,P)
%
%  INPUT
%    FileName: File name of the settings text file to save
%    P: Default settings. If empty, will use it's own defaults as shown.
%      P.Species = '';
%      P.Strain = '';
%      P.Ddirection = 'all';
%      P.Vfunction = 'all';
%      P.DevPerc = 5;
%      P.FileType = '';
%      P.Delimiter = ';';
%      P.CheckSeqDir = 'n';
%      P.Chain = 'H';
%
%  OUTPUT
%    A text file of setting with the name FileName (with .txt extension).
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
%       Chaing = 'HL';
%     End of [FileName].txt------------------------------
%
%  See also readSettingFile

function makeSettingFile(FileName,P)
%Make sure P is a structure
if ~isstruct(P)
    return
end

%Setting up default structure
if isempty(P)
    P = BRILIA('getinput');
end

%Extract the setting that exists
[FID,Msg] = fopen(FileName,'w');
if FID < 0
    disp(Msg);
    return;
end

SettingNames = fieldnames(P);
%Write the date first
fprintf(FID,'%s = ''%s'';\r\n','Date',date);
for j = 1:length(SettingNames)
    Name = SettingNames{j};
    Value = P.(Name);
    
    if ~isnumeric(Value) && ~ischar(Value)
        continue; %Probably some handle object
    elseif isnumeric(Value); %Convert to string if its a matrix or number
        if length(Value) > 1
           Value = mat2str(Value);
        else
           Value = num2str(Value);
        end
    end
    fprintf(FID,'%s = ''%s'';\r\n',Name,Value);
end

fclose(FID);
