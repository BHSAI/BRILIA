%refactorVar will search through all *.m files in a directory and replace
%variable VAR1 to VAR2, saving a backup .m file is set by default, as this
%is not undoable.
%
%  refactorVar(VAR1,VAR2...) will replace variables VAR1 to VAR2 in current
%  direcotry
%
%  refactorVar(VAR1,VAR2,'InDir',DIR,'Backup','y','Subfolder','y') will
%  search inside the specficed DIR directory for refactoring, and will
%  backup the file before making the edits. Backed up files are saved as
%  'BACKEDUP_[FILENAME].m', which helps to prevent overlapping function
%  calls with the new files.
% 
%  Use function "undoRefactor" to undo refactorization, BUT this works only
%  if you used the backup feature. 
%
%  Note: Only full word exact matches are allowed.
%        EX:
%           VAR1 = 'makeString'
%           VAR2 = 'makeStr'
%           If a code line has the following string.
%           S = 'makeString, (makeString), [makeString], {makeString}, A.makeString, A.makeString.B';
%           The new code will have the following new string.
%           S = 'makeStr, (makeStr), [makeStr], {makeStr}, A.makeString, A.makeString.B'
function refactorVar(VAR1,VAR2,varargin)
%Ask user about warning
ConfirmThis = input('WARNING: This will replace ALL variables. Are you sure? y or n: ','s');
if ~isempty(ConfirmThis)
    if ~strcmpi(ConfirmThis(1),'y')
        return
    end
elseif isempty(ConfirmThis)
    return
end

%Prepare the parameter inputs
P = inputParser;
addParameter(P,'InDir','',@ischar);
addParameter(P,'Backup','yes',@ischar);
addParameter(P,'Subfolder','no',@ischar);
parse(P,varargin{:});

InDir = P.Results.InDir;
if isempty(InDir)
    InDir = cd;
end
SlashLoc = regexp(InDir,'\\|\/');
SlashType = InDir(SlashLoc(end));
if ~strcmp(InDir(end),SlashType)
    InDir = [cd SlashType];
end

%Get any subfolders too
if strcmpi(P.Results.Subfolder(1),'y')
    SubDir = genpath(InDir);
    InDir = regexp(SubDir,';','split')';
    if isempty(InDir{end})
        InDir(end) = [];
    end
    %Add the '\'
    for j = 1:length(InDir)
        if ~strcmp(InDir{j}(end),SlashType)
            InDir{j} = [InDir{j} SlashType];
        end
    end
else
    InDir = {InDir};
end

for d = 1:length(InDir)
    %Obtain the *.m file lists
    ListMfiles = dir([InDir{d} '*.m']);
    DelThis = zeros(length(ListMfiles),1,'logical');
    for f = 1:length(ListMfiles)
        if strfind(ListMfiles(f).name,'BACKUP_')
            DelThis(f) = 1;
        end
    end
    ListMfiles(DelThis) = [];

    %Make a backup of the file lists
    if strcmpi(P.Results.Backup(1),'y')
        BackupPre = 'BACKUP_';
        for f = 1:length(ListMfiles)
            copyfile([InDir{d} ListMfiles(f).name],[InDir{d} BackupPre ListMfiles(f).name],'f')
        end
    end

    %For each file, open and replace and save;
    for f = 1:length(ListMfiles)
        FID = fopen([InDir{d} ListMfiles(f).name],'r');
        FileStr = cell(1000,1); %Initialize cell;
        t = 1;
        while feof(FID) == 0
            if t > size(FileStr,1)
                FileStr = [FileStr; cell(1000,1)]; %Add more;
            end
            GetText = fgetl(FID);
            GetText = repStrCode(GetText,VAR1,VAR2);
            FileStr{t,1} = GetText;
            t = t+1;
        end
        fclose(FID);
        FileStr(t:end) = []; %Delete unused ones.

        FID = fopen([InDir{d} ListMfiles(f).name], 'w');
        for i = 1:length(FileStr)
            fprintf(FID,'%s\n', FileStr{i});
        end
        fclose(FID);
    end
end
