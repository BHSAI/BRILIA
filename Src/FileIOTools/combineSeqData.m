%combineSeqFile will combine multiple BRILIA output files into one large
%file.
%
%  combineSeqData(FileList, OutputFile)
%
%  combineSeqData(FileList, OutputFile, 'append')
%
%  INPUT
%    FileList: Nx1 cell or struct of full file names to combine
%    OutputFile: full file name to write file
%    'append': append to exising output file. Otherwise, overwrite.
%
%  OUTPUT
%    Outputs a combined BRILIA data file, saved according to the OutputFile
%      name.

function combineSeqData(FileList, OutputFile, varargin)
%Detemine if there are multiple files to combine
if ischar(FileList) 
    if ~contains(FileList, '\*')
        FileList = dir(FileList);
    else
        warning('%s: No action taken.', mfilename);
        return
    end
end
if isstruct(FileList) %Convert to a cell of file name strings
    if length(FileList) <= 1
        warning('%s: Need more than one file. No action taken.', mfilename);
        return;
    else
        FileListT = FileList;
        FileList = cell(length(FileListT),1);
        for j = 1:length(FileListT)
            FileList{j} = FileListT(j).name;
        end
    end
elseif iscell(FileList) && length(FileList) <= 1
    warning('%s: Need more than one file. No action taken.', mfilename);
    return;
end

%Make sure OutputFile is not the same as one of the InputFile
FileNames = cell(size(FileList));
for k = 1:length(FileNames)
    [~, FileNames{k}, ~] = fileparts(FileList{k});
end
[~, OutName, ~] = fileparts(OutputFile);
if ismember(OutName, FileNames)
    warning('%s: The output "%s" cannot also be in the input files.', mfilename, OutputFile);
    Attempt = 0;
    while 1
        Ask = input('Exclude the output file from the input file and combine the rest to output?\n   y or n :', 's');
        if strcmpi(Ask, 'n')
            fprintf('%s: Exiting due to cancelation.\n', mfilename)
            return
        elseif strcmpi(Ask, 'y')
            fprintf('%s: Continuing to combine rest to output.\n', mfilename)
            FileList(ismember(FileNames, OutName)) = [];
            break
        else
            Attempt = Attempt + 1;
            if Attempt >= 5
                error('%s: Exiting due to invalid option', mfilename);
            end
        end
    end
end


%Determine if you want to append to existing output file
AppendThis = 'n';
if ~isempty(varargin) && ischar(varargin{1}) && contains(varargin{1}, 'append')
    AppendThis = 'y';
end

%Determine if you must delete old file AFTER you have the OK to do so. 
if AppendThis == 'n'
    if exist(OutputFile, 'file') > 0
        delete(OutputFile);
    end
end

if isempty(OutputFile)
    error('%s: No output file specified', mfilename);
end

%See if the folder exists
[OutFilePath, ~, ~] = parseFileName(OutputFile);
if ~isempty(OutFilePath) && ~exist(OutFilePath, 'dir')
    mkdir(OutFilePath)
end

%Append/create file
if strcmpi(AppendThis, 'y')
    Mode = 'a';
else 
    Mode = 'w';
end
FID_0 = fopen(OutputFile, Mode);
if FID_0 < 0
    error('%s: Could not create file %s', mfilename, OutputFile)
end

for f = 1:length(FileList)
    %fprintf('Combining %d of %d files.\n', f, length(FileList));
    try
       FID_1 = fopen(FileList{f});
       if f > 1
           fgetl(FID_1); %Skip the first line header
       end
       text = fgetl(FID_1);
       while ischar(text)
           fprintf(FID_0, '%s\n', text);           
           text = fgetl(FID_1);
       end
       fclose(FID_1);
    catch
       warning('%s: Could not add this file: %s', mfilename, FileList{f});
    end
end
fclose(FID_0);
