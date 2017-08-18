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
    warning('%s: Input a cell of file name strings. No action taken.', mfilename);
    return;
elseif isstruct(FileList) %Convert to a cell of file name strings
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

%Determine if you want to append to existing output file
AppendThis = 'n';
if ~isempty(varargin) && ischar(varargin{1}) && strcmpi(varargin{1}, 'append')
    AppendThis = 'y';
end

%Determine if you must delete old file AFTER you have the OK to do so. 
if AppendThis == 'n'
    if exist(OutputFile, 'file') > 0
        delete(OutputFile);
    end
end

%See if the folder exists
[OutFilePath, ~, ~] = parseFileName(OutputFile);
if ~exist(OutFilePath, 'dir')
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
