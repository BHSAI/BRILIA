%writeDlmFile will take a cell array and convert it into a string array,
%and then save it into a delimited file. Works better than dlmwrite and
%csvwrite.
%
%  writeDlmFile(CellData,OutputFile,Delimiter)
%
%  INPUT
%    CellData: cell array (but must not have an embeded cell in a cell);
%    OutputFile: full name of the file to be saved to
%    Delimiter [';' ',' '\t]: Separater marker between data. Default ';'.
%
%  See also readDlmFile

function writeDlmFile(CellData,OutputFile,varargin)
%Convert all cell inputs as strings only
for j = 1:size(CellData,1)
    for k = 1:size(CellData,2)
        CurVar = CellData{j,k};
        if isnumeric(CurVar)
            if isnan(CurVar) 
                CellData{j,k} = 0;
            else
                CellData{j,k} = num2str(CurVar);
            end
        elseif isempty(CurVar)
            CellData{j,k} = '';
        end
    end
end

%Create the format
Delimiter = ';';
if length(varargin) == 1
    Delimiter = varargin{1};
end
DataFormat = ['%s' Delimiter];
TxtForm = [repmat(DataFormat,1,size(CellData,2)-1) '%s\n'];

%Save the file
FID = fopen(OutputFile,'w');
for j = 1:size(CellData,1)
    fprintf(FID,TxtForm,CellData{j,:});
end
fclose(FID);
