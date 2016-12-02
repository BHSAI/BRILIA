%writeDlmFile will take a cell array and convert it into a string array,
%and then save it into a delimited file. Works better than dlmwrite and
%csvwrite.
%
%  writeDlmFile(CellData,OutputFile,D)
%    CellData = cell array (without an embeded cell in a cell);
%    OutputFile = full name of the file to be saved to
%
%    
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
D = '\t';
if length(varargin) == 1
    D = varargin{1};
end
DataFormat = ['%s' D];
TxtForm = [repmat(DataFormat,1,size(CellData,2)-1) '%s\n'];
FID = fopen(OutputFile,'w');
for j = 1:size(CellData,1)
    fprintf(FID,TxtForm,CellData{j,:});
end
fclose(FID);