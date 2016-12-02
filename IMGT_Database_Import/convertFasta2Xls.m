%convertFasta2Xls takes a fasta file and converts it xlsx file.
%
%  XlsName = Fasta2Xls()  Open a dialogue box for FASTA file finding
%  XlsName = Fasta2Xls(S) Where S is the Fasta file or structure variable.

function XlsName = convertFasta2Xls(varargin)
%Opening FASTA file.
if isempty(varargin); %If there's no input, look for FASTA file and open.
    [filename, pathname] = uigetfile('*.*','Open FASTA File.');
    S = fastaread([pathname filename]);
elseif length(varargin) == 1 %Determine what S is.
    if ischar(varargin{1}) %This is the original FASTA file. Convert to struct first.
        S = fastaread(varargin{1});
        filename = varargin{1};
    elseif isstruct(varargin{1}) 
        S = varargin{1};
        filename = '.xlsx';
    else
        error('File must be FASTA file or stucture variable from fastaread.');
    end
else
    error('Too many input variables');
end
C = struct2cell(S)';%Convert S into cell table.

%Set the savename correctly;
finddot = find(filename == '.');
if isempty(finddot) == 0
    suggestname = filename(1:finddot(end)-1);
else
    suggestname = filename;
end
[savename, pathname] = uiputfile('*.xlsx','Save Excel File.',suggestname);
if (length(savename) < 5) || (strcmp(savename(end-4:end),'.xlsx') == 0)
    savename = [savename '.xlsx'];
end

%Save the file to Excel file
xlswrite([pathname savename],C);
XlsName = savename;