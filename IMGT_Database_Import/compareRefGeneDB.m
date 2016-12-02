%compareRefGeneDB will compare two sets of reference sequence to determine
%if they are the same, and also to identify which genes are missing from
%which database.

%  [SameNames, MissingFrom1, MissingFrom2] = compareRefGeneDB;
%    SameNames are the common names found
%    MissingFrom1 are the missing names from DB1 not common
%    MissingFrom2 are the missing names from DB2 not common


function varargout = compareRefGeneDB

[FileName1, FilePath1] = uigetfile('*.csv;*.xlsx;*.fa;*.txt;*.fasta','Select 1st database');
[FileName2, FilePath2] = uigetfile('*.csv;*.xlsx;*.fa;*.txt;*.fasta','Select 2nd database');

[~,~,FileExt1] = parseFileName([FilePath1 FileName1]);
[~,~,FileExt2] = parseFileName([FilePath2 FileName2]);

if strcmpi(FileExt1,'.csv')
    Data1 = readDlmFile([FilePath1 FileName1]);
elseif strcmpi(FileExt1,'.xlsx')
    [~,~,Data1] = xlsread([FilePath1 FileName1]);
else
    Data1= struct2cell(fastaread([FilePath1 FileName1]))';
end

if strcmpi(FileExt2,'.csv')
    Data2 = readDlmFile([FilePath2 FileName2]);
elseif strcmpi(FileExt1,'.xlsx')
    [~,~,Data2] = xlsread([FilePath2 FileName2]);
else
    Data2= struct2cell(fastaread([FilePath2 FileName2]))';
end

%This should be the raw files from IMGT, with the headers with the |
%seperators

%Extracting all the gene names
Name1 = cell(size(Data1,1),1);
for j = 1:size(Name1,1)
    Header = Data1{j,1};
    IGHLoc = regexp(Header,'IGH[VDJ]','once');
    SepLoc = regexp(Header,'\|');
    SepLoc = SepLoc(SepLoc>IGHLoc);
    SepLoc = SepLoc(1);
    GeneName = Header(IGHLoc:SepLoc-1);
    Name1{j} = strrep(GeneName,' ','');   
end

Name2 = cell(size(Data2,1),1);
for j = 1:size(Name2,1)
    Header = Data2{j,1};
    IGHLoc = regexp(Header,'IGH[VDJ]','once');
    SepLoc = regexp(Header,'\|');
    SepLoc = SepLoc(SepLoc>IGHLoc);
    SepLoc = SepLoc(1);
    GeneName = Header(IGHLoc:SepLoc-1);
    Name2{j} = strrep(GeneName,' ','');   
end

%Determine the intersect of the two names
[SameName, Idx1, Idx2] = intersect(Name1,Name2);

varargout{1} = SameName;
Name1(Idx1) = [];
varargout{2} = Name1;
Name2(Idx2) = [];
varargout{3} = Name2;

DotLoc = find(FileName1 == '.');
SuggestName = [FileName1(1:DotLoc(end)-1) 'Common'];
[SaveName, SavePath] = uiputfile('*.xlsx','Saving the common names, using DB1 format',SuggestName);
xlswrite([SavePath,SaveName],Data1(Idx1,:));

