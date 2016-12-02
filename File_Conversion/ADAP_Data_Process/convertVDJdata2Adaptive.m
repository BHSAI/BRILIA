%convertDataFormat will take VDJdata and output the Adaptive Biotech data
%format. This is because most data come from Adaptive, and our statistics
%programs works with these formats. The final format actually contains both
%convertVDJdata2Adaptive will convert the VDJdata format to Adaptive
%format.

function convertVDJdata2Adaptive(varargin)
KeepHeader = 'Y'; %Default, keep header in file.
if ~isempty(varargin)
    for h = 1:length(varargin)
        if ischar(varargin{h})
            if strcmpi(varargin{h},'NoHeader')
                KeepHeader = 'N';
                varargin(h) = [];
                break
            end
        end
    end
end


%Get the header information
HeaderFile = 'Headers_Adaptive.xlsx';
[~, ~, HeaderData] = xlsread(HeaderFile);
Header2 = HeaderData(2:end,1)';

%Assign similar locations name to VDJdata format
SeqLoc2 = findHeader(Header2,'nucleotide');
RefSeqLoc2 = findHeader(Header2,'RefSeq');
SeqNumLoc2 = findHeader(Header2,'SeqNum');
GrpNumLoc2 = findHeader(Header2,'GroupNum');
TemplateLoc2 = findHeader(Header2,'count (templates)'); 

LengthLoc2 = findHeader(Header2,{'vAlignLength','n2Insertion','dAlignLength','n1Insertion','jAlignLength'});
DelLoc2 = findHeader(Header2,{'vDeletion','d5Deletion','d3Deletion','jDeletion'});
SHMLoc2 = findHeader(Header2,{'vAlignSubstitutionCount','dAlignSubstitutionCount','jAlignSubstitutionCount'});

FunctLoc2 = findHeader(Header2,'sequenceStatus');
CDR3Loc2 = findHeader(Header2,{'aminoAcid' 'cdr3Length'});

VnameLoc = findHeader(Header2,{'vMaxResolved','vFamilyName','vGeneName','vGeneAllele'});
DnameLoc = findHeader(Header2,{'dMaxResolved','dFamilyName','dGeneName','dGeneAllele'});
JnameLoc = findHeader(Header2,{'jMaxResolved','jFamilyName','jGeneName','jGeneAllele'});

%Select the file names and begin filling in the data
if isempty(varargin)
    [FileNames, FilePath] = uigetfile('*.xlsx;*.csv','Open the VDJdata format file','multiselect','on');
else
    FileNames = varargin{1};
    if length(varargin) == 2
        FilePath = varargin{2};
    else
        if ispc
            FilePath = [cd '\'];
        else
            FilePath = [cd '/'];
        end
    end
end
if ischar(FileNames)
    FileNames = {FileNames};
end

%Do this for all files
for f = 1:length(FileNames)
    [VDJdata,NewHeader,FileName,FilePath] = openSeqData([FilePath FileNames{f}]);
    getHeaderVar;
    
    SaveData = cell(size(VDJdata,1),length(Header2));
    
    %Transfering over direct matches
    VDJdataLoc = [SeqLoc  SeqNumLoc  GrpNumLoc  LengthLoc  DelLoc  TemplateLoc  SHMLoc   FunctLoc  CDR3Loc(1:2) RefSeqLoc ];
    SaveDataLoc =[SeqLoc2 SeqNumLoc2 GrpNumLoc2 LengthLoc2 DelLoc2 TemplateLoc2 SHMLoc2  FunctLoc2 CDR3Loc2     RefSeqLoc2];   
    DelThese = VDJdataLoc == 0 | SaveDataLoc == 0; %Make sure there is a corresponding column. Otherwise remove.
    VDJdataLoc(DelThese) = [];
    SaveDataLoc(DelThese) = [];
    
    SaveData(:,SaveDataLoc) = VDJdata(:,VDJdataLoc);
    SaveData(:,CDR3Loc2(2)) = num2cell(cell2mat(SaveData(:,CDR3Loc2(2)))*3); %VDJdata save the amino acid count, wherease adap saves nt count in CDR3.
    
    %Getting VDJdata names
    [Vmap,Dmap,Jmap] = getCurrentDatabase;
    Xmap = {Vmap Dmap Jmap};
    for j = 1:size(VDJdata,1)
        for k = 1:3
            Xnum = VDJdata{j,FamNumLoc(k)};
            if ~isempty(Xnum) && min(~isnan(Xnum) == 1)
                Xname = Xmap{k}(Xnum(1),3:6);
            else
                Xname = repmat({'unresolved'},1,4);
            end
            
            switch k
                case 1
                    SaveData(j,VnameLoc) = Xname;
                case 2
                    SaveData(j,DnameLoc) = Xname;
                case 3
                    SaveData(j,JnameLoc) = Xname;
            end
        end
    end 
    
    %Select output file name
    DotLoc = find(FileName == '.');
    DotLoc = DotLoc(end);
    OutputFilePre = [FileName(1:DotLoc-1) '.Adap'];
    
    %Save to excel or csv file, depending on OS
    if strcmpi(KeepHeader,'Y')
        SaveData = cat(1,Header2,SaveData);
    end
%     if ispc
%         xlswrite([FilePath OutputFilePre '.xlsx'],SaveData);
%     else
        writeDlmFile(SaveData,[FilePath OutputFilePre '.csv'],'\t');
%     end    
end