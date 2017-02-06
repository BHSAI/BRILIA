%convertDataFormat will take VDJdata and output the Adaptive Biotech data
%format. This is because most data come from Adaptive, and our statistics
%programs works with these formats. The final format actually contains both
%VDJdata info + Adaptive format info. Look at Headers.xlsx file to see the
%relationship.

[FileName, FilePath] = uigetfile('*.xlsx;*.csv','Find VDJdata format file');
[~, ~, FileExt] = fileparts([FilePath FileName]);
if strcmpi(FileExt,'.csv');
    VDJdata = readDlmFile([FilePath FileName],'delimiter',';');
elseif strcmpi(FileExt,'.xlsx')
    [~, ~, VDJdata] = xlsread([FilePath FileName]);
else
    error('Unrecognized file extension. csv or xlsx needed');
end
[VDJdata, Header] = filterHeader(VDJdata);

HeaderFile = 'Headers.xlsx';
[~, ~, HeaderData] = xlsread(HeaderFile);
[HeaderData, TempName] = filterHeader(HeaderData,1);
AdapHeaderLoc = findHeader(TempName,'Adap');
Header2 = HeaderData(:,AdapHeaderLoc)';
FullData = cell(size(VDJdata,1), size(Header2,2));

%First fill in the obvious direct matches
SetLoc = findHeader(Header2,Header);
GetLoc = 1:length(Header);
DelZero = SetLoc == 0;
SetLoc(DelZero) = [];
GetLoc(DelZero) = [];
FullData(:,SetLoc) = VDJdata(:,GetLoc);

%Now carefully process the translation from VDJdata to AdapData

for j = 1:size(VDJdata,1)
    %1) Fill out the index location of the VMDNJ
    Map1 = findHeader(Header,{'vAlignLength' 'n2AlignLength' 'dAlignLength' 'n1AlignLength' 'jAlignLength'});
    Map2 = findHeader(Header2,{'vIndex' 'n2Index' 'dIndex' 'n1Index' 'jIndex'});
    TempVal = VDJdata(j,Map1);
    if ischar(TempVal{1});
        for p = 1:length(TempVal)
            TempVal{p} = str2double(TempVal{p});
        end
    end
    VMDNJ = cell2mat(TempVal);
    FullData(j,Map2) = num2cell([1 cumsum(VMDNJ(1:end-1))+1]);
    
    Map3 = findHeader(Header2,{'n2Insertion','n1Insertion'});
    FullData(j,Map3) = num2cell(VMDNJ(1,[2 4]));
    
    clear Map1 Map2 Map3

    %2) Fill out the gene family names
    Map1 = findHeader(Header,{'vMaxResolved' 'dMaxResolved' 'jMaxResolved'}); %([10; 16; 22]); %Gene name location
    Map2(1,1:4) = findHeader(Header2,{'vMaxResolved' 'vFamilyName' 'vGeneName' 'vGeneAllele'}); %[10:13; 17:20; 24:27]; %MaxRes GeneFam GeneNum GeneAllele location
    Map2(2,1:4) = findHeader(Header2,{'dMaxResolved' 'dFamilyName' 'dGeneName' 'dGeneAllele'}); %[10:13; 17:20; 24:27]; %MaxRes GeneFam GeneNum GeneAllele location
    Map2(3,1:4) = findHeader(Header2,{'jMaxResolved' 'jFamilyName' 'jGeneName' 'jGeneAllele'}); %[10:13; 17:20; 24:27]; %MaxRes GeneFam GeneNum GeneAllele location
    for k = 1:length(Map1)
        ParsedName = parseGeneName(VDJdata{j,Map1(k)});
        %Determine the max resolution data.
        if size(ParsedName,1) > 1
            NewParsedName = cell(1,4);
            %Look for common gene portion
            StopRes = 1;
            for r = 1:3 %counter marks where there is no single resolution
                UnqName = unique(ParsedName(:,r+1));
                if length(UnqName) == 1
                    NewParsedName(1,r+1) = UnqName;
                else
                    StopRes = r;
                    FullName = UnqName{1};
                    for d = 2:length(UnqName)
                        FullName = [FullName '| ' UnqName{d}];
                    end
                    NewParsedName{1,r+1} = FullName;
                end
            end
            if StopRes <= 2;
                StopRes = 3;
            end
            NewParsedName(1,1) = NewParsedName(1,StopRes-1);
            FullData(j,Map2(k,:)) = NewParsedName;
        else
            FullData(j,Map2(k,:)) = ParsedName;
        end
    end
    clear Map1 Map2
    
    %3) Fill out the gene substitution count
    Map1 = findHeader(Header,{'vMatchCount' 'dMatchCount' 'jMatchCount'}); %# of NTs that matched the alignment
    Map2 = findHeader(Header2,{'vAlignSubstitutionCount' 'dAlignSubstitutionCount' 'jAlignSubstitutionCount'}); %# of Nts that MISmatched that alignment
    q = 1;
    for r = 1:3
        FullData{j,Map2(r)} = VMDNJ(q) - VDJdata{j,Map1(r)};
        q = q+2;
    end
    clear Map1 Map2
end

%Select output file name
DotLoc = find(FileName == '.');
DotLoc = DotLoc(end);
OutputFilePre = FileName(1:DotLoc-1);
xlswrite([OutputFilePre '.Adap.xlsx'],cat(1,Header2,FullData))
