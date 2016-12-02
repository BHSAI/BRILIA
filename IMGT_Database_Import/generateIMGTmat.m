%generateIMGTmat converts a 3-sheet xlsx files containing the
%fasta-2-xlsx-converted IMGT databases into a IMGT_Mouse_VDJgene.mat file
%that is used as a lookup map. The xlsx file must be formatted as such:
%   1) Go to http://www.imgt.org/vquest/refseqh.html and download V D J
%      fasta files.
%   2) Run the convertFasta2Xls.m function to convert the fasta files into
%      xlsx files. 
%   3) Run the parseField.m function to parse the fasta file headers into 
%      new columns. 
%   4) Combine the V D J xlsxs files into one file called
%      "IMGT_Species_VDJ.xlsx" that have 3 sheets: Vgene, Dgene, Jgene 
%   5) Run this function genereateIMGTmat.m to create the lookup matrix for 
%      NT-AA-GeneName.
%
%      DefaultFile = generateIMGTmat()  Looks for the 
%      IMGT_Mouse_VDJgene.xlsx file, otherwise asks the user to look for
%      it. Returns the DefaultFile name, saved the the current directory.
%
%      DefaultFile = 'IMGT_Mouse_VDJgene.mat';


function GeneFile = generateIMGTmat(varargin)
%Look for the Ref Seq IMGT xlsx file with the Vgene, Dgene, Jgene sheets.
DefaultPath = mfilename('fullpath');
SlashLoc = regexp(DefaultPath,'\\|\/');
FilePath = DefaultPath(1:SlashLoc(end));
[FileName1,FilePath] = uigetfile('*.xlsx','Find the IMGT gene database file');
[FileName2,FilePath] = uigetfile('*.xlsx','Find the IMGT gene table file');
[SaveName,SavePath] = uiputfile('*.xlsx','Save the files as');
if strcmpi(SaveName(end-4:end),'.xlsx') == 0
    SaveName = [SaveName '.xlsx'];
end

%Select species
SpeciesList = {'Mouse','Human'};
for k = 1:length(SpeciesList)
    disp([num2str(k) ') ' SpeciesList{k}]);
end
while 1
    SpeciesNum = input('What species is this for? ');
    if SpeciesNum >=1 && SpeciesNum <= length(SpeciesList)
        break
    end
end
Species = lower(SpeciesList{SpeciesNum});

%Read in the informations
[~, ~, VgeneDB] = xlsread([FilePath, FileName1],'Vgene');
[~, ~, DgeneDB] = xlsread([FilePath, FileName1],'Dgene');
[~, ~, JgeneDB] = xlsread([FilePath, FileName1],'Jgene');
[~, ~, Vtable] = xlsread([FilePath, FileName2],'Vgene');
[~, ~, Dtable] = xlsread([FilePath, FileName2],'Dgene');
[~, ~, Jtable] = xlsread([FilePath, FileName2],'Jgene');
Vtable = removeNAN(Vtable);
Dtable = removeNAN(Dtable);
Jtable = removeNAN(Jtable);
TableHeader = Vtable(1,:);
StrainLoc = findHeader(TableHeader,'Strain','all');
NameLoc = findHeader(TableHeader,'IMGT full name');

%Setup the NT, AA, Name map.
Header = {'nucleotide' 'GeneNameIMGT' 'STDgeneName' 'geneSubgroup' 'geneName' 'geneAllele' 'function' 'readingFrame' 'isolatedHost' 'keyNTloc'};
ColNum = length(Header);

Vmap = cell(size(VgeneDB,1),ColNum);
Dmap = cell(size(DgeneDB,1)*2,ColNum); %Need to double this for fwd/rv direction
Jmap = cell(size(JgeneDB,1),ColNum); 

for j = 1:size(Vmap,1)
    TempNT = upper(VgeneDB{j,2});
    TempNT(TempNT == '.') = []; %Delete the "..." annotation in IMGT ref seq.
    Vmap{j,1} = TempNT;    
    Vmap{j,2} = VgeneDB{j,4}; %IMGT default gene name format   
    Vmap(j,3:6) = parseGeneName(VgeneDB{j,4},Species); %Renamed for easy sorting formating    
    Vmap{j,7} = VgeneDB{j,6}; %Functionality
    if strcmpi(Species,'Mouse')
        Vmap{j,9} = findStrainMatch(Vmap{j,2},Vtable,StrainLoc,NameLoc);
    end
    RF = VgeneDB{j,10};
    if ischar(RF); RF = 0; end
    Vmap{j,8} = RF; %Reading Frame (0) if there' is none
%     [Cloc,~,~] = findVgeneC(Vmap(j,1));
%     Vmap{j,10} = Cloc;
end
Vmap = findVgeneC(Vmap);
Vmap = sortrows(Vmap,4);

%Create D gene map, Dmap. NOTE: Dmap can go in reverse, so double this.
k = 1;
for j = 1:2:size(DgeneDB,1)*2
    Dmap{j,1} = upper(DgeneDB{k,2});
    Dmap{j,2} = DgeneDB{k,4};
    Dmap(j,3:6) = parseGeneName(DgeneDB{k,4},Species);  
    Dmap{j,7} = DgeneDB{k,6};   
    Dmap{j,8} = 1; %D RF doesn't matter
    if strcmpi(Species,'Mouse')
        Dmap{j,9} = findStrainMatch(Dmap{j,2},Dtable,StrainLoc,NameLoc);
    end
    Dmap{j,10} = 0;
    
    Dmap{j+1,1} = seqrcomplement(Dmap{j,1});
    TempName = Dmap{j,2};
    DashLoc = strfind(TempName,'-');
    if isempty(DashLoc)
        DashLoc = strfind(TempName,'*');
    end
%    Dmap{j+1,2} = [TempName(1:DashLoc-1) 'r' TempName(DashLoc:end)];    
    Dmap{j+1,2} = ['r' TempName];    
    Dmap(j+1,3:6) = parseGeneName(Dmap{j+1,2},Species);  
    Dmap{j+1,7} = Dmap{j,7};   
    Dmap{j+1,8} = 1; %D RF doesn't matter
    Dmap{j+1,9} = Dmap{j,9};
    Dmap{j+1,10} = 0;
    
    k = k+1;
end
Dfwd = sortrows(Dmap(1:2:end,:),4);
Drev = sortrows(Dmap(2:2:end,:),4);
Dmap(1:2:end,:) = Dfwd;
Dmap(2:2:end,:) = Drev;


for j = 1:size(Jmap,1)
    TempNT = upper(JgeneDB{j,2});
    TempNT(TempNT == '.') = []; %Delete the "..." annotation in IMGT ref seq.
    Jmap{j,1} = TempNT;    
    Jmap{j,2} = JgeneDB{j,4}; %IMGT default gene name format   
    Jmap(j,3:6) = parseGeneName(JgeneDB{j,4}); %Renamed for easy sorting formating    
    Jmap{j,7} = JgeneDB{j,6}; %Functionality
    if strcmpi(Species,'Mouse')
        Jmap{j,9} = findStrainMatch(Jmap{j,2},Jtable,StrainLoc,NameLoc);
    end
    
    RF = JgeneDB{j,10};
    if ischar(RF); RF = 0; end
    Jmap{j,8} = RF; %Reading Frame (0) if there' is none
    [Wloc,~,~] = findJgeneWF(Jmap(j,1));
    Jmap{j,10} = Wloc;
end
Jmap = sortrows(Jmap,4);


%Go through all seq, replace "N" with "X" because that's what convolveSeq
%can intake for variable nts.
for k = 1:size(Vmap,1)
    Seq = Vmap{j,1};
    Nloc = regexpi(Seq,'N');
    if ~isempty(Nloc)
        Seq(Nloc) = 'X';
        Vmap{j,1} = Seq;
        k
    end
end

for k = 1:size(Dmap,1)
    Seq = Dmap{j,1};
    Nloc = regexpi(Seq,'N');
    if ~isempty(Nloc)
        Seq(Nloc) = 'X';
        Dmap{j,1} = Seq;
        k
    end
end

for k = 1:size(Jmap,1)
    Seq = Jmap{j,1};
    Nloc = regexpi(Seq,'N');
    if ~isempty(Nloc)
        Seq(Nloc) = 'X';
        Jmap{j,1} = Seq;
        k
    end
end


%Save as the default file name in current directory
DotLoc = regexp(SaveName,'\.');
GeneFile = [SaveName(1:DotLoc(end)-1) '.mat'];
save([SavePath GeneFile],'Vmap','Dmap','Jmap','Header');

if ispc
    xlswrite([SavePath SaveName],Vmap,'Vmap');
    xlswrite([SavePath SaveName],Dmap,'Dmap');
    xlswrite([SavePath SaveName],Jmap,'Jmap');
    xlswrite([SavePath SaveName],Header,'Header');
end
