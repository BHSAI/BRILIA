%generateVDJseq will generate VDJ genes by randomly selecting germline VDJ
%genes, randomly deleting nts from gene edges, randomly creating N regions,
%and then combining them into a single sequence. Somatic hyper mutation
%simulations are performed if SHMrate > 0. A mutated sequence is added to
%the dataset only if it reaches SHMrate level of mutation. Only productive
%VDJ rearrangements are stored since these are the one that can go undergo
%SHM. We are not considering unproductive VDJ genes in the 2nd chromosome
%that could also undergo SHM if there is a productive gene.
%
%  [VDJdata, NewHeader] = generateVDJseq(Param1,Value1,...)
%
%  INPUT
%    Param           Value                   This sets:
%    ------------    ----------------------  -------------------------
%    SaveFileName    [FileName].csv          Place to store VDJdata. If
%                                              empty, asks users to choose.
%    CloneCount      1 <= N <= 10000         # of germline VDJ seq to make
%    BranchLength    1 <= N <= 20            # of linear descendants per
%                                              branch per group
%    BranchCount     1 <= N <= 20            # of branches per group
%    SHMrate         0 <= N <= 100           % of SHMs per descendant
%    TDTon           'y' 'n'                 Include TDT-based N regions
%    Species         'human' 'mouse' etc     VDJ database by species
%    Strain          'all' 'C57BL' etc       VDJ database by strain
%    Ddirection      'all' 'fwd' 'inv'       Allowed D gene direction
%    Vfunction       'all' 'f' 'p' 'orf'     Allowed V gene functions
%
%  OUTPUT
%    VDJdata: cell matrix storing all information about each VDJgene
%    NewHeader: cell matrix storing data field names for VDJdata
%  
%  NOTE
%    The total number of sequences generated will be CloneCount *
%    (AddSize+1) * BranchSize.
%
%    Only a linear phylogeny tree is created, starting with the germline
%    VDJ seq, followed by (CloneSize-1) number of descendant sequences.
%    
%    Simulation of SHMs does not account for SHM hot spots.
%
%    VDJ recombination event does not account for bias in certain gene
%    usage. Instead, selects gene uniformly.
%
%    VDJ gene edge deletions follow exponential distrbution with mean
%    deletion counts determined by literature data. Max deletion is capped
%    at 13 to prevent complete gene deletion issues.
%
%    Nregions are synthesized according to the behavior of Terminal
%    Deoxynucleotidyl Transferase. which mostly adds A's and G's, or C's
%    and T's.
%
%    Only genes with well-defined C and W locations can be used.
function [VDJdata, NewHeader] = generateVDJseq(varargin)
%--------------------------------------------------------------------------
%Parse inputs to this simulator
P = inputParser;
addParameter(P,'SaveFileName','',@ischar);
addParameter(P,'CloneCount',[],@isnumeric);
addParameter(P,'BranchLength',[],@isnumeric);
addParameter(P,'BranchCount',[],@isnumeric);
addParameter(P,'SHMrate',[],@(x) isnumeric(x) && (x>=0) && (x<=100)); %For clustering purposes. Set empty to force user to input it later.
addParameter(P,'Species','',@(x) ischar(x) && ismember(lower(x),{'human','mouse',''})); %Note, anything empty for DB filter will prompt user input.
addParameter(P,'Strain','',@ischar);
addParameter(P,'Ddirection','',@(x) ischar(x) && ismember(lower(x),{'all','fwd','inv',''}));
addParameter(P,'Vfunction','',@(x) ischar(x) && min(ismember(regexpi(lower(x),',','split'),{'all','f','p','orf',''}))==1);
addParameter(P,'TDTon','',@(x) ischar(x) && ismember(lower(x),{'y','n',''}));
parse(P,varargin{:});
P = P.Results;
SaveFileName = P.SaveFileName;
CloneCount = P.CloneCount;
BranchLength = P.BranchLength;
BranchCount = P.BranchCount;
SHMrate = P.SHMrate;
TDTon = P.TDTon;
Species = P.Species;
Strain = P.Strain;
Ddirection = P.Ddirection;
Vfunction = P.Vfunction;

%--------------------------------------------------------------------------
%Fill in the missing info by asking the user. Round values if needed.

%Setup the save file name
if isempty(SaveFileName)
    [SaveName, SavePath] = uiputfile('*.csv','Save simulated sequences as');
    SaveFileName = [SavePath SaveName];
end
DotLoc = find(SaveFileName == '.');
SaveFileName_Fasta = [SaveFileName(1:DotLoc(end)-1) '.fa'];
SaveFileName_Setting = [SaveFileName(1:DotLoc(end)-1) '.txt'];

%Set up the species at least.
Species = lower(Species);
SpeciesList = {'mouse','human'};
while isempty(Species) || ~ismember(Species,SpeciesList);
    dispList(SpeciesList);
    Selection = round(input('Choose the species.  '));
    if ~isempty(Selection) && Selection > 0 && Selection <= length(SpeciesList)
        Species = SpeciesList{Selection};
    end        
end

%Set how many unique VDJ germline seq to make
CloneCount = round(CloneCount);
while isempty(CloneCount) || (CloneCount <= 0) %Must check isempty first before checking odd value of CloneCount.
    CloneCount = round(input('How many productive VDJ seq, or B cell clone, do  you want? Default 1000.  '));
    if isempty(CloneCount)
        CloneCount = 1000;
    end
end

%Set up the clonally-related sequences per group, including germline VDJ
BranchLength = round(BranchLength);
while isempty(BranchLength) || (BranchLength <= 0) %Must check isempty first before checking odd value of CloneCount.
    BranchLength = round(input('How many linear descendant per lineage branch? Default 1.  '));
    if isempty(BranchLength)
        BranchLength = 1;
    end
end

%Set up how many branches per group to make
BranchCount = round(BranchCount);
while isempty(BranchCount) || (BranchCount <= 0) %Must check isempty first before checking odd value of CloneCount.
    BranchCount = round(input('How many lineage branch to per clonal gorup? Default 1.  '));
    if isempty(BranchCount)
        BranchCount = 1;
    end
end

SeqCount = (BranchCount * BranchLength + 1) * CloneCount; %Total number of sequences to generate

%Setup SHMrate
while isempty(SHMrate) || (SHMrate < 0 || SHMrate > 100)
    SHMrate = lower(input('What SHM mutation % to simulate? Default 5   '));
    if isempty(SHMrate)
        SHMrate = 5;
    end
end
if SHMrate > 20
    fprintf('Warning: SHMrate is set too high. This mean %0.1f%% of sequence mutates per descendant. \n',SHMrate);
end

%Setup TDT activity and N region max size
TDTon = lower(TDTon);
while isempty(TDTon) || (TDTon ~= 'n' && TDTon ~= 'y')
    TDTon = lower(input('Do you want TDT to be added? y or n. Default y   ','s'));
    if isempty(TDTon)
        TDTon = 'y';
    end
end

%--------------------------------------------------------------------------
%Prepare the database for generating germline sequences

%Load databases and filter reference genese according to specifications
[Vmap, Dmap, Jmap] = getCurrentDatabase('change',Species);
[Vmap, Dmap, Jmap, FiltOption] = filterRefGene(Vmap,Dmap,Jmap,'Strain',Strain,'Ddirection',Ddirection,'Vfunction',Vfunction,'KeepThis','yes');
Strain = FiltOption.Strain;
Vfunction = FiltOption.Vfunction;
Ddirection = FiltOption.Ddirection;

%Create the valid map idx 
Vidx = 1:size(Vmap,1);
for j = 1:length(Vidx)
    if isempty(Vmap{j,1}) || Vmap{j,10} <= 0
        Vidx(j) = 0;
    end
end
Vidx(Vidx==0) = [];

%Create the valid map idx 
Didx = 1:size(Dmap,1);
for j = 1:length(Didx)
    if isempty(Dmap{j,1})
        Didx(j) = 0;
    end
end
Didx(Didx==0) = [];

%Create the valid map idx 
Jidx = 1:size(Jmap,1);
for j = 1:length(Jidx)
    if isempty(Jmap{j,1}) || Jmap{j,10} <= 0
        Jidx(j) = 0;
    end
end
Jidx(Jidx==0) = [];

%--------------------------------------------------------------------------
%Prepared all the parameters for simulations based on inputs

%Species-specfic deletion and N region lengths.
if strcmpi(Species,'mouse')
    %These were collected from the C57BL6 mice data set from A Collins,
    %2015, processed with old version of BRILIA, v1.9.0.
    MeanV3del = 1;
    MeanD5del = 4.5;
    MeanD3del = 3.4;
    MeanJ5del = 3.9;
    MeanMlen = 3.8;
    MeanNlen = 2.9;
elseif strcmpi(Species,'human')
    %Souto-Carneiro, M.M., et al., Characterization of the Human Ig Heavy
    %Chain Antigen Binding Complementarity Determining Region 3 Using a
    %Newly Developed Software Algorithm, JOINSOLVER. The Journal of
    %Immunology, 2004. 172(11): p. 6790-6802.
    MeanV3del = 2.1;
    MeanD5del = 4.5;
    MeanD3del = 5.0;
    MeanJ5del = 6.8;
    MeanMlen = 7;
    MeanNlen = 7;
end
DelCap = 13; %Maximum deletion set for all species, for now.

%TDT insertion and flipping bias (based on mouse data)
ACGTprob = [0.25 0.08 0.60 0.07]; %Prob of A, C, G, T, respectively.
FLIPprobVD = 0.25; %Probability of flipping Nvd side.
FLIPprobDJ = 0.45; %Probability of flipping Ndj side.
if TDTon == 'n'
    InsCap = 0;  %No tdt
else
    InsCap = 13; %Maximum N length
end

%--------------------------------------------------------------------------
%Begin creating and filling in simulated VDJdata entries

%Create the VDJdata default matrix
HeaderData = readDlmFile('Headers_BRILIA.csv','Delimiter',';'); %Obtain the VDJdata header info for output format
NewHeader = HeaderData(2:end,1)';
getHeaderVar;
VDJdata = cell(SeqCount,length(NewHeader));
VDJdata(:,TemplateLoc) = num2cell(ones(SeqCount,1)); %Always initialize TempCt column with 1.
VDJdata(:,SeqNumLoc) = num2cell(1:SeqCount); %Always assign a unique numbering order for VDJdata
VDJdata(:,GrpNumLoc) = num2cell(1:SeqCount); %Always assign a unique numbering order for VDJdata

SeqNum = 1;
GrpNum = 1;
while SeqNum <= SeqCount
    %Select a V, D, J gene
    Vsel = ceil(rand(1)*length(Vidx));
    Dsel = ceil(rand(1)*length(Didx));
    Jsel = ceil(rand(1)*length(Jidx));
    Vnum = Vidx(Vsel);
    Dnum = Didx(Dsel);
    Jnum = Jidx(Jsel);
    
    %Skip pseudo and orf. Only functionals are count. (might filter before)
    if ~strcmpi(Vmap{Vnum,7},'f'); continue; end

    %Determine NT deletions and N1/N2 lengths
    V3del = round(-MeanV3del*log(rand(1)));
    D5del = round(-MeanD5del*log(rand(1)));
    D3del = round(-MeanD3del*log(rand(1)));
    J5del = round(-MeanJ5del*log(rand(1)));
    Mlen = round(-MeanMlen*log(rand(1)));
    Nlen = round(-MeanNlen*log(rand(1)));
    
    %Make sure to cap the deletion and insertion lengths
    if V3del > DelCap; V3del = DelCap;  end
    if D5del > DelCap; D5del = DelCap;  end
    if D3del > DelCap; D3del = DelCap;  end
    if J5del > DelCap; J5del = DelCap;  end
    if Mlen  > InsCap; Mlen  = InsCap;  end
    if Nlen  > InsCap; Nlen  = InsCap;  end

    %Make sure deletions avoid 104C and 118W codons, and whole D gene
    Cloc3 = Vmap{Vnum,10}-2; %Distance from 3' edge of Vgene to 3rd nt of 104 codon
    Wloc1 = Jmap{Jnum,10}; %Distance from 5' edge of Jgene to 1st nt of 118 codon
    if Cloc3 <= 0 || Wloc1 <= 0; continue; end %This will not allow pseudogenes without W/C locations to work.
    if V3del >= Cloc3; continue; end 
    if J5del >= Wloc1; continue; end
    if (D5del + D3del) >= length(Dmap{Dnum,1}); continue; end %for simulation, prevent full deletion of D, though it could happen in reality.

    %Generate the sequence segments and store length information
    Vseq = Vmap{Vnum,1}(1:end-V3del); %Vsegment
    Mseq = generateNregion(Mlen,ACGTprob,FLIPprobVD); %Nvd segment
    Dseq = Dmap{Dnum,1}(D5del+1:end-D3del); %D segment
    Nseq = generateNregion(Nlen,ACGTprob,FLIPprobDJ); %Ndj segment
    Jseq = Jmap{Jnum,1}(J5del+1:end); %Jsegment
    
    Seq = sprintf('%s%s%s%s%s',Vseq,Mseq,Dseq,Nseq,Jseq);
    VMDNJ = [length(Vseq) length(Mseq) length(Dseq) length(Nseq) length(Jseq)]; %All segment lengths
    
    %Check now if you'll get an in-frame CDR3
    Cloc = VMDNJ(1) + V3del - Vmap{Vnum,10} + 1;
    Wloc = sum(VMDNJ(1:4)) - J5del + Jmap{Jnum,10} + 2;
    if mod(Wloc-Cloc+1,3) ~= 0; continue; end %Not in-frame
    CDR3seq = nt2aa(Seq(Cloc:Wloc),'acgtonly','false');
    CDR3len = length(CDR3seq);
    if ~isempty(regexpi(CDR3seq,'\*','once')); continue; end %Stop codon   
    
    %Save information temporarily to Tdata
    Tdata = cell(1,length(NewHeader)); %For temporary single entries;
    Tdata([SeqNumLoc GrpNumLoc SeqLoc RefSeqLoc SeqNameLoc]) = {SeqNum, GrpNum, Seq Seq num2str(SeqNum)}; %Note that Seq and RefSeq are same.
    Tdata(LengthLoc) = num2cell(VMDNJ);
    Tdata(DelLoc) = {V3del D5del D3del J5del};
    Tdata(FamLoc) = {Vmap{Vnum,2} Dmap{Dnum,2} Jmap{Jnum,2}};
    Tdata(FamNumLoc) = {Vnum, Dnum, Jnum};
    Tdata(CDR3Loc) = {CDR3seq CDR3len Cloc Wloc};
    Tdata(FunctLoc) = {'Y'};

    %Perform SHM here.
    if BranchLength > 1
        AllowRepeatMut = 'y'; %Allow SHM in same position as prior SHM?
        Tdata = generateSHMseq(Tdata,NewHeader,SHMrate,BranchLength,BranchCount,AllowRepeatMut);
    end
    
    %Update to VDJdata
    VDJdata(SeqNum:SeqNum+size(Tdata,1)-1,:) = Tdata;
    SeqNum = SeqNum + size(Tdata,1);
    GrpNum = GrpNum + 1;
end 

%Perform final updates and renumber sequences
VDJdata = updateVDJdata(VDJdata,NewHeader,Vmap,Dmap,Jmap);
VDJdata = renumberVDJdata(VDJdata,NewHeader,'seq','grp');

%Create seq names and update
SeqNum = cell2mat(VDJdata(:,SeqNumLoc));
GrpNum = cell2mat(VDJdata(:,SeqNumLoc));
MaxSeqDigit = floor(log10(max(SeqNum)))+1;
MaxGrpDigit = floor(log10(max(GrpNum)))+1;
StrPat = ['Seq%0' num2str(MaxSeqDigit) 'd_Grp%0' num2str(MaxGrpDigit) 'd_%s_%s_%s_Len(%d|%d|%d|%d|%d)_Del(%d|%d|%d|%d)'];
SeqNames = cell(size(SeqNum));
for s = 1:length(SeqNames)
    SeqNames{s} = sprintf(StrPat,SeqNum(s),GrpNum(s),VDJdata{s,FamLoc},VDJdata{s,LengthLoc},VDJdata{s,DelLoc});
end
VDJdata(:,SeqNameLoc) = SeqNames;

%Save the simulated VDJdata
saveSeqData(SaveFileName,VDJdata,NewHeader,'Delimiter',';');

%Save the simulate sequence only as fasta
fastawrite(SaveFileName_Fasta,VDJdata(:,SeqNameLoc),VDJdata(:,SeqLoc));

%Save the settings used to simulate this, but convert matrix to str.
P.CloneCount = CloneCount;
P.BranchLength = BranchLength;
P.BranchCount = BranchCount;
P.SHMrate = SHMrate;
P.TDTon = TDTon;
P.Species = Species;
P.Strain = Strain;
P.Ddirection = Ddirection;
P.Vfunction = Vfunction;
P.MeanV3del = MeanV3del;
P.MeanD5del = MeanD5del;
P.MeanD3del = MeanD3del;
P.MeanJ5del = MeanJ5del;
P.MeanMlen = MeanMlen;
P.MeanNlen = MeanNlen;
P.InsCap = InsCap;
P.DelCap = DelCap;
P.TDT_ACGTprob = ACGTprob;
P.TDT_FLIPprobVD = FLIPprobVD;
P.TDT_FLIPprobDJ = FLIPprobDJ;
P.SaveFileName = SaveFileName;
P.SeqCount = SeqCount;

makeSettingFile(SaveFileName_Setting,P);