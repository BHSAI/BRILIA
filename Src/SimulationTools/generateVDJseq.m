%generateVDJseq will generate VDJ genes by randomly selecting germline VDJ
%genes, randomly deleting nts from gene edges, randomly creating N regions, 
%and then combining them into a single sequence. Somatic hyper mutation
%simulations are performed if SHMrate > 0. A mutated sequence is added to
%the dataset only if it reaches SHMrate level of mutation. Only productive
%VDJ rearrangements are stored since these are the one that can go undergo
%SHM. We are not considering unproductive VDJ genes in the 2nd chromosome
%that could also undergo SHM if there is a productive gene.
%
%  [VDJdata, VDJheader] = generateVDJseq(Param, Value, ...)
%
%  INPUT
%    Param           Value                 This sets:
%    ------------    --------------------  -------------------------
%    SaveAs          [FileName].csv        Place to store VDJdata. If
%                                            empty, asks users to choose.
%    CloneCount      1 <= N <= 10000       # of germline VDJ seq to make
%    BranchLength    1 <= N <= 20          # of linear descendants per
%                                            branch per group
%    BranchWidth     1 <= N <= 20          # of branches per group
%    SHMperc         0 <= N <= 100         % of SHMs per descendant
%    TDTon           'y' 'n'               Include TDT-based N regions
%    Species         'human' 'mouse' etc   VDJ database by species
%    Strain          'all' 'C57BL' etc     VDJ database by strain
%    Ddirection      'all' 'fwd' 'inv'     Allowed D gene direction
%    Vfunction       'all' 'f' 'p' 'orf'   Allowed V gene functions
%
%  OUTPUT
%    VDJdata: cell matrix storing all information about each VDJgene
%    VDJheader: cell matrix storing data field names for VDJdata
%  
%  NOTE
%    The total number of sequences generated will be CloneCount *
%    BranchLength * BranchCount.
%
%    Simulation of SHMs does not account for SHM hot spots yet.
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
function varargout = generateVDJseq(varargin)
P = inputParser;
addParameter(P, 'SaveAs', '', @ischar);
addParameter(P, 'CloneCount', [], @isnumeric);
addParameter(P, 'BranchLength', [], @isnumeric);
addParameter(P, 'BranchWidth', [], @isnumeric);
addParameter(P, 'SHMperc', [], @(x) isnumeric(x) && (x>=0) && (x<=100)); %For clustering purposes. Set empty to force user to input it later.
addParameter(P, 'Species', '', @(x) ischar(x) && ismember(lower(x), {'human', 'mouse', ''})); %Note, anything empty for DB filter will prompt user input.
addParameter(P, 'Strain', '', @ischar);
addParameter(P, 'Ddirection', '', @(x) ischar(x) && ismember(lower(x), {'all', 'fwd', 'inv', ''}));
addParameter(P, 'Vfunction', '', @(x) ischar(x) && min(ismember(regexpi(lower(x), ', ', 'split'), {'all', 'f', 'p', 'orf', ''}))==1);
addParameter(P, 'TDTon', '', @(x) ischar(x) && ismember(lower(x), {'y', 'n', ''}));
addParameter(P, 'Chain', '', @(x) ischar(x) && ismember(lower(x), {'h', 'l', 'hl'}));

[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end
SaveAs = Ps.SaveAs;
CloneCount = Ps.CloneCount;
BranchLength = Ps.BranchLength;
BranchWidth = Ps.BranchWidth;
SHMperc = Ps.SHMperc;
Species = Ps.Species;
Strain = Ps.Strain;
Ddirection = Ps.Ddirection;
Vfunction = Ps.Vfunction;
TDTon = Ps.TDTon;
Chain = Ps.Chain;

%Fill in the missing info by asking the user. Round values if needed.

%Set up the species at least.
Species = lower(Species);
SpeciesList = getGeneDatabase('getlist');
while isempty(Species) || ~ismember(Species, SpeciesList);
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
while isempty(BranchLength) || (BranchLength < 0) %Must check isempty first before checking odd value of CloneCount.
    BranchLength = round(input('How many linear descendant per lineage branch? Default 0.  '));
    if isempty(BranchLength)
        BranchLength = 0;
    end
end

%Set up how many branches per group to make
BranchWidth = round(BranchWidth);
while isempty(BranchWidth) || (BranchWidth < 0) %Must check isempty first before checking odd value of CloneCount.
    BranchWidth = round(input('How many lineage branch to per clonal gorup? Default 0.  '));
    if isempty(BranchWidth)
        BranchWidth = 0;
    end
end

SeqCount = (BranchWidth * BranchLength + 1) * CloneCount; %Total number of sequences to generate

%Setup SHM rate
while isempty(SHMperc) || (SHMperc < 0 || SHMperc > 100)
    SHMperc = lower(input('What SHM mutation % to simulate? Default 2   '));
    if isempty(SHMperc)
        SHMperc = 2;
    end
end
if SHMperc > 20
    fprintf('Warning: SHMperc is set too high. This mean %0.1f%% of sequence mutates per descendant. \n', SHMperc);
end

%Setup TDT activity and N region max size
TDTon = lower(TDTon);
while isempty(TDTon) || (TDTon ~= 'n' && TDTon ~= 'y')
    TDTon = lower(input('Do you want TDT to be added? y or n. Default y   ', 's'));
    if isempty(TDTon)
        TDTon = 'y';
    end
end

%Setup the Chain (Might as well use both now)
Chain = lower(Chain);
while isempty(Chain) || ~ismember(Chain, {'h', 'l', 'hl'})
    Chain = lower(input('What IG chain do you want? HL, H, L ', 's'));
    if isempty(Chain)
        Chain = 'hl';
    end
end

%Setup the Vfunction
Vfunction = lower(Vfunction);
while isempty(Vfunction) || ~ismember(lower(Vfunction), {'all', 'f', 'p', 'orf'})
    Vfunction = lower(input('What V function? F, P, ORF, All ', 's'));
    if isempty(Vfunction)
        Vfunction = 'F';
    end
end

%Setup the Ddirection
Ddirection = lower(Ddirection);
while isempty(Ddirection) || ~ismember(lower(Ddirection), {'all', 'dfwd', 'drev'})
    Ddirection = lower(input('What D direction? All, Dfwd, Drev ', 's'));
    if isempty(Ddirection)
        Ddirection = 'all';
    end
end

%--------------------------------------------------------------------------
%Prepare the database for generating germline sequences

%Load databases and filter reference genese according to specifications
DB = getGeneDatabase(Species);
[DB, Pfilt] = filterGeneDatabase(DB, 'Strain', Strain, 'Ddirection', Ddirection, 'Vfunction', Vfunction);

%Update P based on Pfilt
PfiltFields = fieldnames(Pfilt);
for q = 1:length(PfiltFields)
    Ps.(PfiltFields{q}) = Pfilt.(PfiltFields{q});
end

%Delete invalid map entries\
M = getMapHeaderVar(DB.MapHeader);
FieldNames = {'Vmap' 'Dmap' 'Jmap' 'Vkmap' 'Jkmap' 'Vlmap' 'Jlmap'};
for f = 1:length(FieldNames)
    KeepIdx = ones(size(DB.(FieldNames{f}), 1), 1, 'logical');
    Xmap = DB.(FieldNames{f});
    for j = 1:size(Xmap, 1)
        if isempty(Xmap{j, M.SeqLoc})
            KeepIdx(j) = 0;
        end
    end
    DB.(FieldNames{f}) = Xmap(KeepIdx, :);
end

%--------------------------------------------------------------------------
%Prepared all the parameters for simulations based on inputs

%Species-specfic deletion and N region lengths.
if strcmpi(Species, 'mouse')
    %These were collected from the C57BL6 mice data set from A Collins, 
    %2015, processed with old version of BRILIA, v1.9.0.
    MeanV3del = 1;
    MeanD5del = 4.5;
    MeanD3del = 3.4;
    MeanJ5del = 3.9;
    MeanMlen = 3.8;
    MeanNlen = 2.9;
else%if strcmpi(Species, 'human')
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
[VDJdata, VDJheader] = getBlankDataTable(SeqCount, Chain);
H = getHeavyHeaderVar(VDJheader);
L = getLightHeaderVar(VDJheader);
VDJdata(:, H.TemplateLoc) = num2cell(ones(SeqCount, 1)); %Always initialize TempCt column with 1.
VDJdata(:, H.SeqNumLoc) = num2cell(1:SeqCount); %Always assign a unique numbering order for VDJdata
VDJdata(:, H.GrpNumLoc) = num2cell(1:SeqCount); %Always assign a unique numbering order for VDJdata

SeqNum = 1;
GrpNum = 1;
while SeqNum <= SeqCount
    Tdata = VDJdata(SeqNum, :); %For temporary single entries;
    
    if ~isempty(strfind(Chain, 'h'))
        ValidH = 0;
        while ValidH == 0

            %Select a V, D, J gene
            Vnum = ceil(rand(1)*size(DB.Vmap, M.SeqLoc));
            Dnum = ceil(rand(1)*size(DB.Dmap, M.SeqLoc));
            Jnum = ceil(rand(1)*size(DB.Jmap, M.SeqLoc));

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
            Cloc3 = DB.Vmap{Vnum, M.AnchorLoc}-2; %Distance from 3' edge of Vgene to 3rd nt of 104 codon
            Wloc1 = DB.Jmap{Jnum, M.AnchorLoc}; %Distance from 5' edge of Jgene to 1st nt of 118 codon
            if Cloc3 <= 0 || Wloc1 <= 0; continue; end %This will not allow pseudogenes without W/C locations to work.
            if V3del >= Cloc3; continue; end 
            if J5del >= Wloc1; continue; end
            if (D5del + D3del) >= length(DB.Dmap{Dnum, M.SeqLoc}); continue; end %for simulation, prevent full deletion of D, though it could happen in reality.

            %Generate the sequence segments and store length information
            Vseq = DB.Vmap{Vnum, M.SeqLoc}(1:end-V3del); %Vsegment
            Mseq = generateNregion(Mlen, ACGTprob, FLIPprobVD); %Nvd segment
            Dseq = DB.Dmap{Dnum, M.SeqLoc}(D5del+1:end-D3del); %D segment
            Nseq = generateNregion(Nlen, ACGTprob, FLIPprobDJ); %Ndj segment
            Jseq = DB.Jmap{Jnum, M.SeqLoc}(J5del+1:end); %Jsegment

            Seq = sprintf('%s%s%s%s%s', Vseq, Mseq, Dseq, Nseq, Jseq);
            VMDNJ = [length(Vseq) length(Mseq) length(Dseq) length(Nseq) length(Jseq)]; %All segment lengths

            %Check now if you'll get an in-frame CDR3
            Cloc = VMDNJ(1) + V3del - DB.Vmap{Vnum, M.AnchorLoc} + 1;
            Wloc = sum(VMDNJ(1:4)) - J5del + DB.Jmap{Jnum, M.AnchorLoc} + 2;
            if mod(Wloc-Cloc+1, 3) ~= 0; continue; end %Not in-frame
            if ((Wloc - Cloc + 1) / 3) < 5; continue; end %Too short CDR3
            ReadFrame = mod(Cloc + 2, 3) + 1;
            AASeq = convNT2AA(Seq, 'frame', ReadFrame);
            if ~isempty(regexpi(AASeq, '\*', 'once')); continue; end %Stop codon   
            CDR3seq = convNT2AA(Seq(Cloc:Wloc));
            ValidH = 1;
        end

        %Save information temporarily to Tdata
        Tdata([H.SeqNumLoc H.GrpNumLoc H.SeqLoc H.RefSeqLoc H.SeqNameLoc]) = {SeqNum, GrpNum, Seq Seq num2str(SeqNum)}; %Note that Seq and RefSeq are same.
        Tdata(H.LengthLoc) = num2cell(VMDNJ);
        Tdata(H.DelLoc) = {V3del D5del D3del J5del};
        Tdata(H.GeneNameLoc) = {DB.Vmap{Vnum, M.GeneLoc} DB.Dmap{Dnum, M.GeneLoc} DB.Jmap{Jnum, M.GeneLoc}};
        Tdata(H.GeneNumLoc) = {Vnum Dnum Jnum}; %{DB.Vmap{Vnum, M.MiscLoc} DB.Dmap{Dnum, M.MiscLoc} DB.Jmap{Jnum, M.MiscLoc}};
        Tdata(H.CDR3Loc) = {CDR3seq length(CDR3seq) Cloc Wloc};
        Tdata(H.FunctLoc) = {'Y'};
    end
        
    if ~isempty(strfind(Chain, 'l'))
        ValidL = 0;
        while ValidL == 0

            %Select a locus and then V and J pair
            Loc = ceil(2*rand(1));
            if Loc == 1 %kappa
                Vxmap = DB.Vkmap;
                Jxmap = DB.Jkmap;
            else %lambda
                Vxmap = DB.Vlmap;
                Jxmap = DB.Jlmap;
            end
            Vnum = ceil(rand(1)*size(Vxmap, 1));
            Jnum = ceil(rand(1)*size(Jxmap, 1));

            %Determine NT deletions and N1/N2 lengths
            V3del = round(-MeanV3del*log(rand(1)));
            J5del = round(-MeanJ5del*log(rand(1)));
            Nlen = round(-MeanNlen*log(rand(1)));

            %Make sure to cap the deletion and insertion lengths
            if V3del > DelCap; V3del = DelCap;  end
            if J5del > DelCap; J5del = DelCap;  end
            if Nlen  > InsCap; Nlen  = InsCap;  end

            %Make sure deletions avoid 104C and 118W codons, and whole D gene
            Cloc3 = Vxmap{Vnum, M.AnchorLoc} - 2; %Distance from 3' edge of Vgene to 3rd nt of 104 codon
            Wloc1 = Jxmap{Jnum, M.AnchorLoc}; %Distance from 5' edge of Jgene to 1st nt of 118 codon
            if Cloc3 <= 0 || Wloc1 <= 0; continue; end %This will not allow pseudogenes without W/C locations to work.
            if V3del >= Cloc3; continue; end 
            if J5del >= Wloc1; continue; end

            %Generate the sequence segments and store length information
            Vseq = Vxmap{Vnum, M.SeqLoc}(1:end-V3del); %Vsegment
            Nseq = generateNregion(Nlen, ACGTprob, FLIPprobDJ); %N segment
            Jseq = Jxmap{Jnum, M.SeqLoc}(J5del+1:end); %Jsegment

            Seq = sprintf('%s%s%s', Vseq, Nseq, Jseq);
            VNJ = [length(Vseq) length(Nseq) length(Jseq)]; %All segment lengths

            %Check now if you'll get an in-frame CDR3
            Cloc = VNJ(1) + V3del - Vxmap{Vnum, M.AnchorLoc} + 1;
            Wloc = sum(VNJ(1:2)) - J5del + Jxmap{Jnum, M.AnchorLoc} + 2;
            if mod(Wloc-Cloc+1, 3) ~= 0; continue; end %Not in-frame
            if ((Wloc - Cloc + 1) / 3) < 5; continue; end %Too short
            ReadFrame = mod(Cloc + 2, 3) + 1;
            AASeq = convNT2AA(Seq, 'frame', ReadFrame);
            if ~isempty(regexpi(AASeq, '\*', 'once')); continue; end %Stop codon   
            CDR3seq = convNT2AA(Seq(Cloc:Wloc));
            
            ValidL = 1;
        end

        %Save information temporarily to Tdata
        Tdata([L.SeqNumLoc L.GrpNumLoc L.SeqLoc L.RefSeqLoc L.SeqNameLoc]) = {SeqNum, GrpNum, Seq, Seq, num2str(SeqNum)}; %Note that Seq and RefSeq are same.
        Tdata(L.LengthLoc) = num2cell(VNJ);
        Tdata(L.DelLoc) = {V3del J5del};
        Tdata(L.GeneNameLoc) = {Vxmap{Vnum, M.GeneLoc} Jxmap{Jnum, M.GeneLoc}};
        Tdata(L.GeneNumLoc) = {Vnum Jnum}; %{Vxmap{Vnum, M.MiscLoc}, Jxmap{Jnum, M.MiscLoc}};
        Tdata(L.CDR3Loc) = {CDR3seq length(CDR3seq) Cloc Wloc};
        Tdata(L.FunctLoc) = {'Y'};
    end
    
    %Perform SHM here.
    if BranchLength > 0 && BranchWidth > 0
        AllowRepeatMut = 'y'; %Allow SHM in same position as prior SHM?
        Tdata = generateSHMseq(Tdata, VDJheader, SHMperc, BranchLength, BranchWidth, AllowRepeatMut);
    end
    
    %Update to VDJdata
    VDJdata(SeqNum:SeqNum+size(Tdata, 1)-1, :) = Tdata;
    SeqNum = SeqNum + size(Tdata, 1);
    GrpNum = GrpNum + 1;
end

%Perform final updates and renumber sequences
VDJdata = updateVDJdata(VDJdata, VDJheader, DB);
VDJdata = renumberVDJdata(VDJdata, VDJheader, 'seq', 'grp');

%Create seq names and update
SeqNum = cell2mat(VDJdata(:, H.SeqNumLoc));
GrpNum = cell2mat(VDJdata(:, H.SeqNumLoc));
MaxSeqDigit = floor(log10(max(SeqNum)))+1;
MaxGrpDigit = floor(log10(max(GrpNum)))+1;

if strcmpi(Chain, 'H')
    StrPat = ['Seq%0' num2str(MaxSeqDigit) 'd_Grp%0' num2str(MaxGrpDigit) 'd_%s_%s_%s_Len(%d|%d|%d|%d|%d)_Del(%d|%d|%d|%d)'];
    SeqNames = cell(size(SeqNum));
    for s = 1:length(SeqNames)
        SeqNames{s} = sprintf(StrPat, SeqNum(s), GrpNum(s), VDJdata{s, H.GeneNameLoc}, VDJdata{s, H.LengthLoc}, VDJdata{s, H.DelLoc});
    end
    VDJdata(:, H.SeqNameLoc) = SeqNames;
elseif strcmpi(Chain, 'L')
    StrPat = ['Seq%0' num2str(MaxSeqDigit) 'd_Grp%0' num2str(MaxGrpDigit) 'd_%s_%s_Len(%d|%d|%d)_Del(%d|%d)'];
    SeqNames = cell(size(SeqNum));
    for s = 1:length(SeqNames)
        SeqNames{s} = sprintf(StrPat, SeqNum(s), GrpNum(s), VDJdata{s, L.GeneNameLoc}, VDJdata{s, L.LengthLoc}, VDJdata{s, L.DelLoc});
    end
    VDJdata(:, L.SeqNameLoc) = SeqNames;
else
    StrPat = ['Seq%0' num2str(MaxSeqDigit) 'd_Grp%0' num2str(MaxGrpDigit) 'd_%s_%s_%s_%s_%s'];
    SeqNames = cell(size(SeqNum));
    for s = 1:length(SeqNames)
        SeqNames{s} = sprintf(StrPat, SeqNum(s), GrpNum(s), VDJdata{s, H.GeneNameLoc}, VDJdata{s, L.GeneNameLoc});
    end
    VDJdata(:, H.SeqNameLoc) = SeqNames;
end

%--------------------------------------------------------------------------
%Setup the save file name
if isempty(SaveAs)
    SaveName = sprintf('SIM_%s_%s_%s_SHM%0.0f.csv', Species, Ps.Strain, Chain, SHMperc);
    [SaveName, SavePath] = uiputfile(SaveName, 'Save simulated sequences as');
    SaveAs = [SavePath SaveName];
end
[SavePath, SaveName, ~] = parseFileName(SaveAs);
DotLoc = find(SaveName == '.');
SaveNamePre = SaveName(1:DotLoc(end)-1);%Get the prefix

%Save the BRILIA annotation
VDJdata = findCDR1(VDJdata, VDJheader, DB);
VDJdata = findCDR2(VDJdata, VDJheader, DB);
VDJdata = findCDR3(VDJdata, VDJheader, DB, 'imgt');
saveSeqData(SaveAs, VDJdata, VDJheader); 

%Save just the input file format only
if strcmpi(Chain, 'HL')
    SaveLoc = [H.SeqNameLoc H.SeqLoc L.SeqLoc H.TemplateLoc];
elseif strcmpi(Chain, 'H')
    SaveLoc = [H.SeqNameLoc H.SeqLoc H.TemplateLoc];
else
    SaveLoc = [L.SeqNameLoc L.SeqLoc L.TemplateLoc];
end
SaveName_CSV = [SavePath SaveNamePre 'Input.csv'];
saveSeqData(SaveName_CSV, VDJdata(:, SaveLoc), VDJheader(SaveLoc)); %BRILIA annotation

% %Save just the fasta file format only
% if H.SeqLoc > 0
%     SaveName_FastaH = [SavePath SaveNamePre 'H.fa'];
%     fastawrite(SaveName_FastaH, VDJdata(:, H.SeqNameLoc), VDJdata(:, H.SeqLoc));
% end
% if L.SeqLoc > 0
%     SaveName_FastaL = [SavePath SaveNamePre 'L.fa'];
%     fastawrite(SaveName_FastaL, VDJdata(:, L.SeqNameLoc), VDJdata(:, L.SeqLoc));
% end

%Save the setting file, include more P parameters
Ps.MeanV3del = MeanV3del;
Ps.MeanD5del = MeanD5del;
Ps.MeanD3del = MeanD3del;
Ps.MeanJ5del = MeanJ5del;
Ps.MeanMlen = MeanMlen;
Ps.MeanNlen = MeanNlen;
Ps.InsCap = InsCap;
Ps.DelCap = DelCap;
Ps.TDT_ACGTprob = ACGTprob;
Ps.TDT_FLIPprobVD = FLIPprobVD;
Ps.TDT_FLIPprobDJ = FLIPprobDJ;
SaveName_Setting = [SavePath SaveNamePre '.txt'];
makeSettingFile(SaveName_Setting, Ps);

if nargout >= 1
    varargout{1} = VDJdata;
    if nargout >= 2
        varargout{2} = VDJheader;
    end
end        
