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
addParameter(P, 'Species', '', @(x) ischar(x) && ismember(lower(x), getGeneDatabase('getlist'))); %Note, anything empty for DB filter will prompt user input.
addParameter(P, 'Strain', '', @ischar);
addParameter(P, 'Ddirection', '', @(x) ischar(x) && ismember(lower(x), {'all', 'fwd', 'inv', ''}));
addParameter(P, 'Vfunction', '', @(x) ischar(x) && min(ismember(regexpi(lower(x), ', ', 'split'), {'all', 'f', 'p', 'orf', ''}))==1);
addParameter(P, 'TDTon', '', @(x) ischar(x) && ismember(lower(x), {'y', 'n', ''}));
addParameter(P, 'Chain', '', @(x) ischar(x) && ismember(lower(x), {'h', 'l', 'hl'}));

[P, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {P, Pu, ExpPs, ExpPu};
   return
end

%Fill in the missing info by asking the user. Round values if needed.

%Set up the species at least.
P.Species = lower(P.Species);
SpeciesList = getGeneDatabase('getlist');
while isempty(P.Species) || ~ismember(P.Species, SpeciesList)
    dispList(SpeciesList);
    Selection = round(input('Choose the species.  '));
    if ~isempty(Selection) && Selection > 0 && Selection <= length(SpeciesList)
        P.Species = SpeciesList{Selection};
    end        
end

%Set how many unique VDJ germline seq to make
P.CloneCount = round(P.CloneCount);
while isempty(P.CloneCount) || (P.CloneCount <= 0) %Must check isempty first before checking odd value of P.CloneCount.
    P.CloneCount = round(input('How many productive VDJ seq, or B cell clone, do  you want? Default 10.  '));
    if isempty(P.CloneCount)
        P.CloneCount = 10;
    end
end

%Set up the clonally-related sequences per group, including germline VDJ
P.BranchLength = round(P.BranchLength);
while isempty(P.BranchLength) || (P.BranchLength < 0) %Must check isempty first before checking odd value of P.CloneCount.
    P.BranchLength = round(input('How many linear descendant per lineage branch? Default 5.  '));
    if isempty(P.BranchLength)
        P.BranchLength = 5;
    end
end

%Set up how many branches per group to make
P.BranchWidth = round(P.BranchWidth);
while isempty(P.BranchWidth) || (P.BranchWidth < 0) %Must check isempty first before checking odd value of P.CloneCount.
    P.BranchWidth = round(input('How many lineage branch to per clonal gorup? Default 2.  '));
    if isempty(P.BranchWidth)
        P.BranchWidth = 2;
    end
end

SeqCount = (P.BranchWidth * P.BranchLength + 1) * P.CloneCount; %Total number of sequences to generate

%Setup SHM rate
while isempty(P.SHMperc) || (P.SHMperc < 0 || P.SHMperc > 100)
    P.SHMperc = lower(input('What SHM mutation % to simulate? Default 2   '));
    if isempty(P.SHMperc)
        P.SHMperc = 2;
    end
end
if P.SHMperc > 20
    fprintf('Warning: P.SHMperc is set too high. This mean %0.1f%% of sequence mutates per descendant. \n', P.SHMperc);
end

%Setup TDT activity and N region max size
P.TDTon = lower(P.TDTon);
while isempty(P.TDTon) || (P.TDTon ~= 'n' && P.TDTon ~= 'y')
    P.TDTon = lower(input('Do you want TDT to be added? y or n. Default y   ', 's'));
    if isempty(P.TDTon)
        P.TDTon = 'y';
    end
end

%Setup the P.Chain (Might as well use both now)
P.Chain = lower(P.Chain);
while isempty(P.Chain) || ~ismember(P.Chain, {'h', 'l', 'hl'})
    P.Chain = lower(input('What IG chain do you want? H (default), L, HL', 's'));
    if isempty(P.Chain)
        P.Chain = 'h';
    end
end

%Setup the P.Vfunction
P.Vfunction = lower(P.Vfunction);
while isempty(P.Vfunction) || ~ismember(lower(P.Vfunction), {'all', 'f', 'p', 'orf'})
    P.Vfunction = lower(input('What V function? F (default), P, ORF, All ', 's'));
    if isempty(P.Vfunction)
        P.Vfunction = 'F';
    end
end

%Setup the P.Ddirection
P.Ddirection = lower(P.Ddirection);
while isempty(P.Ddirection) || ~ismember(lower(P.Ddirection), {'all', 'dfwd', 'drev'})
    P.Ddirection = lower(input('What D direction? All (default), Dfwd, Drev ', 's'));
    if isempty(P.Ddirection)
        P.Ddirection = 'all';
    end
end

%--------------------------------------------------------------------------
%Prepare the database for generating germline sequences

%Load databases and filter reference genese according to specifications
DB = getGeneDatabase(P.Species);
[DB, Pfilt] = filterGeneDatabase(DB, 'P.Strain', P.Strain, 'P.Ddirection', P.Ddirection, 'P.Vfunction', P.Vfunction);

%Update P based on Pfilt
PfiltFields = fieldnames(Pfilt);
for q = 1:length(PfiltFields)
    P.(PfiltFields{q}) = Pfilt.(PfiltFields{q});
end

%Delete invalid map entries\
M = getMapHeaderVar(DB.MapHeader);
FieldNames = {'Vmap' 'Dmap' 'Jmap' 'Vkmap' 'Jkmap' 'Vlmap' 'Jlmap'};
for f = 1:length(FieldNames)
    KeepIdx = ones(size(DB.(FieldNames{f}), 1), 1, 'logical');
    Xmap = DB.(FieldNames{f});
    for j = 1:size(Xmap, 1)
        if isempty(Xmap{j, M.Seq})
            KeepIdx(j) = 0;
        end
    end
    DB.(FieldNames{f}) = Xmap(KeepIdx, :);
end

%--------------------------------------------------------------------------
%Prepared all the parameters for simulations based on inputs

%P.Species-specfic deletion and N region lengths.
if strcmpi(P.Species, 'mouse')
    %These were collected from the C57BL6 mice data set from A Collins, 
    %2015, processed with old version of BRILIA, v1.9.0.
    MeanV3del = 1;
    MeanD5del = 4.5;
    MeanD3del = 3.4;
    MeanJ5del = 3.9;
    MeanMlen = 3.8;
    MeanNlen = 2.9;
else%if strcmpi(P.Species, 'human')
    %Souto-Carneiro, M.M., et al., Characterization of the Human Ig Heavy
    %P.Chain Antigen Binding Complementarity Determining Region 3 Using a
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
if P.TDTon == 'n'
    InsCap = 0;  %No tdt
else
    InsCap = 13; %Maximum N length
end

%--------------------------------------------------------------------------
%Begin creating and filling in simulated VDJdata entries

%Create the VDJdata default matrix
[VDJdata, VDJheader] = getBlankDataTable(SeqCount, P.Chain);
H = getHeavyHeaderVar(VDJheader);
L = getLightHeaderVar(VDJheader);
VDJdata(:, H.TemplateLoc) = num2cell(ones(SeqCount, 1)); %Always initialize TempCt column with 1.
VDJdata(:, H.SeqNumLoc) = num2cell(1:SeqCount); %Always assign a unique numbering order for VDJdata
VDJdata(:, H.GrpNumLoc) = num2cell(1:SeqCount); %Always assign a unique numbering order for VDJdata

SeqNum = 1;
GrpNum = 1;
while SeqNum <= SeqCount
    Tdata = VDJdata(SeqNum, :); %For temporary single entries;
    
    if ~isempty(strfind(P.Chain, 'h'))
        ValidH = 0;
        while ValidH == 0

            %Select a V, D, J gene
            Vnum = ceil(rand(1)*size(DB.Vmap, M.Seq));
            Dnum = ceil(rand(1)*size(DB.Dmap, M.Seq));
            Jnum = ceil(rand(1)*size(DB.Jmap, M.Seq));

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
            Cloc3 = DB.Vmap{Vnum, M.Anchor}-2; %Distance from 3' edge of Vgene to 3rd nt of 104 codon
            Wloc1 = DB.Jmap{Jnum, M.Anchor}; %Distance from 5' edge of Jgene to 1st nt of 118 codon
            if Cloc3 <= 0 || Wloc1 <= 0; continue; end %This will not allow pseudogenes without W/C locations to work.
            if V3del >= Cloc3; continue; end 
            if J5del >= Wloc1; continue; end
            if (D5del + D3del) >= length(DB.Dmap{Dnum, M.Seq}); continue; end %for simulation, prevent full deletion of D, though it could happen in reality.

            %Generate the sequence segments and store length information
            Vseq = DB.Vmap{Vnum, M.Seq}(1:end-V3del); %Vsegment
            Mseq = generateNregion(Mlen, ACGTprob, FLIPprobVD); %Nvd segment
            Dseq = DB.Dmap{Dnum, M.Seq}(D5del+1:end-D3del); %D segment
            Nseq = generateNregion(Nlen, ACGTprob, FLIPprobDJ); %Ndj segment
            Jseq = DB.Jmap{Jnum, M.Seq}(J5del+1:end); %Jsegment

            Seq = sprintf('%s%s%s%s%s', Vseq, Mseq, Dseq, Nseq, Jseq);
            VMDNJ = [length(Vseq) length(Mseq) length(Dseq) length(Nseq) length(Jseq)]; %All segment lengths

            %Check now if you'll get an in-frame CDR3
            Cloc = VMDNJ(1) + V3del - DB.Vmap{Vnum, M.Anchor} + 1;
            Wloc = sum(VMDNJ(1:4)) - J5del + DB.Jmap{Jnum, M.Anchor} + 2;
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
        Tdata(H.GeneNameLoc) = {DB.Vmap{Vnum, M.Gene} DB.Dmap{Dnum, M.Gene} DB.Jmap{Jnum, M.Gene}};
        Tdata(H.GeneNumLoc) = {Vnum Dnum Jnum}; %{DB.Vmap{Vnum, M.MiscLoc} DB.Dmap{Dnum, M.MiscLoc} DB.Jmap{Jnum, M.MiscLoc}};
        Tdata(H.CDR3Loc) = {CDR3seq length(CDR3seq) Cloc Wloc};
        Tdata(H.FunctLoc) = {'Y'};
    end
        
    if ~isempty(strfind(P.Chain, 'l'))
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
            Cloc3 = Vxmap{Vnum, M.Anchor} - 2; %Distance from 3' edge of Vgene to 3rd nt of 104 codon
            Wloc1 = Jxmap{Jnum, M.Anchor}; %Distance from 5' edge of Jgene to 1st nt of 118 codon
            if Cloc3 <= 0 || Wloc1 <= 0; continue; end %This will not allow pseudogenes without W/C locations to work.
            if V3del >= Cloc3; continue; end 
            if J5del >= Wloc1; continue; end

            %Generate the sequence segments and store length information
            Vseq = Vxmap{Vnum, M.Seq}(1:end-V3del); %Vsegment
            Nseq = generateNregion(Nlen, ACGTprob, FLIPprobDJ); %N segment
            Jseq = Jxmap{Jnum, M.Seq}(J5del+1:end); %Jsegment

            Seq = sprintf('%s%s%s', Vseq, Nseq, Jseq);
            VNJ = [length(Vseq) length(Nseq) length(Jseq)]; %All segment lengths

            %Check now if you'll get an in-frame CDR3
            Cloc = VNJ(1) + V3del - Vxmap{Vnum, M.Anchor} + 1;
            Wloc = sum(VNJ(1:2)) - J5del + Jxmap{Jnum, M.Anchor} + 2;
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
        Tdata(L.GeneNameLoc) = {Vxmap{Vnum, M.Gene} Jxmap{Jnum, M.Gene}};
        Tdata(L.GeneNumLoc) = {Vnum Jnum}; %{Vxmap{Vnum, M.MiscLoc}, Jxmap{Jnum, M.MiscLoc}};
        Tdata(L.CDR3Loc) = {CDR3seq length(CDR3seq) Cloc Wloc};
        Tdata(L.FunctLoc) = {'Y'};
    end
    
    %Perform SHM here.
    if P.BranchLength > 0 && P.BranchWidth > 0
        AllowRepeatMut = 'y'; %Allow SHM in same position as prior SHM?
        Tdata = generateSHMseq(Tdata, VDJheader, P.SHMperc, P.BranchLength, P.BranchWidth, AllowRepeatMut);
    end
    
    %Update to VDJdata
    VDJdata(SeqNum:SeqNum+size(Tdata, 1)-1, :) = Tdata;
    SeqNum = SeqNum + size(Tdata, 1);
    GrpNum = GrpNum + 1;
end

%Perform final updates and renumber sequences
Map = getVDJmapper(VDJheader);
VDJdata = updateVDJdata(VDJdata, Map, DB);
VDJdata = renumberVDJdata(VDJdata, VDJheader, 'seq', 'grp');

%Create seq names and update
SeqNum = cell2mat(VDJdata(:, H.SeqNumLoc));
GrpNum = cell2mat(VDJdata(:, H.SeqNumLoc));
MaxSeqDigit = floor(log10(max(SeqNum)))+1;
MaxGrpDigit = floor(log10(max(GrpNum)))+1;

if strcmpi(P.Chain, 'H')
    StrPat = ['Seq%0' num2str(MaxSeqDigit) 'd_Grp%0' num2str(MaxGrpDigit) 'd_%s_%s_%s_Len(%d|%d|%d|%d|%d)_Del(%d|%d|%d|%d)'];
    SeqNames = cell(size(SeqNum));
    for s = 1:length(SeqNames)
        SeqNames{s} = sprintf(StrPat, SeqNum(s), GrpNum(s), VDJdata{s, H.GeneNameLoc}, VDJdata{s, H.LengthLoc}, VDJdata{s, H.DelLoc});
    end
    VDJdata(:, H.SeqNameLoc) = SeqNames;
elseif strcmpi(P.Chain, 'L')
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
if isempty(P.SaveAs)
    SaveName = sprintf('SIM_%s_%s_%s_SHM%0.0f.csv', P.Species, P.P.Strain, P.Chain, P.SHMperc);
    [SaveName, SavePath] = uiputfile(SaveName, 'Save simulated sequences as');
    P.SaveAs = [SavePath SaveName];
end
[SavePath, SaveName, ~] = parseFileName(P.SaveAs);
DotLoc = find(SaveName == '.');
SaveNamePre = SaveName(1:DotLoc(end)-1);%Get the prefix

%Save the BRILIA annotation
VDJdata = findCDR1(VDJdata, Map, DB);
VDJdata = findCDR2(VDJdata, Map, DB);
VDJdata = findCDR3(VDJdata, Map, DB, 'imgt');
saveSeqData(P.SaveAs, VDJdata, VDJheader); 

%Save just the input file format only
if strcmpi(P.Chain, 'HL')
    SaveLoc = [H.SeqNameLoc H.SeqLoc L.SeqLoc H.TemplateLoc];
elseif strcmpi(P.Chain, 'H')
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
SaveName_Setting = [SavePath SaveNamePre '.txt'];
makeSettingFile(SaveName_Setting, P);

if nargout >= 1
    varargout{1} = VDJdata;
    if nargout >= 2
        varargout{2} = VDJheader;
    end
end        
