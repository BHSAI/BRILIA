%fixDegenVDJ will remove pseudogene (P) and open reading frame (ORF)
%annotations if a functional gene annotations are also suggested by BRILIA.
%
%  VDJdata = fixPseudoVDJ(VDJdata, VDJheader, DB)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    DB: Gene database structure (getCurrentDatabase.m)
%
%  OUTPUT
%    VDJdata: modified VDJdata where pseudogenes and ORF suggestion are
%      removed if amongst function gene suggestions.

function VDJdata = fixDegenVDJ(VDJdata, Map, DB)
%Determine chain and extract key locations
if strcmpi(Map.Chain, 'HL')
    VDJ = {'V', 'D', 'J', 'Vk', 'Jk', 'Vl', 'Jl'};
elseif strcmpi(Map.Chain, 'H')
    VDJ = {'V', 'D', 'J'};
elseif strcmpi(Map.Chain, 'L')
    VDJ = {'Vk', 'Jk', 'Vl', 'Jl'};
end

%Determine what DB maps are available
M = getMapHeaderVar(DB.MapHeader);
MapNames = fieldnames(DB);
ValidMap = findCell(MapNames, 'map', 'MatchWord', 'partial');
MapNames = MapNames(ValidMap(ValidMap > 0));

%Remove Pseudogenes and ORF annotations if amongst a functional gene
UpdateIdx = zeros(size(VDJdata, 1), 1, 'logical');
for x = 1:length(VDJ)
    %Determine location of gene name and map number
    switch VDJ{x}
        case 'V'
            GeneNameLoc = Map.hGeneName(1);
            GeneNumLoc  = Map.hGeneNum(1);
        case 'D'
            GeneNameLoc = Map.hGeneName(2);
            GeneNumLoc  = Map.hGeneNum(2);
        case 'J'
            GeneNameLoc = Map.hGeneName(3);
            GeneNumLoc  = Map.hGeneNum(3);
        case 'Vk'
            GeneNameLoc = Map.lGeneName(1);
            GeneNumLoc  = Map.lGeneNum(1);
        case 'Jk'
            GeneNameLoc = Map.lGeneName(2);
            GeneNumLoc  = Map.lGeneNum(2);
        case 'Vl'
            GeneNameLoc = Map.lGeneName(1);
            GeneNumLoc  = Map.lGeneNum(1);
        case 'Jl'
            GeneNameLoc = Map.lGeneName(2);
            GeneNumLoc  = Map.lGeneNum(2);
    end
    
    %Make sure the database exists and is not empty
    DBname = sprintf('%smap', VDJ{x});
    DBidx = findCell(MapNames, DBname);
    if DBidx(1) == 0; continue; end %DB doesn't exist, do not try
    Xmap = DB.(DBname);
    
    %Determine the Locus of the database
    if length(VDJ{x}) == 1
        Locus = 'H';
    else
        Locus = upper(VDJ{x}(2));
    end

    %Identify location of pseudo genes
    PseudoIdx = zeros(size(Xmap, 1), 1, 'logical');
    for v = 1:size(Xmap, 1)
        if ~isempty(regexpi(Xmap{v, M.FunctLoc}, 'P|ORF'))
            PseudoIdx(v) = 1;
        end
    end

    %Iteratively find groups
    GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
    UnqGrpNum = unique(GrpNum);
    for y = 1:length(UnqGrpNum)
        %Find the sequences per grp, and the Xnum
        IdxLoc = find(UnqGrpNum(y) == GrpNum);   

        Xnum = VDJdata{IdxLoc(1), GeneNumLoc};
        Xname = VDJdata{IdxLoc(1), GeneNameLoc};
        if length(Xnum) <= 1; continue; end
        if isempty(Xname); continue; end
        
        %Determine this sequence's locus. Locus must match.
        if Locus ~= 'H'
            LocusLoc = regexpi(Xname, 'IG[HKL]', 'once');
            SeqLocus = Xname(LocusLoc+2);
            if SeqLocus ~= Locus %Locus mismatch, so skip
                continue
            end
        end

        %Look for any pseudogene
        DelThis = zeros(length(Xnum), 1, 'logical');
        for k = 1:length(Xnum)
            if PseudoIdx(Xnum) == 1
                DelThis(k) = 1;
            end
        end
        HowManyPseudo = sum(DelThis);

        %Delete pseudo or ORF ONLY IF it is not the only one
        if HowManyPseudo >= 1 && (HowManyPseudo < length(Xnum)) %Must have at least 1 instance, but not all instances
            %Need to adjust the names
            Xnum(DelThis) = [];
            RepPat = repmat('%s|', 1, length(Xnum));
            RepPat(end) = [];
            Xname = sprintf(RepPat, Xmap{:, Xnum});

            %Update VDJdata family number and name
            VDJdata(IdxLoc, [GeneNumLoc GeneNameLoc]) = repmat({Xnum Xname}, length(IdxLoc), 1);
            UpdateIdx(IdxLoc) = 1;
        end
    end
end

%Update those that have changed
if max(UpdateIdx) > 0
    VDJdata(UpdateIdx, :) = buildRefSeq(VDJdata(UpdateIdx, :), Map, DB, Map.Chain, 'germline', 'first'); %must do first seq of all cluster
    VDJdata(UpdateIdx, :) = updateVDJdata(VDJdata(UpdateIdx, :), Map, DB);
    fprintf('%s: Corrected %d pseudogene entries \n', mfilename, sum(UpdateIdx));
end
