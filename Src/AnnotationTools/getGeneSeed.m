%getGeneSeed will look through the database V or J genes, and extract N
%number of nts of the framework regions, including the conserved C and W
%codon set. 
%
%  SeedSeq = getGeneSeed(DB,X,Nleft,Nright,Alphabet)
%
%  INPUT
%    DB: gene database structure(getCurrentDatabase.m)
%    X ['V' 'Vk' 'Vl' 'J' 'Jk' 'Jl']: Specificy what database gene is used.
%    Nleft: number of NT or AA left of anchor
%    Nright: number of NT or AA right of anchor
%    Alphabet ['nt' or 'aa']: Choose NT or AA letters
%
%  OUTPUT
%    SeedSeq: cell array of NT or AA sequences that can be used for
%      alignment purposes. 
%
%  NOTE
%    Genes without a solid anchor point C or W are REMOVED from SeedSeq
%    output! This is to prevent confusion pinpointing where the seed
%    sequence marks the position of the C or W.

function SeedSeq = getGeneSeed(DB,X,Nleft,Nright,Alphabet)
%Determine which database to use
M = getMapHeaderVar(DB.MapHeader);
Fields = fieldnames(DB);
DBidx = findCell(Fields,[X 'map']);
Xmap = DB.(Fields{DBidx});

%Make sure the map isn't empty
if isempty(Xmap)
    SeedSeq = {''};
    return;
end

%Remove all sequences without conserved residue or empty sequences
DelThese = zeros(size(Xmap,1),1,'logical');
for j = 1:size(Xmap,1)
    if Xmap{j,M.AnchorLoc} < 3
        DelThese(j) = 1;
    elseif isempty(Xmap{j,M.SeqLoc})
        DelThese(j) = 1;
    elseif ~strcmpi(Xmap{j,M.FunctLoc},'F')
        DelThese(j) = 1;
    end   
end
Xmap(DelThese,:) = [];

%Make sure the map isn't empty
if isempty(Xmap)
    SeedSeq = {''};
    return;
end

%Determine if this is a V or J
X = '';
for j = 1:size(Xmap,1)
    GeneName = Xmap{j,M.GeneLoc};
    StrPat = 'IG[HKL][VDJ]';
    Xidx = regexpi(GeneName,StrPat,'end');
    if ~isempty(Xidx)
        X = GeneName(Xidx);
        break
    end
end
if isempty(X)
    SeedSeq = {''};
    error('getGeneSeed: Unable to determine gene V or J');
    return;
end

%Determine how many N number of nts to keep
if strcmpi(Alphabet,'aa')
    Nleft = Nleft*3;
    Nright = Nright*3+2; % The 2 comes from the anchor codon 2nd-3rd nt.
end

%Trim all sequence to contain the conserved residue
SeedSeq = cell(size(Xmap,1),1);
if strcmpi(X,'V')
    for j = 1:size(Xmap,1)
        SeedSeq{j,1} = padtrimSeq(Xmap{j,M.SeqLoc}, length(Xmap{j,M.SeqLoc})-Xmap{j,M.AnchorLoc}+1, Nleft, Nright);
    end
elseif strcmpi(X,'J')
    for j = 1:size(Xmap,1)
        SeedSeq{j,1} = padtrimSeq(Xmap{j,M.SeqLoc}, Xmap{j,M.AnchorLoc}, Nleft, Nright);
    end
end

%Translate set to AA if needed
if strcmpi(Alphabet,'aa')
    SeedSeq = nt2aa(SeedSeq,'ACGTOnly','false');
end

%Get only unique seeds
SeedSeq = unique(SeedSeq);
