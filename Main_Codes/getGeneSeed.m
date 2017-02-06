%getGeneSeed will look through the database V or J genes, and extract N
%number of nts of the framework regions, including the conserved C and W
%codon set. 
%
%  SeedSeq = getGeneSeed(Xmap,X,Nleft,Nright,Alphabet)
%
%  INPUT
%    Xmap: either Vmap or Dmap from getCurrentDatabase code
%    X ['V' or 'J']: Specifiy if it's a V or J gene
%    Nleft: number of nucleotides or amino acids left of anchor
%    Nright: number of nucleotides or amino acids right of anchor
%    Alphabet ['nt' or 'aa']: Choose NT or AA letters.
%
%  OUTPUT
%    SeedSeq: cell array of NT or AA sequences that can be used for
%    alignment purposes. 
%
%  NOTE
%    Genes without a solid anchor point C or W are REMOVED from SeedSeq
%    output! This is to prevent confusion pinpointing where the seed
%    sequence marks the position of the C or W.

function SeedSeq = getGeneSeed(Xmap,X,Nleft,Nright,Alphabet)
KeyNTloc = 10; %Column in Xmap that store the position of the 104C 1st nt from the 5' end, or 118W 3rd nt from the 3' end.
FuncLoc = 7;
H.SeqLoc = 1; 

%Remove all sequences without conserved residue or empty sequences
DelThese = zeros(size(Xmap,1),1,'logical');
for j = 1:size(Xmap,1)
    if Xmap{j,KeyNTloc} < 3
        DelThese(j) = 1;
    elseif isempty(Xmap{j,H.SeqLoc})
        DelThese(j) = 1;
    elseif ~strcmpi(Xmap{j,FuncLoc},'F')
        DelThese(j) = 1;
    end   
end
Xmap(DelThese,:) = [];

%Determine how many N number of nts to keep
if strcmpi(Alphabet,'aa')
    Nleft = Nleft*3;
    Nright = Nright*3+2; % The 2 comes from the anchor codon 2nd-3rd nt.
end

%Trim all sequence to contain the conserved residue
SeedSeq = cell(size(Xmap,1),1);
if strcmpi(X,'V')
    for j = 1:size(Xmap,1)
        SeedSeq{j,1} = padtrimSeq(Xmap{j,H.SeqLoc}, length(Xmap{j,H.SeqLoc})-Xmap{j,KeyNTloc}+1, Nleft, Nright);
    end
elseif strcmpi(X,'J')
    for j = 1:size(Xmap,1)
        SeedSeq{j,1} = padtrimSeq(Xmap{j,H.SeqLoc}, Xmap{j,KeyNTloc}, Nleft, Nright);
    end
end
%Translate set to AA if needed
if strcmpi(Alphabet,'aa')
    SeedSeq = nt2aa(SeedSeq,'ACGTOnly','false');
end

%Get only unique seeds
SeedSeq = unique(SeedSeq);
