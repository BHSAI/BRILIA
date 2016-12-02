%findPnucleotides will check to see if there are p nucleotides based on the
%alignment and deletion info. 
%
%  Plen = findPnucleotides(Xalign,Xdel,Xmax,Xdir) 
%    Xalign is the VDJdata alignment results, 3xN string array
%    X0 is the position in the alignment where ref(bottom) seq is checked.
%    Xmax is the max # of NTs allowed in p-nucleotide match
%    Xdir is the direction check. 5 or 3 for 5' left side or 3' right side.
%    Xdel is # of NTs deleted. Mainly should be 0, to help bypass search.

function Plen = findPnucleotides(Xalign,X0,Xmax,Xdir,Xdel)
%Determine conditions that do not qualify for a search
if Xmax == 0
    Plen = 0;
    return
elseif Xdel > 0
    Plen = 0;
    return
elseif Xdir == 5 && X0 == 1
    Plen = 0;
    return
elseif Xdir == 3 && X0 == size(Xalign,2)
    Plen = 0;
    return    
end

%If direction is 5' , easier to just flip the alignment and position, calc
%Plen same way
if Xdir == 5
    Xalign = fliplr(Xalign);
    X0 = size(Xalign,2) - X0+1;
end

%Determine the sample and reference NT evaluation regions
SeqEnd = X0 + Xmax;
if SeqEnd > size(Xalign,2)
    SeqEnd = size(Xalign,2);
end

RefStart = X0 - Xmax + 1;
if RefStart < 1
    RefStart = 1;
end

%Extract the NT evaluation region, but ensuring it's just ACGTU
SeqEval = Xalign(1,X0+1:SeqEnd);
RefEval = seqrcomplement(Xalign(3,RefStart:X0));
MaxValidMatch = min([length(RefEval) length(SeqEval) (regexpi(SeqEval,'[^ACGTU]*','once') - 1)]); %Find shortest distance to non-NT value. -1 to get to max valid NT position.
if isempty(MaxValidMatch) == 0
    Xmax = MaxValidMatch;
else
    Xmax = MaxValidMatch;
end
SeqEval = SeqEval(1:Xmax);
RefEval = RefEval(1:Xmax);

%Find continuous matches in the p region.
Plen = 0;
while Plen < length(RefEval)
    if strcmpi(RefEval(Plen+1),SeqEval(Plen+1))
        Plen = Plen+1;
    else
        break
    end
end    