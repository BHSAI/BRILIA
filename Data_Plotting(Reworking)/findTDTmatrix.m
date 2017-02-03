%findTDTmatrix will look at the N1 and N2 compositions, and then extract
%the 4x4 adjacency matrix of the nucleotide. Each entry tells what is the
%probability C is next A,C,G,T. 
%
%  [TDTmatrix, TDTcomp] = findTDTmatrix(VDJdata,NewHeader) will return the
%  matrix revealing the likelihood X nt is next to Y nt. TDTmatrix is a
%  normalized, diagonally symmetric matrix. TDTcomp reveals frequency count
%  of a each nt.
%
%  findTDTmatrix(VDJdata,NewHeader,'divide') does the same thing as above,
%  but will attempt to divide the N regions into a 5' and 3' side, and then
%  find the stats using 5' and complement 3' strand N region. This is
%  because TDT could add to complement strand G and A's which translates to
%  C and T's on the coding strand. 
% 
%  Note: Avoids P nucleotides by default, and will attempt to find
%  consensus N regions for each group of sequence.

function [TDTmatrix, varargout] = findTDTmatrix(varargin)
%Parse the input
P = inputParser;
addOptional(P,'VDJdata',{},@iscell);
addOptional(P,'NewHeader',{},@iscell);
addParameter(P,'ScanThese','MN',@ischar);
addParameter(P,'LeftNTs','AG',@ischar);
addParameter(P,'RightNTs','CT',@ischar);
addParameter(P,'Mode','single',@(x) any(validatestring(x,{'single','divide','flip'})));
parse(P,varargin{:});

VDJdata = P.Results.VDJdata;
NewHeader = P.Results.NewHeader;
ScanThese = P.Results.ScanThese;
LeftNTs = P.Results.LeftNTs;
RightNTs = P.Results.RightNTs;
Mode = P.Results.Mode;

%Open file if needed
if isempty(VDJdata)
    [VDJdata,NewHeader] = openSeqData;
end
getHeaderVar;

GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
GrpNumUnq = unique(GrpNum);
IdxMap = 1:size(GrpNum,1);

TDTmatrix = zeros(4); %Probability X nt and Y nt are adjacent, within the N regions.
TDTcomp = zeros(1,4); %Composition of nts in the N region as ACGT
Nleft = 0;
Nright = 0;
Ndiv = 0;
TotFlip = 0;
UncertainFlip = 0;
for y = 1:length(GrpNumUnq)
    %Always process per group, otherwise you'll get clonal expansion bias issues
    GrpLoc = (GrpNumUnq(y) == GrpNum);
    IdxLoc = IdxMap(GrpLoc);
    j = IdxLoc(1); %First one is fine
    
    %Extract relevant information from VDJdata
    RefSeq = VDJdata{IdxLoc(1),RefSeqLoc};
    Classifier = VDJdata{IdxLoc(1),FormClassLoc(2)};
    
    for q = 1:length(ScanThese)   
        Mloc = regexpi(Classifier,ScanThese(q));
        Mseq = RefSeq(Mloc);
        BadNTloc = regexpi(Mseq,'[^ACGTU]');
        Mseq(BadNTloc) = [];
        Mct = length(Mseq);
        if Mct == 0
            continue %No N region.
        end

        
        switch lower(Mode)
            case 'divide'
                [LeftSide, RightSide] = divideNregion(Mseq,LeftNTs,RightNTs);
                if ~isempty(RightSide)
                    Mseq = [LeftSide seqcomplement(RightSide)];
                    Nleft = Nleft + length(LeftSide);
                    Nright = Nright + length(RightSide); 
                    Ndiv = Ndiv + 1;
                end                
            case 'flip'
                [Mseq, Flipped] = flipNregion(Mseq,LeftNTs,RightNTs); %Flipped = -1 for unknoqn, 0 for no flip, and 1 for flipped.
                if Flipped < 0 %Equal number of LeftNTs and RightNTs
                    UncertainFlip = UncertainFlip + abs(Flipped);
                else
                    TotFlip = TotFlip + Flipped;
                end
                Ndiv = Ndiv + 1;
            case 'single' %Don't do anything
        end
        
        Mint = nt2int(Mseq);
        TDTcomp = TDTcomp + histc(Mint,1:4);

        if Mct >= 2
            AdjCoor = [Mint(1:end-1); Mint(2:end)];
            for n = 1:size(AdjCoor,2)
                %You only want a diagonal matrix, so shove lower number top.
                C1 = AdjCoor(1,n);
                C2 = AdjCoor(2,n);
                if C1 <= 4 && C2 <= 4
                    TDTmatrix(C1,C2) = TDTmatrix(C1,C2) + 1;
                    TDTmatrix(C2,C1) = TDTmatrix(C2,C1) + 1;
                end
            end
        end       
    end
end

if nargout >= 2
    varargout{1} = TDTcomp;
    if nargout == 3
        if strcmpi(Mode,'divide')
            varargout{2} = [Nleft/Ndiv Nright/Ndiv];
        elseif strcmpi(Mode,'flip')
            varargout{2} = [1-TotFlip/Ndiv TotFlip/Ndiv UncertainFlip/Ndiv];
        else
            varargout{2} = 0;
        end
    end
end