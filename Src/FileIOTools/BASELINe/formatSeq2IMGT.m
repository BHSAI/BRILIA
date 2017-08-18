%formatSeq2IMGT will convert a sequence to IMGT format if the CDR regions
%are defined. Will fill in unknown nts as N.
%
%  SeqIMGT = formatSeq2IMGT(VDJdata, VDJheader)
%
%  SeqIMGT = formatSeq2IMGT(Seq, CDR1s, CDR1e, CDR2s, CDR2e, CDR3s, CDR3e)

function [SeqIMGT, RefSeqIMGT] = formatSeq2IMGT(varargin)

if nargin == 2 && iscell(varargin{1}) && iscell(varargin{2}) %Has VDJdata dn VDJheader
    VDJdata = varargin{1};
    VDJheader = varargin{2};
    [H, L, Chain] = getAllHeaderVar(VDJheader);
    SeqIMGT = cell(size(VDJdata,1), length(Chain));
    RefSeqIMGT = cell(size(VDJdata,1), length(Chain));
    for k = 1:length(Chain)
        if strcmpi(Chain(k),'H')
            B = H;
        else
            B = L;
        end
        SeqLoc = B.SeqLoc;
        RefSeqLoc = B.RefSeqLoc;
        CDRLoc = zeros(1,6);
        if length(B.CDR1Loc) == 4
            CDRLoc(1) = B.CDR1Loc(3);
            CDRLoc(2) = B.CDR1Loc(4);
        end
        if length(B.CDR2Loc) == 4
            CDRLoc(3) = B.CDR2Loc(3);
            CDRLoc(4) = B.CDR2Loc(4);
        end
        if length(B.CDR3Loc) == 4
            CDRLoc(5) = B.CDR3Loc(3);
            CDRLoc(6) = B.CDR3Loc(4);
        end
        for j = 1:size(VDJdata,1)
            InputCDRCell = num2cell(zeros(1, 6));
            for q = 1:6
                if CDRLoc(q) ~= 0
                    InputCDRCell{q} = VDJdata{j, CDRLoc(q)};
                end
            end
            SeqIMGT{j, k} = format_Seq2IMGT(VDJdata{j, SeqLoc}, InputCDRCell{:});
            RefSeqIMGT{j, k} = format_Seq2IMGT(VDJdata{j, RefSeqLoc}, InputCDRCell{:});
        end
    end
elseif nargin == 7
    Seq = varargin{1};
    CDR1s = varargin{2};
    CDR1e = varargin{3};
    CDR2s = varargin{4};
    CDR2e = varargin{5};
    CDR3s = varargin{6};
    CDR3e = varargin{7};
    SeqIMGT = format_Seq2IMGT(Seq, CDR1s, CDR1e, CDR2s, CDR2e, CDR3s, CDR3e);
    RefSeqIMGT = {};
else
    error('%s: Wrong number of input.', mfilename);
end

function FormattedSeq = format_Seq2IMGT(Seq, CDR1s, CDR1e, CDR2s, CDR2e, CDR3s, CDR3e)

%IMGT AA position format
imgtCDR1s = 27;
imgtCDR1e = 38;
imgtCDR2s = 56;
imgtCDR2e = 65;
imgtCDR3s = 105;

%Start with CDR3s
FormattedSeq = repmat('.', 1, 104*3+(CDR3e-CDR3s+1));
if CDR3s > 0 && CDR3e > CDR3s
    s1 = (imgtCDR3s - 1)*3 + 1 ;
    s2 = (imgtCDR3s - 1)*3 + (CDR3e - CDR3s + 1);
    FormattedSeq(s1:s2) = Seq(CDR3s:CDR3e);
end

if CDR3s > 0 && CDR2e > 0
    s1 = imgtCDR2e*3 + 1;
    s2 = imgtCDR2e*3 + (CDR3s-1) - (CDR2e+1) + 1;
    FormattedSeq(s1:s2) = Seq(CDR2e+1:CDR3s-1);
else
    s1 = imgtCDR2e*3 + 1;
    s2 = imgtCDR2e*3 + (CDR3s-1);
    FormattedSeq(s1:s2) = Seq(1:CDR3s-1);
    return;
end

if CDR2s > 0 && CDR2e > CDR2s
    s1 = (imgtCDR2s - 1)*3 + 1 ;
    s2 = (imgtCDR2s - 1)*3 + (CDR2e - CDR2s + 1);
    FormattedSeq(s1:s2) = Seq(CDR2s:CDR2e);
else
    s1 = (imgtCDR2s - 1)*3 + 1 ;
    s2 = (imgtCDR2s - 1)*3 + CDR2e;
    FormattedSeq(s1:s2) = Seq(1:CDR2e);
    return;
end

if CDR2s > 0 && CDR1e > 0
    s1 = imgtCDR1e*3 + 1;
    s2 = imgtCDR1e*3 + (CDR2s-1) - (CDR1e+1) + 1;
    FormattedSeq(s1:s2) = Seq(CDR1e+1:CDR2s-1);
else
    s1 = imgtCDR1e*3 + 1;
    s2 = imgtCDR1e*3 + (CDR2s-1);
    FormattedSeq(s1:s2) = Seq(1:CDR2s-1);
    return;
end

if CDR1s > 0 && CDR1e > CDR1s
    s1 = (imgtCDR1s - 1)*3 + 1 ;
    s2 = (imgtCDR1s - 1)*3 + (CDR1e - CDR1s + 1);
    FormattedSeq(s1:s2) = Seq(CDR1s:CDR1e);
else
    s1 = (imgtCDR1s - 1)*3 + 1 ;
    s2 = (imgtCDR1s - 1)*3 + CDR1e;
    FormattedSeq(s1:s2) = Seq(1:CDR1e);
    return;
end

if CDR1s > 0
    s1 = 1;
    s2 = CDR1s - 1;
    FormattedSeq(s1:s2) = Seq(1:CDR1s-1);
end
