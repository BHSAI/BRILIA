%saveSeqDataForBaseline will convert a VDJdata to a readable fasta file
%for BASELINe[REF].
%
%  saveSeqDataForBaseline(VDJdata, VDJheader, )
%
%  Fasta file with >>> sequence is a group separator, >> is the reference
%  sequence, and > is the clone sequence. 
%
%  Example: 
%    >>> Group 1
%    >> Germline1 Sequence
%    ACGCTGA
%    > Seq 1
%    ACGTTGA
%    > Seq 2
%    AGTTTGA
%
%REF: Gur Yaari; Mohamed Uduman; Steven H. Kleinstein. Quantifying
%selection in high-throughput Immunoglobulin sequencing data sets. Nucleic
%Acids Res. 2012 May 27.
%BaseLine v1.3 at http://selection.med.yale.edu/baseline/


function saveSeqDataForBaseline(FullSaveName, VDJdata, VDJheader, varargin)
[H, L, Chain] = getAllHeaderVar(VDJheader);
[SeqIMGT, RefSeqIMGT] = formatSeq2IMGT(VDJdata, VDJheader);
GrpNum = cell2mat(VDJdata(:, H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);
Data(1:(size(VDJdata,1)+length(UnqGrpNum))) = struct('Sequence','','Header','');

q = 1;
for y = 1:length(UnqGrpNum)
    Idx = find(GrpNum == UnqGrpNum(y));
    
    %Assemble a name for the germline 
    RefName = sprintf('>>Grp%d', UnqGrpNum(y));
    Data(q).Header = RefName;
    Data(q).Sequence = strrep(RefSeqIMGT{Idx(1), 1}, '.', 'N');
    q = q + 1;
    
    %Now add the rest
    for k = 1:length(Idx)
        SeqName = VDJdata{Idx(k), H.SeqNameLoc};
        if isnumeric(SeqName)
            SeqName = num2str(SeqName);
        end
        Data(q).Header = SeqName;
        Data(q).Sequence = strrep(SeqIMGT{k, 1}, '.', 'N');
        q = q + 1;
    end
end

fastawrite(FullSaveName, Data);
