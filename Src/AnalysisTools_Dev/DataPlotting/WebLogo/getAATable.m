%getAATable outputs the AA hydropathy index (HPI) and Van der Waals Volume
%(VdVW) data table. It also returns the RGB color code matrix as a 2nd
%output.
%
%  AATable = getAATable
%
%  [AATable, AAColor] = getAATable 
%
%  OUTPUT
%    AATable: list of AA and their properties. Each column store the
%      following:
%          Col 1:  3-letter AA
%          Col 2:  1-letter AA
%          Col 3:  Hydropathy Index (unitless)
%          Col 4:  Van der Waals Volume (A^3)
%          Col 5:  AA property group
%          Col 6:  Codons that encode the AA
%
%    AAColor: a 20x3 RGB matrix. Each AA is color coded by its property
%      group number, as discussed below:
%           Num     Color        Property        Residues
%           ----  -------- ------------------- ------------- 
%            1     Blue     Polar Posistive     H K R
%            2     Red      Polar Negative      D E
%            3     Green    Polar Neutral       N Q S T 
%            4     Gray     NonPolar Aliphatic  A I L M V
%            5     Magenta  NonPolar Aromatic   F W Y
%            6     Brown    Unique              G P
%            7     Cyan     Cysteine Bond       C
%
%  Ref for HPI : J Kyte & RF Doolittle. 1982. J. Mol. Biol. 157:105-132.
%  Ref for VdwV: FM Richards. 1974. J. Mol. Biol. 82:1-14.
function varargout = getAATable

AATable = {
    'Ala'    'A'    [ 1.8000]    [ 67]    [4]   {'GCT' 'GCC' 'GCA' 'GCG'}
    'Arg'    'R'    [-4.5000]    [148]    [1]   {'CGT' 'CGC' 'CGA' 'CGG' 'AGA' 'AGG'}
    'Asn'    'N'    [-3.5000]    [ 96]    [3]   {'AAT' 'AAC'}
    'Asp'    'D'    [-3.5000]    [ 91]    [2]   {'GAT' 'GAC'}
    'Cys'    'C'    [ 2.5000]    [ 86]    [7]   {'TGT' 'TGC'}
    'Gln'    'Q'    [-3.5000]    [114]    [3]   {'CAA' 'CAG'}
    'Glu'    'E'    [-3.5000]    [109]    [2]   {'GAA' 'GAG'}
    'Gly'    'G'    [-0.4000]    [ 48]    [6]   {'GGT' 'GGC' 'GGA' 'GGG'}
    'His'    'H'    [-3.2000]    [118]    [1]   {'CAT' 'CAC'}
    'Ile'    'I'    [ 4.5000]    [124]    [4]   {'ATT' 'ATC' 'ATA'}
    'leu'    'L'    [ 3.8000]    [124]    [4]   {'TTA' 'TTG' 'CTT' 'CTC' 'CTA' 'CTG'}
    'Lys'    'K'    [-3.9000]    [135]    [1]   {'AAA' 'AAG'}
    'Met'    'M'    [ 1.9000]    [124]    [4]   {'ATG'}
    'Phe'    'F'    [ 2.8000]    [135]    [5]   {'TTT' 'TTC'}
    'Pro'    'P'    [ 1.6000]    [ 90]    [6]   {'CCT' 'CCC' 'CCA' 'CCG'}
    'Ser'    'S'    [-0.8000]    [ 73]    [3]   {'TCT' 'TCC' 'TCA' 'TCG' 'AGT' 'AGC'}
    'Thr'    'T'    [-0.7000]    [ 93]    [3]   {'ACT' 'ACC' 'ACA' 'ACG'}
    'Trp'    'W'    [-0.9000]    [163]    [5]   {'TGG'}
    'Tyr'    'Y'    [-1.3000]    [141]    [5]   {'TAT' 'TAC'}
    'Val'    'V'    [ 4.2000]    [105]    [4]   {'GTT' 'GTC' 'GTA' 'GTG'}
    };

%Re-sort AA prop to match MatLab order: ARNDCQEGHILKMFPSTWYV
AAorder = cell2mat(AATable(:,2))';
AAsortnum = aa2int(AAorder);
AATable(AAsortnum',:) = AATable;

if nargout == 2
    %Assigning a color scheme to each group, [R G B]
    CLR(1,:) = [0.0 0.0 1.0]; %Blue
    CLR(2,:) = [1.0 0.0 0.0]; %Red 
    CLR(3,:) = [0.0 0.6 0.0]; %Green
    CLR(4,:) = [0.3 0.3 0.3]; %Gray
    CLR(5,:) = [1.0 0.0 1.0]; %Magenta
    CLR(6,:) = [0.6 0.3 0.2]; %Brown
    CLR(7,:) = [0.0 0.8 0.8]; %Cyan
    GrpN = cell2mat(AATable(:,5));
    AAColor = CLR(GrpN,:);
end

%Return the output
if nargout >= 1
    varargout{1} = AATable;
    if nargout >= 2
        varargout{2} = AAColor;
    end
end