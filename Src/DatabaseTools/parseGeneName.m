%parseGeneName will take in a gene name, and return the gene family, gene
%number, and gene allele as a Mx3 cell matrix.
%
%  [FamNum, GeneNum, AlleleNum] = parseGeneName(AllGeneNames)
%
%  INPUT
%    AllGeneNames: char or Mx1 cell matrix of gene names
%  
%  OUTPUT
%    FamNum: the IG[HLK][VDJ]## number. Ex: 'IGHV1';
%    GeneNum: the string between the FamNum and AlleleNum. EX: '-12'
%    AlleleNum: the string after any asterick. EX: '*01'

function varargout = parseGeneName(AllGeneNames)
if ischar(AllGeneNames)
    AllGeneNames = {AllGeneNames};
end

FamNum = cell(size(AllGeneNames));
GeneNum = cell(size(AllGeneNames));
AlleleNum = cell(size(AllGeneNames));

for j = 1:length(AllGeneNames)
    GeneName = AllGeneNames{j};

    %Determine the IG family number
    Pat1 = '\w*IG[HLK][VDJ]\d*';
    FamNumT = regexpi(GeneName, Pat1, 'match');
    FamNumT = FamNumT{end};
    GeneName = strrep(GeneName, FamNumT, '');

    %Get the allele number + any appended custom name (EX IGHV1-1*01-mouse)
    Pat2 = '\*.+';
    AlleleNumT = regexpi(GeneName, Pat2, 'match');
    if isempty(AlleleNumT)
        AlleleNumT = {''};
    end
    AlleleNumT = AlleleNumT{end};

    %Get the gene number, which is whatever is left
    GeneNumT = strrep(GeneName, AlleleNumT, '');
    
    %Add all variables to the cell
    FamNum{j} = FamNumT;
    GeneNum{j} = GeneNumT;
    AlleleNum{j} = AlleleNumT;   
end

if nargout >= 1
    varargout{1} = FamNum;
    if nargout >= 2
        varargout{2} = GeneNum;
        if nargout >= 3
            varargout{3} = AlleleNum;
        end
    end
end
