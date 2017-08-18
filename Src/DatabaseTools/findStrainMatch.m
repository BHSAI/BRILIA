%findStrainMatch will compare a gene name and identify all host strains
%that have this gene, based on IMGT's gene table that can be found at
%http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire
%=genetable&species=Mus_musculus&group=IGHD, but pasted into Excel and
%formatted with reformatGeneTable.m. 
%
%  AllStrains = findStrainMatch(GeneName,GeneTable)
%
%  INPUT
%    GeneName: The full IMGT gene name that you are looking for strain info
%    GeneTable: The gene table downloaded from IMGT, inlcuding the header
%
%  OUTPUT
%    AllStrains: string of all strains of a host that has the GeneName,
%      sepeartaed by ";"
%
%  NOTE
%    To get the gene table, go the gene table website in IMGT (see above)
%    and then copy and paste the table. Use reformatGeneTable.m to format the
%    tables to a more standardized format.
%
%  See also reformatGeneTable

function AllStrains = findStrainMatch(GeneName, GeneTable)
NameLoc = findCell(GeneTable(1, :), {'IMGT full name', 'IMGT allele name'}, 'MatchCase', 'Any');
StrainLoc = findCell(GeneTable(1, :), 'Strain');

%Look through the gene table for host where nts were isolated from
EvalLoc = findCell(GeneTable(:, NameLoc), GeneName);
EvalLoc(EvalLoc == 1) = []; %Can't have header as the selection.
if EvalLoc(1) > 0
    %Get all the strains
    Strains = cell(1, length(EvalLoc)*length(StrainLoc)); %Since there are > 1 column for strains
    DelLoc = zeros(size(Strains), 'logical');
    q = 1;
    for k = 1:length(StrainLoc)
        for w = 1:length(EvalLoc)
            Strains(q) = GeneTable(EvalLoc(w), StrainLoc(k));
            if isempty(Strains{q})
                DelLoc(q) = 1;
            end
            q = q+1;
        end            
    end
    Strains(DelLoc) = [];
    
    %Get the unique strains, or say Cloned otherwise
    UnqStrains = unique(Strains);
    if isempty(UnqStrains)
        UnqStrains = {'Cloned'};
    end
    
    %Combine all strains
    RepStr = repmat('%s;', 1, length(UnqStrains));
    RepStr(end) = [];
    AllStrains = sprintf(RepStr, UnqStrains{:});
else
    AllStrains = 'None';
end
