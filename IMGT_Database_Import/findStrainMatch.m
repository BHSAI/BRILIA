
function AllStrains = findStrainMatch(GeneName,Xtable,StrainLoc,NameLoc)
%Look through the gene table for host where nts were isolated from
EvalLoc = findHeader(Xtable(:,NameLoc),GeneName,'all');
if EvalLoc(end) > 0
    Strains = cell(1,length(EvalLoc)*length(StrainLoc));
    DelLoc = zeros(size(Strains))>1;
    q = 1;
    for k = 1:length(StrainLoc)
        for w = 1:length(EvalLoc)
            Strains(q) = Xtable(EvalLoc(w),StrainLoc(k));
            if isempty(Strains{q})
                DelLoc(q) = 1;
            end
            q = q+1;
        end            
    end
    Strains(DelLoc) = [];
    UnqStrains = unique(Strains);
    if isempty(UnqStrains)
        UnqStrains = {'Cloned'};
    end
    %Combine all strains
    AllStrains = UnqStrains{1};
    for r = 2:length(UnqStrains)
        AllStrains = [AllStrains ';' UnqStrains{r}];
    end
else
    AllStrains = 'None';
end