%convertGermline2Fasta will get the current VDJmaps and output the fasta
%format.
function convertGermline2Fasta
[Vmap, Dmap, Jmap] = getCurrentDatabase;
[Vmap, Dmap, Jmap] = filterRefGene(Vmap,Dmap,Jmap);

SeqLoc = 1;
HeaderLoc = 2;
[SaveName, SavePath] = uiputfile('*.fa','Select the save name');
DotLoc = regexp(SaveName,'\.');
SaveNamePre = SaveName(1:DotLoc(end)-1);

for j = 1:3
    switch j
        case 1
            Seq = Vmap(:,SeqLoc);
            Header = Vmap(:,HeaderLoc);
            FileNameAddon = 'IGHV';
        case 2
            Seq = Dmap(:,SeqLoc);
            Header = Dmap(:,HeaderLoc);
            FileNameAddon = 'IGHD';
        case 3
            Seq = Jmap(:,SeqLoc);
            Header = Jmap(:,HeaderLoc);
            FileNameAddon = 'IGHJ';
    end
    %Delete the empty ones
    DelThis = zeros(length(Seq),1)>1;
    for k = 1:length(Seq)
        if isempty(Seq{k})
            DelThis(k) = 1;
        end
    end
    Seq(DelThis) = [];
    Header(DelThis) = [];
    
    DataStruct = cell2struct([Header Seq],{'Header' 'Sequence'},2);
        
    SaveName = [SaveNamePre FileNameAddon '.fa'];
    fastawrite([SavePath SaveName],DataStruct)  
    
    clear Seq Header
end