%fixPseudoVDJ will look through V annotations, and then fix those that
%have multiple V gene suggestions, BUT, some are pseudogenes. Pseudogenes
%are removed by default, and the first non-pseudogene annotation gets
%taken. This does not work on single-suggestions. Only for V gene, as it's
%possible for a D J pseudogene to become a normal gene after enough NT
%deletions.
%
%This will also look for pseudo reference seq and

function VDJdata = fixPseudoVDJ(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end       
getHeaderVar;

%Header for the Vmap;
Header = {'nucleotide' 'GeneNameIMGT' 'STDgeneName' 'geneSubgroup' 'geneName' 'geneAllele' 'function' 'readingFrame' 'isolatedHost' 'keyNTloc'};
MapFunctLoc = findHeader(Header,'function');

%Identify location of pseudo genes
PseudoLoc = zeros(size(Vmap,1),1,'logical');
for v = 1:size(Vmap,1)
    if ~isempty(regexpi(Vmap{v,MapFunctLoc},'P|ORF'))
        PseudoLoc(v) = 1;
    end
end

%Iteratively find groups
GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    %Find the sequences per grp, and the Vnum
    IdxLoc = find(UnqGrpNum(y) == GrpNum);   
    Vnum = VDJdata{IdxLoc(1),FamNumLoc(1)};
    
    %Look for any pseudo V (ignore pseudo D and J b/c too uncertain)
    DelThis = zeros(1,length(Vnum),'logical');
    for k = 1:length(Vnum)
        if PseudoLoc(Vnum) == 1
            DelThis(k) = 1;
        end
    end
    HavePseudo = max(DelThis);
    HowManyPseudo = sum(DelThis);
    
    if HavePseudo %Figure out what to do
        if HowManyPseudo > 1 && (HowManyPseudo ~= length(Vnum)) %Must have at least 1 instance, but not all instances
            %Need to adjust the name
            KeepThis = DelThis == 0;
            TempCell = cell(sum(KeepThis),6);
            TempCell(:,1) = num2cell(Vnum(KeepThis))';
            TempCell(:,2) = Vmap(Vnum(KeepThis),3);
            TempResults = reduceFamily(TempCell);

            VDJdata(IdxLoc,[FamNumLoc(1) FamLoc(1)]) = repmat(TempResults(1,1:2),length(IdxLoc),1);
            disp('Reduced Pseudo');
        end
    end
end

%Fill in the details now
VDJdata = buildRefSeq(VDJdata,NewHeader,'germline','single'); %must do singles, since group therapy not done.
VDJdata = buildVDJalignment(VDJdata,NewHeader,Vmap,Dmap,Jmap); %Alignment Info
VDJdata = makeClassifier(VDJdata,NewHeader); %Classifier + FormattedSeq
VDJdata = appendMutCt(VDJdata,NewHeader); %SHM infor on the VMDNJ segments
VDJdata = findCDR3(VDJdata,NewHeader); %Get the CDR3 seq and info 
VDJdata = labelNonprodVDJ(VDJdata,NewHeader);            