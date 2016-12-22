%getUnqStrain will look at the VDJ database to look for unique mouse
%strains. Returns a cell with  unique strains.
%
%  UnqStrain = getUnqStrain(Vmap,Dmap,Jmap)
%
%  UnqStrain = getUnqStrain(Vmap,Dmap,Jmap,'condense',4) will condense the
%  list if the first N number of letter are the same accross the strains
%  EX:  {'BALB.K' 'BALB/b' 'BALB/c'} will become a single entry for
%  'BALB.K,BALB/b,BALB/c'

function UnqStrain = getUnqStrain(Vmap,Dmap,Jmap,varargin)
P = inputParser;
addParameter(P,'Condense',0,@isnumeric);
parse(P,varargin{:});
Condense = P.Results.Condense;

[~,~,~,Header] = getCurrentDatabase;
HostLoc = findCell(Header,'isolatedHost');

%concantenate string
StrainName = '';
for k = 1:3
    switch k
        case 1
            Xmap = Vmap;
        case 2
            Xmap = Dmap;
        case 3
            Xmap = Jmap;
    end
    for x = 1:size(Xmap,1)
        StrainName = [StrainName Xmap{x,HostLoc} ';'];
    end
end
StrainCell = regexpi(StrainName,';','split');
UnqStrain = unique(StrainCell)';

%Remove any empty entries, which should be first
if isempty(UnqStrain{1})
    UnqStrain(1) = [];
end

if Condense > 0
    CurStr = ''; 
    DelThis = zeros(length(UnqStrain),1,'logical');
    for j = 1:length(UnqStrain)
        if length(UnqStrain{j}) < Condense
            CurStr = '';
            continue
        elseif length(UnqStrain{j}) >= Condense && isempty(CurStr)
            CurStr = UnqStrain{j};
            continue
        end
            
        if strcmpi(CurStr(1:Condense),UnqStrain{j}(1:Condense))
            UnqStrain{j} = [UnqStrain{j-1} ',' UnqStrain{j}];
            DelThis(j-1) = 1;
        else
            CurStr = UnqStrain{j};
        end
    end
    UnqStrain(DelThis) = [];
end
