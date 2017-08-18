%filterRefGene will process the database DB to retain only those of
%interest, as defined by the filter criteria.
%
%  DB = filterRefGene(DB,Param,Value)
%  
%  INPUT
%    DB: database structure from getCurrentDatabase.m
%    Param-Value paris are as follows
%      Parameter     Value
%      ------------  ------------------------
%      Strain        'All' 'C57BL' 'BALB', ...
%      Ddirection    'All' 'fwd' 'rev' (or 'inv')
%      Vfunction     'All' 'F' 'P' 'ORF' 
%
%  OUTPUT
%    DB: the filtered version of the input DB, where unused gene sequences
%      are set to empty ''. This prevents the shifting of Map numbering
%      schemes required by the software to perform.
 
function [DB, varargout] = filterRefGene(DB, varargin)
%Determine the inputs
P = inputParser;
addParameter(P, 'Strain', '', @ischar); %For mouse, select strain(s)
addParameter(P, 'Ddirection', '', @ischar); %Inv, fwd, or all direction
addParameter(P, 'Vfunction', '', @ischar); %F for functional, P for pseudo, ORF for ORF (EX: ORF, P, F)
parse(P, varargin{:});
P = P.Results;

%Determine map fields only
Fields = fieldnames(DB);
MapLoc = findCell(Fields, 'map', 'MatchWord', 'partial');
Fields = Fields(MapLoc);
M = getMapHeaderVar(DB.MapHeader);

%==========================================================================
%If dealing with mouse database, select by mouse strain.
if ~isempty(regexpi(DB.FilePath, 'mouse'))
    UnqStrain = [{'All'}; getUnqStrain(DB, 4)];

    if isempty(P.Strain) %Ask users to select desired strain
        disp('What mouse strain is it? ')
        dispList(UnqStrain)
        StrainNum = -1;
        while StrainNum < 1 || StrainNum > length(UnqStrain)
            StrainNum = input('Select Option: ');
            if isempty(StrainNum)
                StrainNum = 1;
                break
            end
        end
    else %Determine which unq strain matches with user input
        StrainNum = findCell(UnqStrain, P.Strain, 'MatchCase', 'Any', 'MatchWord', 'Partial');
        if StrainNum == 0
            disp('Could not find match. Setting to ALL strain by default');
            StrainNum = 1;
        end
        if length(StrainNum) > 1 %Use must decide which one to pick
            disp('Multiple strains possible. What mouse strain is it? ')
            dispList(UnqStrain(StrainNum))
            StrainNumT = -1; 
            while StrainNumT < 1 || StrainNumT > length(StrainNum)
                StrainNumT = input('Select Option: ');
                if isempty(StrainNumT)
                    StrainNumT = 1;
                    break
                end
            end
            StrainNum = StrainNum(StrainNumT);
        end
    end
    
    %Filter database by strain
    if StrainNum > 1 
        SearchFor = strrep(UnqStrain{StrainNum}, ';', '|');
        for j = 1:length(Fields)
            Xmap = DB.(Fields{j});
            DelThis = zeros(size(Xmap, 1), 1, 'logical');
            for k = 1:size(Xmap, 1)
                if isempty(regexpi(Xmap{k, M.StrainLoc}, SearchFor, 'once'))
                    DelThis(k) = 1;
                end
            end
            Xmap(DelThis, M.SeqLoc) = {''};
            DB.(Fields{j}) = Xmap;
        end
        FiltOption.Strain = UnqStrain{StrainNum}; %Save the final filter value
    else
        FiltOption.Strain = 'All';  %Save the final filter value
    end
else
    FiltOption.Strain = 'All';  %Save the final filter value
end

%==========================================================================
%Determine if you want to search Dinverse too
DsearchOptions = {'All'; 'fwd'; 'rev'};

if isempty(P.Ddirection) %Ask user to determine which direction D to search
    disp('Which D gene direction do you want?');
    dispList(DsearchOptions);
    SearchDopt = input('Select option: ');
    if isempty(SearchDopt)
        SearchDopt = 1;
    end
    while SearchDopt < 1 || SearchDopt > 3
        SearchDopt = input('Select option: ');
    end
else %Search of user selection
    SearchDopt = findCell(DsearchOptions, P.Ddirection, 'MatchCase', 'any', 'MatchWord', 'partial');
end
FiltOption.Ddirection = DsearchOptions{SearchDopt}; %Save the final filter value

switch SearchDopt
    case 2 %Keep only Dfwd by deleting inverse
        for j = 1:size(DB.Dmap, 1)
            if DB.Dmap{j, M.GeneLoc}(1) == 'r'
                DB.Dmap{j, M.SeqLoc} = '';
            end
        end
    case 3 %Keep only Dinv by deleting forward
        for j = 1:size(DB.Dmap, 1)
            if DB.Dmap{j, M.GeneLoc}(1) ~= 'r'
                DB.Dmap{j, M.SeqLoc} = '';
            end
        end
    otherwise %Don't do anything
end

%==========================================================================
%Determine if you want to include pseudo or orf sequences
FunctOptions = {'All'; 'F'; 'P'; 'ORF'};

if isempty(P.Vfunction) %Ask user to determine which direction D to search
    disp('Which V gene functionality do you want? Can list multiple via comma (EX: 2, 3)');
    dispList(FunctOptions);
    FunctOpt = -1;
    while min(FunctOpt) < 1 || max(FunctOpt) > length(FunctOptions)
        FunctOpt = input('Select option: ', 's');
        if isempty(FunctOpt); FunctOpt = '1'; end
        if min(isstrprop(strrep(FunctOpt, ',', ''), 'digit')) == 1
            FunctOpt = eval(['[' FunctOpt ']']);
        else
            FunctOpt = -1;
        end
    end
else
    Vfunction = regexp(P.Vfunction, ',', 'split');
    FunctOpt = findCell(FunctOptions, Vfunction, 'MatchCase', 'any');
    if max(FunctOpt) == 0
        warning('Invalid P.Vfunction input. Using ''All'' option by default');
        FunctOpt = 1;
    end
end

if min(FunctOpt(1)) > 1 %You have to filter
    FieldIdx = findCell(Fields, 'V', 'MatchWord', 'Partial');
    if FieldIdx(1) == 0; FieldIdx = []; end
    for j = 1:length(FieldIdx)
        Xmap = DB.(Fields{FieldIdx(j)});
        for k = 1:size(Xmap, 1)
            CurFunct = upper(Xmap{k, M.FunctLoc});
            ParLoc = regexp(CurFunct, '\(|\)|\[|\]'); %Delete parenthesis and brackets
            CurFunct(ParLoc) = [];
            if ~ismember({CurFunct}, FunctOptions(FunctOpt))
                Xmap{k, M.SeqLoc} = '';
            end
        end
        DB.(Fields{FieldIdx(j)}) = Xmap;
    end
end

StrPat = repmat('%s, ', 1, length(FunctOpt));
StrPat(end) = [];
FiltOption.Vfunction = sprintf(StrPat, FunctOptions{FunctOpt});

%==========================================================================
if nargout >= 2
    varargout{1} = FiltOption;
end
