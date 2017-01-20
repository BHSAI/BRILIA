%filterRefGene will process the Vmap, Dmap, and Jmap to retain only those
%of interest.
%
%  [Vmap,Dmap,Jmap] = filterRefGene(Vmap,Dmap,Jmap) will prompt users to
%  select what sequences to keep. Discarded entries have the sequence set
%  to '' empty.
%
%  [Vmap,Dmap,Jmap] =
%  filterRefGene(Vmap,Dmap,Jmap,Parameter,Value,Parmaeter2,Value2,...)
%  will set the VDJ maps sequence column to empty '' based on what is
%  being filtered.
%    
%    Parameter      Value (pick one)
%    Strain         all, C57BL, BALB, etc... 
%    Ddirection     all, Dfwd, Dinv
%    Vfunction      all, F, P, ORF
%    KeepThis       yes (will keep these) or no (will delete these)
%
%    IMPORTANT!  
%       - Function F,P,ORF only works for V genes. 
%       - Direction only works for D genes. 
%       - Strain only works for mouse database.

function [Vmap,Dmap,Jmap] = filterRefGene(Vmap,Dmap,Jmap,varargin)
%Determine the inputs
P = inputParser;
addParameter(P,'Strain','',@ischar); %For mouse, select strain(s)
addParameter(P,'Ddirection','',@ischar); %Inv, fwd, or all direction
addParameter(P,'Vfunction','',@ischar); %F for functional, P for pseudo, ORF for ORF (EX: ORF,P,F)
addParameter(P,'KeepThis','yes',@ischar)%yes will keep what the strain info, while no will do opposite and remove them instead.
parse(P,varargin{:});

Strain = P.Results.Strain;
Ddirection = P.Results.Ddirection;
Vfunction = P.Results.Vfunction;
KeepThis = P.Results.KeepThis;

%Check current database for the name and header file
[~,~,~,Header,DBname] = getCurrentDatabase; %Header = {'nucleotide' 'GeneNameIMGT' 'STDgeneName' 'geneSubgroup' 'geneName' 'geneAllele' 'function' 'readingFrame' 'isolatedHost' 'keyNTloc'};
[~,DBname,~] = parseFileName(DBname);

%See if the database is mouse. If it is, users will have to choose
IsMouseDB = 0;
if ~isempty(regexpi(DBname,'mouse'))
    IsMouseDB = 1;
end

%Set up the deletion index matrix
Vdel = zeros(size(Vmap,1),1,'logical');
Ddel = zeros(size(Dmap,1),1,'logical');
Jdel = zeros(size(Jmap,1),1,'logical');

%==========================================================================
%Select by mouse strain, if not specified
if IsMouseDB == 1
    HostLoc = findCell(Header,'isolatedHost');
    UnqStrain = [{'All'}; getUnqStrain(Vmap,Dmap,Jmap,'Condense',4)];
    if isempty(Strain) %Ask users to select desired strain
        disp('What mouse strain is it? ')
        dispList(UnqStrain)
        StrainNum = -1;
        while StrainNum <1 || StrainNum > length(UnqStrain)
            StrainNum = input('Select Option: ');
            if isempty(StrainNum)
                StrainNum = 1;
                break
            end
        end
    else %Determine which unq strain matches with user input
        StrainNum = findCell(UnqStrain,Strain,'MatchCase','Any','MatchWord','Partial');
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

    if StrainNum > 1 %Need to filter, since it's not ALL
        SearchFor = strrep(UnqStrain{StrainNum},',','|');
        for k = 1:3
            switch k
                case 1
                    Xmap = Vmap;
                case 2
                    Xmap = Dmap;
                case 3
                    Xmap = Jmap;
            end

            DelThis = zeros(size(Xmap,1),1,'logical');
            for j = 1:size(Xmap,1)
                if isempty(regexpi(Xmap{j,HostLoc},SearchFor,'once'))
                    DelThis(j) = 1;
                end
            end

            switch k
                case 1
                    Vdel = Vdel | DelThis;
                case 2
                    Ddel = Ddel | DelThis;
                case 3
                    Jdel = Jdel | DelThis;
            end    
        end
    end
end

%==========================================================================
%Determine if you want to search Dinverse too
NameLoc = findCell(Header,'GeneNameIMGT');
DsearchOptions = {'All'; 'Dfwd'; 'Dinv'};

if isempty(Ddirection) %Ask user to determine which direction D to search
    disp('Which D gene direction do you want?');
    dispList(DsearchOptions);
    SearchDopt = input('Select option: ');
    if isempty(SearchDopt)
        SearchDopt = 1;
    end
    while SearchDopt < 1 || SearchDopt > 3
        SearchDopt = input('Select option: ');
    end
else
    SearchDopt = findCell(DsearchOptions,Ddirection,'MatchCase','any','MatchWord','partial');
end

switch SearchDopt
    case 2 %Keep only Dfwd by deleting inverse
        for j = 1:size(Dmap,1)
            if Dmap{j,NameLoc}(1) == 'r'
                Ddel(j) = 1;
            end
        end
    case 3 %Keep only Dinv by deleting forward
        for j = 1:size(Dmap,1)
            if Dmap{j,NameLoc}(1) ~= 'r'
                Ddel(j) = 1;
            end
        end
    otherwise %Don't do anything
end

%==========================================================================
%Determine if you want to include pseudo or orf sequences
FunctLoc = findCell(Header,'function');
FunctOptions = {'All'; 'F'; 'P'; 'ORF'};

if isempty(Vfunction) %Ask user to determine which direction D to search
    disp('Which V gene functionality do you want? Can list multiples, separated by comma (EX: 2,3)');
    dispList(FunctOptions);
    FunctOpt = -1;
    while min(FunctOpt) < 1 || max(FunctOpt) > length(FunctOptions)
        FunctOpt = input('Select option: ','s');
        if min(isstrprop(strrep(FunctOpt,',',''),'digit')) == 1
            FunctOpt = eval(['[' FunctOpt ']']);
        else
            FunctOpt = -1;
        end
    end
else
    Vfunction = regexp(Vfunction,',','split');
    FunctOpt = findCell(FunctOptions,Vfunction);
end

if min(FunctOpt(1)) > 1 && length(FunctOpt) > 1 %You have to filter
    for j = 1:size(Vmap,1)
        CurFunct = Vmap{j,FunctLoc};
        ParLoc = regexp(CurFunct,'\(|\)|\[|\]');
        CurFunct(ParLoc) = [];
        WantThis = 0;
        for k = 1:length(FunctOpt)
            EvalFunct = FunctOptions{FunctOpt(k)};
            if length(CurFunct) == length(EvalFunct) && strcmpi(CurFunct,EvalFunct)
                WantThis = 1;
                break
            end
        end
        if WantThis == 0
            Vdel(j) = 1;
        end
    end
end

%==========================================================================
%Now set all nuceotide columns to empty
SeqLoc = findCell(Header,'nucleotide');

if strcmpi(KeepThis,'no')
    Vdel = Vdel == 0;
    Ddel = Ddel == 0;
    Jdel = Jdel == 0;
end
Vmap(Vdel,SeqLoc) = {''};
Dmap(Ddel,SeqLoc) = {''};
Jmap(Jdel,SeqLoc) = {''};

%Make sure there is at least some entries left for V,D,J
if sum(Vdel) == size(Vmap,1)
    disp('Nothing left in Vmap');
end
if sum(Ddel) == size(Dmap,1)
    disp('Nothing left in Dmap');
end
if sum(Jdel) == size(Jmap,1)
    disp('Nothing left in Jmap');
end
