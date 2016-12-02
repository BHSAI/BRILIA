%filterRefGene will process the Vmap, Dmap, and Jmap to retain only those
%of interest, mainly including functional genes or strain information.

function [Vmap,Dmap,Jmap] = filterRefGene(Vmap,Dmap,Jmap,varargin)
Header = {'nucleotide' 'GeneNameIMGT' 'STDgeneName' 'geneSubgroup' 'geneName' 'geneAllele' 'function' 'readingFrame' 'isolatedHost' 'keyNTloc'};
Vdel = zeros(size(Vmap,1),1) > 1;
Ddel = zeros(size(Dmap,1),1) > 1;
Jdel = zeros(size(Jmap,1),1) > 1;

%Select by strain==========================================================
HostLoc = findHeader(Header,'isolatedHost');

if ~isempty(varargin)
    StrainNum = varargin{1};
else
    StrainNum = -1;
end
while StrainNum < 0 || StrainNum > 8
    %Select the desired strain
    disp('0) Any Strain');
    disp('1) A/J');
    disp('2) BALB/c or .K');
    disp('3) C57BL/6 or 6J or 10 (DEFAULT)');
    disp('4) I/St');
    disp('5) MRL/1pr');
    disp('6) NFS');
    disp('7) NZB');
    disp('8) 129/Sv');
    StrainNum = input('Select the Mus Musculus mouse strain: ');
    if isempty(StrainNum)
        StrainNum = 3;
    end
end

if StrainNum > 0
    %Select the key word to search for when filtering data
    switch StrainNum
        case 1
            SearchFor = 'F';
        case 2
            SearchFor = 'P';
        case 3
            SearchFor = 'C57BL';
        case 4
            SearchFor = 'I';
        case 5
            SearchFor = 'MRL';
        case 6
            SearchFor = 'NFS';
        case 7
            SearchFor = 'NZB';
        case 8
            SearchFor = '129';
    end

    for k = 1:3
        switch k
            case 1
                Xmap = Vmap;
            case 2
                Xmap = Dmap;
            case 3
                Xmap = Jmap;
        end

        DelThis = zeros(size(Xmap,1),1)>1;
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

%Determine if you want to search Dinverse too=============================
NameLoc = findHeader(Header,'GeneNameIMGT');

while 1
    SearchDinv = input('Do you want to search Dinverses too? y(Default) or n ','s');
    if strcmpi(SearchDinv,'y') || isempty(SearchDinv);
        break
    elseif strcmpi(SearchDinv,'n');        
        for j = 1:size(Dmap,1)
            if sum(regexpi(Dmap{j,NameLoc},'r')) > 0
                Ddel(j) = 1;
            end
        end
        break
    end
end

Vmap(Vdel,1) = {''};
Dmap(Ddel,1) = {''};
Jmap(Jdel,1) = {''};     



