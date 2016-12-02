%getCurrentDatabase will look in the Database_Manager folder to look for
%the currently active database. This can also be used to switch the active
%database. Useful since many codes are searching for the current database,
%and this will make it easy to switch between host databases.
%
%  DatabaseFullPath = getCurrentDatabase  will look for the current
%  database folder in the Active directory.
%  
%  [Vmap, Dmap, Jmap] = getCurrentDatabase  will look for the current
%  database folder in the Active directory, and also load the VDJ
%  databases.
%
%  [...] = getCurrentDatabase('change') will prompt user to select a
%  database from the active or dormant folder. Once selected, the database
%  will move the active folder, while the others move the dormant folder.

function varargout = getCurrentDatabase(varargin)
%Find the directory for the active and dormant database *.mat files
Mpath = mfilename('fullpath');
SlashLoc = regexp(Mpath,'\\|\/');
SlashType = Mpath(SlashLoc(end));
FilePath = Mpath(1:SlashLoc(end));
ActiveDir = [FilePath 'Active' SlashType];
DormantDir = [FilePath 'Dormant' SlashType];

%Extract the database files
ActiveMat = dir([ActiveDir '*.mat']);
DormantMat = dir([DormantDir '*.mat']);

%If there is nothing in the active folder, look for database
if isempty(ActiveMat)
    SelectDatabase = 1;
else
    SelectDatabase = 0;
end

%If use specifiy to change database, look for database
if ~isempty(varargin)
    if strcmpi(varargin{1},'change')
        SelectDatabase = 1;
        SelectThisOne = ''; %If you know the exact file, select this

        %Acknowledge IMGT database references
        sprintf('Germline VDJ databases were downloaded from the IMGT, the international ImMunoGeneTics information system \n at http://www.imgt.org (founder and director: Marie-Paule Lefranc, Montpellier, France)')
    else
        SelectDatabase = 1;
        SelectThisOne = varargin{1};
    end
end

%If there is a need to change databse, do so now.
if SelectDatabase == 1
    %Ask use to choose what database to use
    while 1
        if isempty(SelectThisOne)
            if ~isempty(ActiveMat)
                disp(['0 or Enter) ' ActiveMat(1).name]);
            end
            for j = 1:length(DormantMat)
                disp([num2str(j) ') ' DormantMat(j).name]);
            end
            SelectNum = input('Select the database number to use');
            if ~isempty(SelectNum)
                if SelectNum > 0 && SelectNum <= length(DormantMat)
                    disp(['Using: ' DormantMat(SelectNum).name]);
                    break
                end
            elseif isempty(SelectNum)
                SelectNum = 0;
                disp(['Using: ' ActiveMat(1).name]);
                break
            else
                SelectNum = 0;
                disp(['Using: ' ActiveMat(1).name]);
                break            
            end
        else %You know what database to use, select this one.
            if ~isempty(strfind(ActiveMat(1).name,SelectThisOne))
                SelectNum = 0;
                break;
            else
                SelectNum = 0; %Just in case  you can't find one.
                for j = 1:length(DormantMat)
                    if strfind(DormantMat(j).name,SelectThisOne)
                       SelectNum = j;
                       break
                    end
                end
                break;
            end
        end
    end
    
    %Switch the database file now
    if ~isempty(ActiveMat) && SelectNum > 0
        for j = 1:length(ActiveMat) %There should only be 1 active mat. Move anything else, just in case.
            movefile([ActiveDir ActiveMat(j).name],[DormantDir ActiveMat(j).name]);
        end
    end
    
    if ~isempty(DormantMat) && SelectNum > 0
        movefile([DormantDir DormantMat(SelectNum).name],[ActiveDir DormantMat(SelectNum).name]);
    end
    
    ActiveMat = dir([ActiveDir '*.mat']);
end

%Load the database, depending on what the output argument was
if nargout == 1
    varargout{1} = [ActiveDir ActiveMat(1).name];
elseif nargout >= 3
    load([ActiveDir ActiveMat(1).name]);
    varargout{1} = Vmap;
    varargout{2} = Dmap;
    varargout{3} = Jmap;
    if nargout >= 4
        varargout{4} = Header;
        if nargout == 5
            varargout{5} = [ActiveDir ActiveMat(1).name];
        end
    end
end