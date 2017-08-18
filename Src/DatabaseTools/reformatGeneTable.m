%reformateGeneTable will fill in the empty spaces in the gene table file,
%caused by IMGT's merged table format issues. 
%
%  INPUT
%    The input files are excel tables that have undergone the following:
%    1) Copy and paste the IMGT gene table into Excel
%    2) Unmerge all cells
%    3) Move headers with "Strain", "Clone Names", etc to the top row, 
%       overwritting any headers above these. 
%    4) Delete extra header rows
%    5) Save the excel table as [filename]_Raw.xlsx
%
%  OUTPUT
%    Returns a ; delimited formatted table that removes empty columns, and
%    fills in the column with the full gene name, which is required for the
%    next step of extracting strains given a full gene name.
%
%  NOTE
%    Gene table can be found at IMGT, for instance at
%    http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&
%    repertoire=genetable&species=Mus_musculus&group=IGKV.
%
%  See also findStrainMatch

function reformatGeneTable
%Select the VDJ gene table file
[FileNames, FilePath] = uigetfile('*.xlsx', 'Find the Gene Table excel file', 'MultiSelect', 'on');
if ischar(FileNames)
    FileNames = {FileNames};
end

for f = 1:length(FileNames)
    FileName = FileNames{f};
    [~, ~, Vtable] = xlsread([FilePath FileName]);
    FullNameLoc = findCell(Vtable(1, :), {'IMGT full name', 'IMGT allele name'});

    %Remove NANs and then empty columns
    for r = 1:size(Vtable,1)
        for c = 1:size(Vtable,2)
            if isnan(Vtable{r,c})
                Vtable{r,c} = '';
            end
        end
    end
    
    DelCol = ones(1, size(Vtable, 2), 'logical');
    for c = 1:size(Vtable, 2)
        for r = 1:size(Vtable, 1)
            if ~isempty(Vtable{r, c})
                DelCol(c) = 0;
                break;
            end
        end
    end
    Vtable(:, DelCol) = [];

    %Fill in the full names only entries
    CurName = '';
    for j = 2:size(Vtable, 1)
        if isempty(CurName)
            if isempty(Vtable{j, FullNameLoc}); 
                continue
            else 
                CurName = Vtable{j, FullNameLoc};
            end
        else
            if isempty(Vtable{j, FullNameLoc})
                Vtable{j, FullNameLoc} = CurName;
            else
                CurName = Vtable{j, FullNameLoc};
            end
        end
    end

    %Save the new files
    DotLoc = find(FileName == '.');
    SaveName = [FileName(1:DotLoc(end)-4) 'Formatted.csv'];
    writeDlmFile(Vtable, [FilePath SaveName], ';');
end
