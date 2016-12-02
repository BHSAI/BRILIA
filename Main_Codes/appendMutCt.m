%appendMutation will fill in the mismatch count per segment, under the
%SHM_V, _Ndv, _D, _Ndj, and _J columns.

function VDJdata = appendMutCt(varargin)
if isempty(varargin)
    [FileNames, FilePath] = uigetfile('*Fix.xlsx','Selection files','multiselect','on');
    if ischar(FileNames)
        FileNames = {FileNames};
    end

    for j = 1:length(FileNames)
        FileName = FileNames{j};
        [VDJdata, NewHeader, ~, ~] = openSeqData([FilePath FileName]);
        getHeaderVar
        
        [~, ~, ~, ~, VMDNJmutRef] = findMutationFreq(VDJdata,NewHeader,'group');

        if findHeader(NewHeader,'SHM_V') == 0 %This doesn't have columns, so add.
            NewHeader = [NewHeader 'SHM_V' 'SHM_Nvd' 'SHM_D' 'SHM_Ndj' 'SHM_J'];
            VDJdata = [VDJdata num2cell(VMDNJmutRef)];
        else
            MutCtLoc = findHeader(NewHeader,{'SHM_V' 'SHM_Nvd' 'SHM_D' 'SHM_Ndj' 'SHM_J'});
            VDJdata(:,MutCtLoc) = num2cell(VMDNJmutRef);
        end

        %Before saving to xlsx, convert columns with matrix values into char
        VDJdata = reformatAlignment(VDJdata,1);
        for q = 1:size(VDJdata,1)
            for w = 1:3
                VDJdata{q,FamNumLoc(w)} = mat2str(VDJdata{q,FamNumLoc(w)});
            end
        end

        %Save to excel or csv file, depending on OS
        DotLoc = find(FileName == '.');
        if ispc
            xlswrite([FilePath FileName(1:DotLoc(end)-1) '.xlsx'],cat(1,NewHeader,VDJdata));
        else
            writeDlmFile(cat(1,NewHeader,VDJdata),[FilePath FileName(1:DotLoc(end)-1) '.csv'],'\t');
        end
    end
else
    VDJdata = varargin{1};
    NewHeader = varargin{2};
    getHeaderVar

    [~, ~, ~, ~, VMDNJmutRef] = findMutationFreq(VDJdata,NewHeader,'group');

    if findHeader(NewHeader,'SHM_V') == 0 %This doesn't have columns, so add.
        NewHeader = [NewHeader 'SHM_V' 'SHM_Nvd' 'SHM_D' 'SHM_Ndj' 'SHM_J'];
        VDJdata = [VDJdata num2cell(VMDNJmutRef)];
    else
        MutCtLoc = findHeader(NewHeader,{'SHM_V' 'SHM_Nvd' 'SHM_D' 'SHM_Ndj' 'SHM_J'});
        VDJdata(:,MutCtLoc) = num2cell(VMDNJmutRef);
    end
end
