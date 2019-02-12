%addCDR3Prop will add CDR3 properties to the existing data to make it
%similar to Ilja' data sheet. The following properties will be added:
%
%  H-CDR3_Code
%  H-CDR3_MW
%  H-CDR3_VDWV
%  H-CDR3_Hydropathicity
%  H-CDR3_Aliphaticity
%  H-CDR3_Aromaticity
%  H-CDR3_pH7Charge
%  H-CDR3_pH6Charge
%  H-CDR3_pH5Charge
%
%  NOTE:
%  H-CDR3_pI is NOT included, since this doesn't make sense for CDR3 that
%  does NOT have a charge, or if there is only 1 charge (and thus cannot
%  have a net 0 charge for the CDR3 only).

function addCDR3Prop
FileNames = getBriliaFiles();
for f = 1:length(FileNames)
    tic
    [VDJdata, VDJheader] = openSeqData(FileNames{f});
    Map = getVDJmapper(VDJheader);

    AddHeader = {
            'H-CDR3_MW';
            'H-CDR3_VDWV';
            'H-CDR3_Code';
            'H-CDR3_Hydropathicity';
            'H-CDR3_Aromaticity';
            'H-CDR3_Aliphaticity';
            'H-CDR3_pH7Charge';
            'H-CDR3_pH6Charge';
            'H-CDR3_pH5Charge'};

    AddData = cell(size(VDJdata, 1), length(AddHeader));
    CDR3Seq = VDJdata(:, Map.hCDR3(1));
    AddData(:, contains(AddHeader, 'H-CDR3_Code')) = cellfun(@convAA2PropMEX, CDR3Seq, 'unif', false);
    AddData(:, contains(AddHeader, 'H-CDR3_MW'))   = num2cell(calcMW(CDR3Seq));
    AddData(:, contains(AddHeader, 'H-CDR3_VDWV')) = num2cell(calcVDWV(CDR3Seq));
    AddData(:, contains(AddHeader, 'H-CDR3_Hydropathicity')) = num2cell(calcHPI(CDR3Seq));
    AddData(:, contains(AddHeader, 'H-CDR3_Aromaticity'))    = num2cell(calcAromaticity(CDR3Seq));
    AddData(:, contains(AddHeader, 'H-CDR3_Aliphaticity'))   = num2cell(calcAliphaticity(CDR3Seq));
    AddData(:, contains(AddHeader, 'H-CDR3_pH7Charge')) = num2cell(calcCharge(CDR3Seq, 'pH', 7, 'IncludeTerminus', false));
    AddData(:, contains(AddHeader, 'H-CDR3_pH6Charge')) = num2cell(calcCharge(CDR3Seq, 'pH', 6, 'IncludeTerminus', false));
    AddData(:, contains(AddHeader, 'H-CDR3_pH5Charge')) = num2cell(calcCharge(CDR3Seq, 'pH', 5, 'IncludeTerminus', false));

    %Overwrite and/or add the new data
    [~, DelThese] = intersect(VDJheader, AddHeader);
    KeepLoc = ones(1, size(VDJdata, 2), 'logical');
    KeepLoc(DelThese) = 0;
    saveSeqData(FileNames{f}, [VDJdata(:, KeepLoc) AddData], [VDJheader(KeepLoc) AddHeader']);    

    [f length(FileNames) toc]
end