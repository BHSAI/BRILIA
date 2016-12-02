%extractVDJdata will extract specific sequences from VDJdata files.
%
% VDJdataS = extractClones(FileName,FilePath,Property,Value)
% VDJdataS = extractClones(VDJdata,NewHeader,Property,Value)
% VDJdataS = extractClones([],[],Property,Value) %Will ask user to open file.
% Property = 'GrpNum','SeqNum',

function VDJdataS = extractVDJdata(Input1,Input2,varargin)
if ischar(Input1) %It's probably the file name
    [VDJdata,NewHeader] = openSeqData([Input2 Input1]);
elseif isempty(Input1) %You need to select file
    [VDJdata,NewHeader] = openSeqData;
else
    VDJdata = Input1;
    NewHeader = Input2;
end
getHeaderVar;

%Here is the property-to-header conversion table. [PropInputName PropInputLoc IsNumber]
%IsNumber: 0-string, 1-single number, 2-matrix
PropTable = {'Seq'      SeqLoc          0;      
             'RefSeq'   RefSeqLoc       0;
             'SeqNum'   SeqNumLoc       1;
             'GrpNum'   GrpNumLoc       1; 
             'CDR3aa'   CDR3Loc(1)      0;
             'CDR3len'  CDR3Loc(2)      1;
             'CDR3prop' CDR3Loc(1)      0};

%Determine the property values
PropNames = varargin(1:2:end);
PropValue = varargin(2:2:end);
PropColLoc = findHeader(PropTable(:,1),PropNames);

%Make filter for search
KeepThis = ones(size(VDJdata,1),1) == 1;
for k = 1:length(PropNames)
    %Find col. Is it numeric or string?
    IsNumber = PropTable{PropColLoc(k),3};
    VDJColLoc = PropTable{PropColLoc(k),2};
    VDJValue = PropValue{k};
    
    if strcmpi(PropTable{PropColLoc(k),1},'CDR3prop')
        TempCDR = VDJdata(:,CDR3Loc(1)); %Store it temporarily
        VDJdata(:,CDR3Loc(1)) = getPropertyCode(VDJdata(:,CDR3Loc(1)));
    end
    
    if IsNumber == 1
        KeepThis2 = zeros(size(VDJdata,1),1,'logical');
        SearchMat = cell2mat(VDJdata(:,VDJColLoc));
        for c = 1:length(VDJValue)
            KeepThis2 = KeepThis2 | SearchMat == VDJValue(c);
        end
    elseif IsNumber == 0
        for j = 1:size(VDJdata,1)
           if ~isempty(regexpi(VDJdata{j,VDJColLoc},VDJValue))
               KeepThis2(j) = 1;
           end
        end
    end
    KeepThis = KeepThis & KeepThis2;
    
    if strcmpi(PropTable{PropColLoc(k),1},'CDR3prop')
        VDJdata(:,CDR3Loc(1)) = TempCDR;
    end
end

VDJdataS = VDJdata(KeepThis,:);

