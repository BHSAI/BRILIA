%DEPRACATED: Use getAllHeaderVar to determine what kind of property.
%getDataType is a simple function get the type of data that should in a
%certain column of VDJdata.
%
%  DataType = getDataType(HeaderName)
%
%  [DataType, DataTypeTable] = getDataType(HeaderName)
%
%  INPUT
%    HeaderName: name of the VDJheader, or VDJheader{n}.
%
%  OUTPUT
%    DataType: the type of data stored for this HeaderName, which are
%      alphanum: char + num. Default if you can't find a match.
%      num: number or matrix
%      prop: amino acid property code
%      seq: nt or aa letters, but where X is wildcard
%    DataTypeTable: Nx2 cell matrix containing the following data:
%      Column 1: the name of the VDJheader element
%      Column 2: the type of data stored in the column, which are:

function varargout = getDataType(HeaderName)
DataTypeTable = {
    'SeqName'          'alphanum'   ;
    'SeqNum'           'num'        ;
    'GroupNum'         'num'        ;
    'CDR3_Property'    'prop'       ;
    'CDR3_AminoAcid'   'seq'        ;
    'CDR3_Length'      'num'        ;
    'Seq'              'seq'        ;
    'RefSeq'           'seq'        ;
    'TemplateCount'    'num'        ;
    'TreeChildCount'   'num'        ;
    'Functional'       'alphanum'   ;
    'V_GeneName'       'alphanum'   ;
    'V_MapNum'         'num'        ;
    'V_Deletion5'      'num'        ;
    'V_Deletion3'      'num'        ;
    'D_GeneName'       'alphanum'   ;
    'D_MapNum'         'num'        ;
    'D_Deletion5'      'num'        ;
    'D_Deletion3'      'num'        ;
    'J_GeneName'       'alphanum'   ;
    'J_MapNum'         'num'        ;
    'J_Deletion5'      'num'        ;
    'J_Deletion3'      'num'        ;
    'Length_V'         'num'        ;
    'Length_Nvd'       'num'        ;
    'Length_D'         'num'        ;
    'Length_Ndj'       'num'        ;
    'Length_J'         'num'        ;
    'SHM_V'            'num'        ;
    'SHM_Nvd'          'num'        ;
    'SHM_D'            'num'        ;
    'SHM_Ndj'          'num'        ;
    'SHM_J'            'num'        ;
    'TrimmedSeqLeft'   'seq'        ;
    'TrimmedSeqRight'  'seq'        ;
    'Nvd_DDins'        'seq'        ;
    'Ndj_DDins'        'seq'        ;
    'Misc'             'alphanum'   ;
    };

DataType = 'alphanum'; %Default value
DataTypeLoc = findCell(DataTypeTable(:,1),HeaderName);
if ~isempty(DataTypeLoc) && DataTypeLoc > 0
    DataType = DataTypeTable{DataTypeLoc(1),2};
end

if nargout >= 1
    varargout{1} = DataType;
    if nargout == 2
        varargout{2} = DataTypeTable;
    end
end
