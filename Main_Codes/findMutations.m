%findMutations takes a VDJdata and returns the mutations observed when
%compared against the reference gene

%  S = findMutations(VDJdata)  where S is a Nx4 cell array.
%  [Vcmpr, Dcmpr, Jcmpr, VDJcmpr] = findMutations(VDJdata)
function varargout = findMutations(varargin)
if isempty(varargin)
    [FileName, FilePath] = uigetfile('*.xlsx;*.csv','Select VDJdata file');
    [~, ~, VDJdata] = xlsread([FilePath FileName]);
    VDJdata = filterHeader(VDJdata);
    VDJdata = reformatAlignment(VDJdata,3);
end

Vcmpr = cell(size(VDJdata,1),2);
Dcmpr = cell(size(VDJdata,1),2);
Jcmpr = cell(size(VDJdata,1),2);
VDJcmpr = cell(size(VDJdata,1),2);

for j = 1:size(VDJdata,1)
    %Extract information
    VMDNJ = cell2mat(VDJdata(j,4:8));
    Valign = VDJdata{j,14};
    Vref3del = VDJdata{j,13};
    Dalign = VDJdata{j,20};
    Dref5del = VDJdata{j,18};
    Jalign = VDJdata{j,26};
    Jref5del = VDJdata{j,24};
    
    %Find V classifier
    V0 = regexp(Valign(3,:),'\w');
    V0 = V0(end);
    ValignTrim = Valign(:,V0-VMDNJ(1)+1:V0-Vref3del);
    VmisMatch = ValignTrim(2,:) == ' ';
    %Might have to trim left side due to sequencing gap error
    q = 0;
    while VmisMatch(q+1) == 1
        q = q+1;
    end
    VmisMatch(1:q) = 0;
    Vcmpr{j,1} = [ValignTrim(1,VmisMatch); ValignTrim(3,VmisMatch)];
    Vmatch = VmisMatch == 0;
    Vmatch(1:q) = 0;
    Vcmpr{j,2} = [ValignTrim(1,Vmatch); ValignTrim(3,Vmatch)];
    
    %Find D classifier
    D0 = regexp(Dalign(3,:),'\w');
    D1 = D0(1);
    DalignTrim = Dalign(:,D1+Dref5del:D1+Dref5del-1+VMDNJ(3));
    DmisMatch = DalignTrim(2,:) == ' ';
    Dcmpr{j,1} = [DalignTrim(1,DmisMatch); DalignTrim(3,DmisMatch)];
    Dcmpr{j,2} = [DalignTrim(1,DmisMatch==0); DalignTrim(3,DmisMatch==0)];
   
    %Find J classifier
    J0 = regexp(Jalign(3,:),'\w','once');
    JalignTrim = Jalign(:,J0+Jref5del:J0+Jref5del-1+VMDNJ(5));
    JmisMatch = JalignTrim(2,:) == ' ';
    %Might have to trim right side due to sequencing gap error
    q = 0;
    while JmisMatch(end-q) == 1
        q = q+1;
    end
    JmisMatch(end-q:end) = 0;
    Jcmpr{j,1} = [JalignTrim(1,JmisMatch); JalignTrim(3,JmisMatch)];
    Jmatch = JmisMatch == 0;
    Jmatch(end-q:end) = 0;
    Jcmpr{j,2} = [JalignTrim(1,Jmatch); JalignTrim(3,Jmatch)];
    
    VDJcmpr{j,1} = [Vcmpr{j,1} Dcmpr{j,1} Jcmpr{j,1}];
    VDJcmpr{j,2} = [Vcmpr{j,2} Dcmpr{j,2} Jcmpr{j,2}];
end

if nargout == 1
    varargout{1} = cat(2,Vcmpr,Dcmpr,Jcmpr,VDJcmpr);
else
    varargout{1} = Vcmpr;
    varargout{2} = Dcmpr;
    varargout{3} = Jcmpr;
    varargout{4} = VDJcmpr;
end
    