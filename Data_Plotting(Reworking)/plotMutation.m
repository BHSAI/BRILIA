%analyzeMutation will take the VDJdata file, look at the nT mutations, and
%plot the frequence of what NTs are added/deleeted.

function plotMutation(varargin)
FilePre = '';
if nargin == 4
    Vcmpr = varargin{1};
    Dcmpr = varargin{2};
    Jcmpr = varargin{3};
    VDJcmpr = varargin{4};
elseif nargin == 1
    Vcmpr = varargin{1}(:,1:2);
    Dcmpr = varargin{1}(:,3:4);
    Jcmpr = varargin{1}(:,5:6);
    VDJcmpr = varargin{1}(:,7:8);    
elseif nargin == 2 || nargin == 3
    load(varargin{1})
    V1 = Vmatrix;
    D1 = Dmatrix;
    J1 = Jmatrix;
    VDJ1 = VDJmatrix;
    
    load(varargin{2})
    V2 = Vmatrix;
    D2 = Dmatrix;
    J2 = Jmatrix;
    VDJ2 = VDJmatrix;
    
    Vmatrix = V2-V1;
    Dmatrix = D2-D1;
    Jmatrix = J2-J1;
    VDJmatrix = VDJ2-VDJ1;
    
    if nargin == 3
        FilePre = varargin{3};
    end
else
    
    error('Wrong number of inputs');
end

%Create the alphabet to number, number to matrix map.
M = zeros(1,200);
M(double('AGCT')) = 1:4;
AxisLabel = {'A' 'G' 'C' 'T'};

if nargin ~= 2 && nargin ~= 3
    %Frequency matrix fo ACGT going to ACGT. Rows are A,C,G,T from the query,
    %and Cols are A,C,G,T from the origin source. So, row 1 col 2 will tell you
    %how often does C mutate to A.

    Vmatrix = zeros(4,4);
    Vmatrix2 = zeros(4,4);
    for j = 1:size(Vcmpr,1);
        Results = double(upper(Vcmpr{j,1}));
        for k = 1:size(Results,2)
            Vmatrix(M(Results(1,k)),M(Results(2,k))) = Vmatrix(M(Results(1,k)),M(Results(2,k))) + 1;
        end
        Results2 = double(upper(Vcmpr{j,2}));
        for k = 1:size(Results2,2)
            Vmatrix2(M(Results2(1,k)),M(Results2(2,k))) = Vmatrix2(M(Results2(1,k)),M(Results2(2,k))) + 1;
        end
    end
    Vmatrix = Vmatrix + Vmatrix2;
    Vmatrix = Vmatrix./(repmat(sum(Vmatrix,1),4,1));

    Dmatrix = zeros(4,4);
    Dmatrix2 = zeros(4,4);
    for j = 1:size(Dcmpr,1);
        Results = double(upper(Dcmpr{j,1}));
        for k = 1:size(Results,2)
            if M(Results(1,k)) == 0 || M(Results(2,k)) == 0; continue; end            
            Dmatrix(M(Results(1,k)),M(Results(2,k))) = Dmatrix(M(Results(1,k)),M(Results(2,k))) + 1;
        end
        Results2 = double(upper(Dcmpr{j,2}));
        for k = 1:size(Results2,2)
            Dmatrix2(M(Results2(1,k)),M(Results2(2,k))) = Dmatrix2(M(Results2(1,k)),M(Results2(2,k))) + 1;
        end
    end
    Dmatrix = Dmatrix + Dmatrix2;
    Dmatrix = Dmatrix./(repmat(sum(Dmatrix,1),4,1));

    Jmatrix = zeros(4,4);
    Jmatrix2 = zeros(4,4);
    for j = 1:size(Jcmpr,1);
        Results = double(upper(Jcmpr{j,1}));
        for k = 1:size(Results,2)
            if Results(1,k) == Results(2,k); pause; end
            Jmatrix(M(Results(1,k)),M(Results(2,k))) = Jmatrix(M(Results(1,k)),M(Results(2,k))) + 1;
        end
        Results2 = double(upper(Jcmpr{j,2}));
        for k = 1:size(Results2,2)
            Jmatrix2(M(Results2(1,k)),M(Results2(2,k))) = Jmatrix2(M(Results2(1,k)),M(Results2(2,k))) + 1;
        end
    end
    Jmatrix = Jmatrix + Jmatrix2;
    Jmatrix = Jmatrix./(repmat(sum(Jmatrix,1),4,1));


    VDJmatrix = zeros(4,4);
    VDJmatrix2 = zeros(4,4);
    for j = 1:size(VDJcmpr,1);
        Results = double(upper(VDJcmpr{j}));
        for k = 1:size(Results,2)
            if M(Results(1,k)) == 0 || M(Results(2,k)) == 0; continue; end
            VDJmatrix(M(Results(1,k)),M(Results(2,k))) = VDJmatrix(M(Results(1,k)),M(Results(2,k))) + 1;
        end
        Results2 = double(upper(VDJcmpr{j,2}));
        for k = 1:size(Results2,2)
            VDJmatrix2(M(Results2(1,k)),M(Results2(2,k))) = VDJmatrix2(M(Results2(1,k)),M(Results2(2,k))) + 1;
        end
    end
    VDJmatrix = VDJmatrix + VDJmatrix2;
    VDJmatrix = VDJmatrix./(repmat(sum(VDJmatrix,1),4,1));
end

%Plotting results
if nargin == 2 || nargin == 3;
    Crange = [-0.05 0.05];
else
    Crange = [0 1];
end

figure(1)
subplot(2,2,1)
image(Vmatrix,'CDataMapping','scaled')
colormap(jet)
set(gca,'XTick',1:4)
set(gca,'XTickLabel',AxisLabel)
set(gca,'XAxisLocation','top')
set(gca,'YTick',1:4)
set(gca,'YTickLabel',AxisLabel)
set(gca,'OuterPosition',[0 0.5 0.5 0.5]);
set(gca,'Clim',Crange);
title('V gene mutations')
colorbar

subplot(2,2,2)
image(Dmatrix,'CDataMapping','scaled')
colormap(jet)
set(gca,'XTick',1:4)
set(gca,'XTickLabel',AxisLabel)
set(gca,'XAxisLocation','top')
set(gca,'YTick',1:4)
set(gca,'YTickLabel',AxisLabel)
set(gca,'OuterPosition',[0.5 0.5 0.5 0.5]);
set(gca,'Clim',Crange);
title('D gene mutations')
colorbar

subplot(2,2,3)
image(Jmatrix,'CDataMapping','scaled')
colormap(jet)
set(gca,'XTick',1:4)
set(gca,'XTickLabel',AxisLabel)
set(gca,'XAxisLocation','top')
set(gca,'YTick',1:4)
set(gca,'YTickLabel',AxisLabel)
set(gca,'OuterPosition',[0 0 0.5 0.5]);
set(gca,'Clim',Crange);
title('J gene mutations')
colorbar

subplot(2,2,4)
image(VDJmatrix,'CDataMapping','scaled')
colormap(jet)
set(gca,'XTick',1:4)
set(gca,'XTickLabel',AxisLabel)
set(gca,'XAxisLocation','top')
set(gca,'YTick',1:4)
set(gca,'YTickLabel',AxisLabel)
set(gca,'OuterPosition',[0.5 0 0.5 0.5]);
set(gca,'Clim',Crange);
title('VDJ gene mutations')
colorbar

set(gcf,'PaperPosition',[0 0 8 8])
[SaveFile, SavePath] = uiputfile('*.png','Save plot as ',FilePre);
saveas(gcf,[SavePath SaveFile]);
save([SavePath SaveFile '.mat'],'Vmatrix','Dmatrix','Jmatrix','VDJmatrix');
