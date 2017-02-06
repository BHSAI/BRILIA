%plotMutCorrScatter will plot the X0 -> X1 mutation frequencies as
%determined per segment of the full VDJ gene. For instance, it will compare
%A -> G mutation in the V gene segment versus the D and J gene segments.
%This is used mainly to discern how well the annotation results yielded
%consisted SHMs accross the VDJ segments. Returns the slope and Pearson
%correlation value for a 12-point scatter plot, where each point is a
%unique X0 -> X1 mutation.
%
%  [Gx, Ax] = plotMutCorrScatter
%
%  [Gx, Ax] = plotMutCorrScatter(MutData)
%
%  [Gx, Ax] = plotMutCorrScatter(VDJdata,VDJheader)
%
%  [Gx, Ax] = plotMutCorrScatter(...,Param,Value,...)
% 
%  [Gx, Ax] = plotMutCorrScatter(...,P)
%
%  P = plotMutCorrScatter('getinput')
%
%  INPUT
%    MutData: structured data of pairwise mutations, returned by getMutData
%    VDJdata: MxN cell of all VDJ annotaiton data
%    VDJheader: 1xN cell of all header name for VDJdata
%    P: a structure containing P.(Param) = Value information. This input
%       option might be easier for adjusting plot parameters. You can also
%       use the direct param, cell pair input listed below. To make a
%       default P structure, use P = plotMutCorrScatter('GetInput');
%
%    Some plot features can be set using the following Param-Value pairs:
%      Param           Value       Description
%      --------------  ----------  --------------------------------------
%      PairWith        'parent'    Uses SHM data between child and parent
%                      'germline'  Uses SHM data between child and germline
%      AddLegend       'y' 'n'     Show legend
%      AddXtitle       'y' 'n'     Show x axis label
%      AddXlabel       'y' 'n'     Show x axis value
%      AddYtitle       'y' 'n'     Show y axis label
%      AddYlabel       'y' 'n'     Show y axis label
%      FontName        Arial       Font name of X and Y axes
%      FontSize        10          Font size of X and Y axes
%      FigHeight       inches      Height of the whole figure
%      FigWidth        inches      Width of the whole figure
%      LegendFontSize  12          Set the font size of the legend
%      MarkerSize      250         Set the size of the scatter plot dots
%
%    *If you want to save, use these param-value paris:
%      Parmaeter Name  Value       Description
%      --------------- ---------   ----------------------------------------
%      Save            'n','y'     To save or not to save. Default 'n'.
%      SaveName        String      The file name to save everything. The
%                                    file name will append the GrpNum to
%                                    it. EX: if SaveName is 'A', then the
%                                    files will be saved as A.Grp1.tif. If
%                                    empty and Save = 'y', will ask user to
%                                    select folder and file name.
%      Format          'tif','png' Image format type to save the figure as.
%                      'jpg','fig'   Default is 'tif'.
%      DPI             N           Dots per square inch. Default is 300.
%
%  OUTPUT
%    Gx: figure handle
%    Ax: axes handle
%
%  See also getMutData
    
function varargout = plotMutCorrScatter(varargin)
%The special case check for extracting input to plotTreeData
JustGettingInput = 0;
if nargin == 1 && ischar(varargin{1})
    if strcmpi(varargin{1},'getinput')
        JustGettingInput = 1;
        varargin = {}; %Want to get defaults.
    end
end

if JustGettingInput == 0
    %See if user gave VDJdata and NewHeader, or MutData, or nothing
    if isempty(varargin) %Need to find file
        MutData = getMutationData;
    elseif ~isempty(varargin)
        if isempty(varargin{1})
            MutData = getMutationData;
            varargin(1) = [];
        elseif ischar(varargin{1}) %Filename was given
            MutData = getMutationData(varargin{1});
            varargin(1) = [];
        elseif length(varargin) == 2 && iscell(varargin{1}) && iscell(varargin{2}) %VDJdata and VDJheader was given
            MutData = getMutationData(varargin{1:2});
            varargin(1:2) = [];
        elseif isstruct(varargin{1})
            MutData = varargin{1};
            varargin(1) = [];
        end
    else
        error('getMutationData: Check the inputs');
    end
end

%Parse the inputs
P = inputParser;
addParameter(P,'PairWith','parent',@(x) any(validatestring(lower(x),{'parent','germline'})));
%Plotting parameters
addParameter(P,'AddLegend','y',@(x) any(validatestring(lower(x),{'y','n'})));
addParameter(P,'AddXlabel','y',@(x) any(validatestring(lower(x),{'y','n'})));
addParameter(P,'AddXtitle','y',@(x) any(validatestring(lower(x),{'y','n'})));
addParameter(P,'AddYlabel','y',@(x) any(validatestring(lower(x),{'y','n'})));
addParameter(P,'AddYtitle','y',@(x) any(validatestring(lower(x),{'y','n'})));
addParameter(P,'FontName','Arial',@ischar);
addParameter(P,'FontSize',10,@isnumeric);
addParameter(P,'FigHeight',3.3,@isnumeric);
addParameter(P,'FigWidth',3.3,@isnumeric);
addParameter(P,'LegendFontSize',10,@isnumeric)
addParameter(P,'MarkerSize',250,@isnumeric);
%Saving parameters
addParameter(P,'Save','n',@(x) any(validatestring(lower(x),{'y','n'})));
addParameter(P,'SaveName','',@ischar);
addParameter(P,'Format','tif',@(x) any(validatestring(lower(x),{'tif','jpg','png','fig'})));
addParameter(P,'DPI',300,@isnumeric);
parse(P,varargin{:});
P = P.Results;

if JustGettingInput == 1
    varargout{1} = P;
    return;
end

%Extract the appropriate parameters
if strcmpi(P.PairWith,'parent')
    Vmat = MutData.VmutPar;
    Dmat = MutData.DmutPar;
    Jmat = MutData.JmutPar;
else
    Vmat = MutData.VmutGerm;
    Dmat = MutData.DmutGerm;
    Jmat = MutData.JmutGerm;
end

%Normalize the mutations per ACGT bp that was mutated (or per each col)
Vmat = Vmat./repmat(sum(Vmat,1),4,1);
DJmat = Dmat+Jmat;
DJmat = DJmat./repmat(sum(DJmat,1),4,1);

%Linearize and get rid of the 0 diagonal ones
Vx = Vmat(:);
Vx(1:5:end) = [];
DJy = DJmat(:);
DJy(1:5:end) = [];

%Get rid of the T's
Vx2 = Vx(1:end);
DJy2 = DJy(1:end);

%Plot the 12 dots according to the marker styles
MarkerStyles = {'d' 'o' 's' '^'}; %Germline shape
MarkerColors = {[0 0.6 0]; [1 0 0]; [.1 .1 1]; [.3 .3 .3]}; %Mutant color 
Germ2MutIdx = [...
                2 1;
                3 1;
                4 1;
                1 2;
                3 2;
                4 2;
                1 3;
                2 3;
                4 3;
                1 4;
                2 4;
                3 4]; %Left col is mutant, right side is germline.

%Set the figure size
Gx = figure;
set(Gx,'units','pixel');
set(Gx,'position',[50 50 600 600]);

%Setup the plot axes properties
Ax = gca;
xlim(Ax,[0 1]);
ylim(Ax,[0 1]);
set(Ax,'XTick',[0:0.2:1],'YTick',[0:0.2:1]);
set(Ax,'Box','on','TickDir','out','XGrid','on','Ygrid','on','FontSize',16) %Fix grid and other properties
set(Ax,'Position',[0.13 0.17 0.8 0.8]);

if strcmpi(P.AddXtitle,'y')
    HX = xlabel(Ax,'Mut. Freq. of V');
    set(HX,'FontSize',24,'FontName','Arial')
end

if strcmpi(P.AddYtitle,'y')
    HY = ylabel(Ax,'Mut. Freq. of DJ');
    set(HY,'FontSize',24,'FontName','Arial')
end

%Fix the label units
if strcmpi(P.AddXlabel,'y')
    Xticks = get(Ax,'XTickLabels');
    for k = 1:length(Xticks)
        Xticks{k} = sprintf('%0.1f',eval(Xticks{k}));
    end
    set(Ax,'XTickLabels',Xticks);
else
    set(Ax,'XTickLabels',[]);
end

if strcmpi(P.AddYlabel,'y')
    Yticks = get(Ax,'YTickLabels');
    for k = 1:length(Yticks)
        Yticks{k} = sprintf('%0.1f',eval(Yticks{k}));
    end
    set(Ax,'YTickLabels',Yticks);
else
    set(Ax,'YTickLabels',[]);
end

%Extend the plots to fill space
set(Ax,'OuterPosition',[0 0 1 1]);
TightInsets = get(Ax,'TightInset');
Positions = [TightInsets(1) TightInsets(2) 1-TightInsets(3)-TightInsets(1) 1-TightInsets(4)-TightInsets(2)];
set(Ax,'Position',Positions)

%Plot the data
for j = 12:-1:1
    hold(Ax,'on')
    LX = scatter(Vx(j),DJy(j),'fill');
    hold(Ax,'off')
    Mstyle = MarkerStyles{Germ2MutIdx(j,1)};
    Mcolor = MarkerColors{Germ2MutIdx(j,2)};
    set(LX,'Marker',Mstyle,'MarkerFaceColor',Mcolor,'SizeData',P.MarkerSize,'MarkerEdgeColor',[1 1 1],'LineWidth',0.3)
end

%Calculate Rsq
x = Vx2;
y = DJy2;
[Pval,~] = polyfit(x,y,1);
Slope = Pval(1); %Difference in slope
xval = [0:0.1:1];
yval = polyval(Pval,xval);
hold(Ax,'on')
plot(Ax,[0 1],[0 1],'--','color',[0.6 0.6 0.6])
plot(Ax,xval,yval,'-k')

%Calculate the R pearson correlation
R = corr(x,y);
text(0.95,0.2,sprintf('R_{corr} =  %0.2f',R),'FontName','Arial','FontSize',P.LegendFontSize,'HorizontalAlignment','right');
text(0.95,0.1,sprintf('Slope =  %0.2f',Slope),'FontName','Arial','FontSize',P.LegendFontSize,'HorizontalAlignment','right');

%Draw the legend
if strcmpi(P.AddLegend,'y')
    LegXcoor1 = [0.1 0.1 0.1 0.1];
    LegXcoor2 = LegXcoor1 + 0.15;
    LegYcoor = [0.9 0.83 0.76 0.69];
    LegText1 = {'A_{0}','C_{0}','G_{0}','T_{0}'};
    LegText2 = {'A_{1}','C_{1}','G_{1}','T_{1}'};
    hold(Ax,'on')
    for k = 1:4
        text(LegXcoor1(k)+0.03,LegYcoor(k),LegText1{k},'FontName','Arial','FontSize',16,'FontWeight','bold','VerticalAlignment','middle','Color',MarkerColors{k});
        scatter(Ax,LegXcoor2(k),LegYcoor(k),'SizeData',250,'Marker',MarkerStyles{k},'MarkerEdgeColor','k','LineWidth',1);
        text(LegXcoor2(k)+0.03,LegYcoor(k),LegText2{k},'FontName','Arial','FontSize',16,'VerticalAlignment','middle');
    end
    hold(Ax,'off')
end

%Save the figure, if the user wants
SaveNamePre = ''; %Prefix for the file name to save. Initialized with ''.
if strcmpi(P.Save,'y')
    if isempty(SaveNamePre) %Need to establish the prefix, save path, file ext here.
        if isempty(P.SaveName) %Need to ask users to select file name
            [SaveFile, SavePath] = uiputfile('*.tif;*.png;*.jpg;*.fig');
            [SavePath, SaveFile, SaveExt] = parseFileName([SavePath SaveFile]);
        else
            [SavePath, SaveFile, SaveExt] = parseFileName(P.SaveName);
            if isempty(SaveExt)
                SaveExt = ['.' P.Format]; %Use the format option specified in the input. P.Format cannot be empty.
                if SaveExt(2) == '.'; SaveExt(2) = []; end %Prevents ..jpg file format nuissances
            end
        end
        DotLoc = find(SaveFile == '.');
        if isempty(DotLoc)
            SaveNamePre = SaveFile;
        else
            SaveNamePre = SaveFile(1:DotLoc(end)-1);
        end
    end

    %Assemble the full save name, and save depending on file ext
    FullSaveName = sprintf('%s%s.d%s',SavePath,SaveNamePre,SaveExt);
    switch SaveExt
        case '.tif'
            print(Gx,FullSaveName,'-dtiff',['-r' num2str(P.DPI)]);
        case '.jpg'
            print(Gx,FullSaveName,'-djpeg',['-r' num2str(P.DPI)]);
        case '.png'
            print(Gx,FullSaveName,'-dpng',['-r' num2str(P.DPI)]);                
        case '.fig'
            saveas(Gx,FullSaveName);
    end
end  

if nargout >= 1
    varargout{1} = Gx;
    if nargout >= 2
        varargout{2} = Ax;
    end
end
