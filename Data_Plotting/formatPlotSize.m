%formatPlotSize will set the plots sizes will create a template figure with
%specific sizes.

%  formatPlotSize(Ax,FigW,FigH) Ax = plot axes handle, FigW is the desired
%  width in inches, and FigH is the desired height in inches.

function [Gx,Ax] = formatPlotSize(varargin)
%Set the inputs
P = inputParser;
addOptional(P,'Ax',[]);
addOptional(P,'FigW',0,@isnumeric);
addOptional(P,'FigH',0,@isnumeric);
addOptional(P,'Option','outer',@ischar);
addParameter(P,'Unit','inch',@ischar);
addParameter(P,'AddBorder',0,@isnumeric);

parse(P,varargin{:})
Ax = P.Results.Ax;
FigW = P.Results.FigW;
FigH = P.Results.FigH;
Option = P.Results.Option;
Unit = P.Results.Unit;
AddBorder = P.Results.AddBorder; 

%1) Set or get the figure and axis handles
if isempty(Ax)
    Gx = figure();
    Ax = gca;
elseif ishandle(Ax)
    Gx = get(Ax,'parent');
end

%Set everything to the specified units, no normalize.
set(Ax,'Units',Unit);
set(Gx,'Units',Unit);
set(Gx,'PaperUnits',Unit);

if strcmpi(Option,'inner')
    %If you have zerso for FigW or FigH, set it to current inner size;
    if FigW == 0 || FigH == 0
        CurPos = get(Ax,'position');
        if FigW == 0 
            FigW = CurPos(3);
        end
        if FigH == 0
            FigH = CurPos(4);
        end
    end
    
    %You will readjust FigW and FigH such that the inner gca size matches
    %it. FigWi and FigHi are the inner dimension.
    FigWi = FigW;
    FigHi = FigH;

    %Extract the current positions
    TightPos = get(Ax,'TightInset');
    
    %If there is a colorbar,  you need to append TightPos.
    Cx = findall(Gx,'Type','colorbar');
    if ~isempty(Cx)
        Cx.Units = 'inch';        
        CxPos = Cx.Position;
        %Need to figure out how much text is there
        CxText = Cx.TickLabels;
        MaxTextLen = 0;
        MaxText = '';
        for k = 1:length(CxText)
            if length(CxText{k}) > MaxTextLen
                MaxTextLen = length(CxText{k});
                MaxText = CxText{k};
            end
        end
        %Figure out how big a text font is.
        Tx = text(0,0,MaxText,'FontName',Cx.FontName,'FontSize',Cx.FontSize,'units','inch');
        CxPos(3) = CxPos(3) + Tx.Extent(3);
        delete(Tx)

        CurPos = Ax.Position;

        TightT = zeros(1,4);
        TightT(1) = CurPos(1) - CxPos(1);
        TightT(2) = CurPos(2) - CxPos(2);        
        TightT(3) = sum(CxPos([1 3])) - sum(CurPos([1 3]));
        TightT(4) = sum(CxPos([2 4])) - sum(CurPos([2 4]));
        
        TightPos = max([TightPos; TightT],[],1);        
    end

    %Add the border if there is any
    TightPos = TightPos + ones(1,4)*AddBorder; %Thickness of left, bot, right, top border
    
    
    NewAxPos = [TightPos(1) TightPos(2) FigWi FigHi];
    FigW = FigWi + sum(TightPos([1 3]));
    FigH = FigHi + sum(TightPos([2 4]));
    set(Ax,'Position',NewAxPos);

    %Now determine the new plot positions based on tight inset
    GxPos = [0 0 FigW FigH];
    %reset the paper size if it's off
    PaperSize = get(Gx,'PaperSize');
    if FigW+GxPos(1) > PaperSize(1)
        PaperSize(1) = FigW + GxPos(1);
    end
    if FigH+GxPos(2) > PaperSize(2)
        PaperSize(2) = FigH + GxPos(2);
    end
    set(Gx,'PaperSize',PaperSize);

    %Reset the figur position on the paper    
    set(Gx,'Position',GxPos)
    set(Gx,'PaperPosition',GxPos)
else    
    %If you have zerso for FigW or FigH, set it to current inner size;
    if FigW == 0 || FigH == 0
        CurPos = get(Ax,'position');
        TightPos = get(Ax,'TightInset') + ones(1,4)*AddBorder;
        if FigW == 0 
            FigW = CurPos(3) + sum(TightPos([1 3]));
        end
        if FigH == 0
            FigH = CurPos(4) + sum(TightPos([2 4]));
        end
    end
    
    %Now determine the new plot positions based on tight inset
    GxPos = [0 0 FigW FigH];
    %reset the paper size if it's off
    PaperSize = get(Gx,'PaperSize');
    if FigW+GxPos(1) > PaperSize(1)
        PaperSize(1) = FigW + GxPos(1);
    end
    if FigH+GxPos(2) > PaperSize(2)
        PaperSize(2) = FigH + GxPos(2);
    end
    set(Gx,'PaperSize',PaperSize);

    %Reset the figur position on the paper    
    set(Gx,'Position',GxPos)
    set(Gx,'PaperPosition',GxPos)

    %Final step is to rescale the inner gca
    TightPos = get(Ax,'TightInset'); %Thickness of left, bot, right, top border
    
    %If there is a colorbar,  you need to append TightPos.
    Cx = findall(Gx,'Type','colorbar');
    if ~isempty(Cx)
        Cx.Units = 'inch';
        CxPos = Cx.Position;
        
        %Need to figure out how much text is there
        CxText = Cx.TickLabels;
        MaxTextLen = 0;
        MaxText = '';
        for k = 1:length(CxText)
            if length(CxText{k}) > MaxTextLen
                MaxTextLen = length(CxText{k});
                MaxText = CxText{k};
            end
        end
        %Figure out how big a text font is.
        Tx = text(0,0,MaxText,'FontName',Cx.FontName,'FontSize',Cx.FontSize,'units','inch');
        CxPos(3) = CxPos(3) + Tx.Extent(3);
        delete(Tx)
        
        CurPos = Ax.Position;

        TightT = zeros(1,4);
        TightT(1) = CurPos(1) - CxPos(1);
        TightT(2) = CurPos(2) - CxPos(2);        
        TightT(3) = sum(CxPos([1 3])) - sum(CurPos([1 3]));
        TightT(4) = sum(CxPos([2 4])) - sum(CurPos([2 4]));
        
        TightPos = max([TightPos; TightT],[],1);        
    end
    
    NewAxPos = [TightPos(1) TightPos(2) GxPos(3)-sum(TightPos([1,3])) GxPos(4)-sum(TightPos([2,4]))];
    set(Ax,'Position',NewAxPos);      
end