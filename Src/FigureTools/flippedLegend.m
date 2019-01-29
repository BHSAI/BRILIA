%flippedLegend is a modifier for the Matlab's legend, in which the legends
%are drawn in flipper order to matched stacked bar plots. Simply use
%flippedLegend similar to legend.m

function varargout = flippedLegend(varargin)
[Lx1, Lx2] = legend(varargin{:});

Types = arrayfun(@(x) x.Type, Lx2, 'un', 0);
TextIdx = find(strcmpi(Types, 'text'));
HggIdx = find(strcmpi(Types, 'hggroup'));

%Swap the text positions
AllPos = flipud(vertcat(Lx2(TextIdx).Position));
for j = 1:numel(TextIdx)
    Lx2(TextIdx(j)).Position = AllPos(j, :);
end

%Swap the box positions
AllVert = cell(numel(HggIdx), 1);
for j = 1:numel(HggIdx)
    AllVert{j} = get(Lx2(HggIdx(j)).Children, 'Vertices');
end
AllVert = flipud(AllVert);
for j = 1:numel(HggIdx)
    set(Lx2(HggIdx(j)).Children, 'Vertices', AllVert{j});
end

varargout{1} = Lx1;
varargout{2} = Lx2;