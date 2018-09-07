function [f1, ax] = myplot(rect, varargin)
if nargin == 1
    varargin = {rect};
    rect = [];
end
if isempty(rect)
    rect = [500 450];
end
fontName = 'Times New Roman';
f1 = figure('Position', [0, 0, rect]);
b = false;
for itr = 1:numel(varargin)
    if ischar(varargin{itr})
        if strcmpi('linewidth', varargin{itr})
            b =true;
        end
    end
end
if b
    plot(varargin{:})
else
    plot(varargin{:}, 'linewidth', 1);
end
ax = gca;
set(ax, 'FontName',fontName,'FontSize',20 );
grid on
box on

end