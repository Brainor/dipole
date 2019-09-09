function h=draw_plot(y,varargin)
%DRAW_PLOT 画二维图，带有标题和x,y标题, 其中x方向axis为tight
% y为plot的所有参量组成的cell
% varargin分别为Title,Xlabel,Ylabel

if ~iscell(y)
    error('Plot:incorrectType','Error.\nThe plot parameter should be in a cell.');
end
h=plot(y{:});
h=h.Parent;
x=h.Children.XData;
h.XLim=[min(x),max(x)];
h.Box='off';
if ~isempty(varargin)
    if(~isempty(varargin{1})),title(h,varargin{1},'Interpreter','latex');end
    if length(varargin)>1
        if(~isempty(varargin{2})),xlabel(h,['$$',varargin{2},'$$'],'Interpreter','latex');end
        if length(varargin)>2
            if(~isempty(varargin{3})),ylabel(h,['$$',varargin{3},'$$'],'Interpreter','latex');end
            if length(varargin)>3
                set(h,varargin{4:end});
            end
        end
    end
end
% outerpos = h.OuterPosition;
% ti = h.TightInset;
% h.Position = [outerpos(1) + ti(1),outerpos(2) + ti(2), outerpos(3) - ti(1) - ti(3), outerpos(4) - ti(2) - ti(4)];
% h.LooseInset=h.TightInset;
end