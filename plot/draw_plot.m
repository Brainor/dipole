function h=draw_plot(y,varargin)
%DRAW_PLOT ����άͼ�����б����x,y����, ����x����axisΪtight
% yΪplot�����в�����ɵ�cell
% varargin�ֱ�ΪTitle,Xlabel,Ylabel

if ~iscell(y)
    error('Plot:incorrectType','Error.\nThe plot parameter should be in a cell.');
end
plot(y{:});
h=gca;
x=h.Children.XData;
h.XLim=[min(x),max(x)];
h.Box='off';
if ~isempty(varargin)
    title(varargin{1},'Interpreter','latex');
    if length(varargin)>1
        xlabel(['$$',varargin{2},'$$'],'Interpreter','latex');
        if length(varargin)>2
            ylabel(['$$',varargin{3},'$$'],'Interpreter','latex');
            if length(varargin)>3
                set(h,varargin{4:end});
            end
        end
    end
end

end