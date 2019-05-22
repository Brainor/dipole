function h=draw_pcolor(position,c,varargin)
%DRAW_PCOLOR draw the time evolution pcolor image
% position, cell, position parameter for subplot, e.g. [3,1,2]
% c is a cell, refer to C, X, Y in pcolor, currently is C(:,:,iz) if c{1} has 3 dimensions
% varargin·Ö±ğÎªTitle,Xlabel,Ylabel
position=num2cell(position);
h=subplot(position{:});
X=c{2};Y=c{3};C=c{1};
if(size(C,3)>1);C=squeeze(C(:,:,2));end
ph=pcolor(h,X,Y,C);ph.ZData=ph.CData;
axis(h,'tight');
shading(h,'interp');
colorbar(h);
colormap(h,'jet');
if ~isempty(varargin)
    title(varargin{1},'Interpreter','latex');
    if length(varargin)>1
        xlabel(h,['$$',varargin{2},'$$'],'Interpreter','latex');
        if length(varargin)>2
            ylabel(h,['$$',varargin{3},'$$'],'Interpreter','latex');
            if length(varargin)>3
                set(h,varargin{4:end});
            end
        end
    end
end
end
