function h=draw_pcolor(position,c,varargin)
    %DRAW_PCOLOR draw the time evolution pcolor image
    % position, cell, position parameter for subplot, e.g. {3,1,2}, 或者为axes
    % object
    % c is a cell, refer to C, X, Y in pcolor, currently is C(:,:,iz) if c{1} has 3 dimensions
    % varargin分别为Title,Xlabel,Ylabel
    if(iscell(position))
        h=subplot(position{:});
    else
        h=position;
    end
    X=c{2};Y=c{3};C=c{1};
    if(size(C,3)>1);C=squeeze(C(:,:,2));end
    ph=pcolor(h,X,Y,C);ph.ZData=ph.CData;
    axis(h,'tight');
    shading(h,'interp');
    colorbar(h);
    colormap(h,'jet');
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
end
