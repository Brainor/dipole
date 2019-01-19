function r = rms(v,varargin)
%RMS 求方均根
% 可加参量n表示对第n维进行平均
if ~isempty(varargin) && isscalar(varargin{:})
    n=varargin{1};
    r=sqrt(1/size(v,n)*sum(v.^2,n));
else
    r=sqrt(1/length(v)*sum(v.^2));
end
end

