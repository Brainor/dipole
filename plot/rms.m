function r = rms(v,varargin)
%RMS �󷽾���
% �ɼӲ���n��ʾ�Ե�nά����ƽ��
if ~isempty(varargin) && isscalar(varargin{:})
    n=varargin{1};
    r=sqrt(1/size(v,n)*sum(v.^2,n));
else
    r=sqrt(1/length(v)*sum(v.^2));
end
end

