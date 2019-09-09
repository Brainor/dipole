function h=delt(f,varargin)
    % DELT find the turbulent part. Minus the mean value in the y direction
    % ��ȥ�߽�Ĳ���
    % ��Ϊ�������������Ǳ߽�, �ڶ�������Ϊƽ����ά��
    h=zeros(size(f));
    size_expand=ones(1,ndims(f));
    if isempty(varargin)
        size_expand(2)=size(f,2)-2;
        switch(length(size(f)))
            case 3
                h(2:end-1,2:end-1,2:end-1)=f(2:end-1,2:end-1,2:end-1)-repmat(mean(f(2:end-1,2:end-1,2:end-1),2),size_expand);
            case 2
                h(2:end-1,2:end-1)=f(2:end-1,2:end-1)-repmat(mean(f(2:end-1,2:end-1),2),size_expand);
        end
    else
        size_expand(varargin{1})=size(f,varargin{1});
        h=f-repmat(mean(f,varargin{1}),size_expand);
    end
end