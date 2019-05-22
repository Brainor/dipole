function makevideo(varargin)
%MAKEVIDEO Make MP4 files
%   ��ͼƬ���д���mp4�ļ�
%   ��һ������ͼƬ(Ŀ¼+)�ļ���listing, Ĭ��Ϊ*.png, Windows�ļ���֧��?*, UNIX֧��?*[,]
%   �ڶ�������Ϊ��Ƶ�ļ���name, Ĭ��Ϊ'video'
%   ����������Ϊindex, Ĭ��Ϊ1:end
%   ���ĸ�����Ϊframerate, Ĭ��Ϊ2
filenames='*.png';
output_name='video';
index=[];
framerate=2;
if nargin>0
    filenames=varargin{1};
    if nargin>1
        output_name=varargin{2};
        if nargin>2
            index=varargin{3};
            if nargin==4
                framerate=varargin{4};
            end
        end
    end
end
if isunix%��UNIXϵͳ
    [~,listing]=unix(['ls -1 ',filenames]);
    listing=split(listing,newline);%Ϊcell
    listing=listing(1:end-1);
else%��Windows����
    listing=dir(filenames);
    listing=strcat({listing.folder},'\',{listing.name});
end
if ~isempty(index)
    listing=listing(index);
end
aviobj=VideoWriter([output_name,'.mp4'],'MPEG-4');
aviobj.FrameRate=framerate;
open(aviobj);
for i=1:length(listing)
    frames=imread(listing{i});
    frames(size(frames,1)+rem(size(frames,1),2),size(frames,2)+rem(size(frames,2),2),1)=0;%%��frame������ż��    
    writeVideo(aviobj,frames);
end
close(aviobj)
end

