function makevideo(varargin)
%MAKEVIDEO Make MP4 files
%   从图片序列创建mp4文件
%   第一个变量图片(目录+)文件名listing, 默认为*.png, Windows文件名支持?*, UNIX支持?*[,]
%   第二个变量为视频文件名name, 默认为'video'
%   第三个变量为index, 默认为1:end
%   第四个变量为framerate, 默认为2
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
if isunix%是UNIX系统
    [~,listing]=unix(['ls -1 ',filenames]);
    listing=split(listing,newline);%为cell
    listing=listing(1:end-1);
else%是Windows环境
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
    frames(size(frames,1)+rem(size(frames,1),2),size(frames,2)+rem(size(frames,2),2),1)=0;%%把frame长宽变成偶数    
    writeVideo(aviobj,frames);
end
close(aviobj)
end

