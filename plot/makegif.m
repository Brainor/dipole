close all;
%��������
filename= 'plot.gif'; %���gif�ļ�������
gift1=0.5;gift=0.5;%t1:��һ��ͼ��ͣ��ʱ�䣬t����ʱ�������������Ʋ����ٶȣ���λ��
ext = {'\*.jpeg', '\*.jpg', '\*.png', '\*.pgm', '\*.tig', '\*.bmp'};  d = [ ];
for i = 1:length(ext)
    d =[d; dir( [cd,ext{i}] ) ]; % cd:��ǰ·��
end
str = {d.name};
if ~isempty(str)
    [Selection,ok] = listdlg('ListString',str,'name','Choose pictures','PromptString',...
        'Please choose pictures','SelectionMode','Multiple', 'ListSize',[400,200]);
else
    error('No picture find , add filename extension or change path.')
end
set(0,'defaultfigurecolor','w');

for i = 1:length(Selection)
    figure(i)
    imshow((imread(str{Selection(i)})),'InitialMagnification','fit','Border','tight')% Or :  d(Selection(i)).name == str{Selection}
    %      title(str(Selection(i)));
    text(0,0.95,num2str(i),'units','normalized');
    frame=getframe(i);
    im=frame2im(frame);%����gif�ļ���ͼ�������index����ͼ��
    [I,map]=rgb2ind(im,256);
    k=i-0;
    if k==1
        imwrite(I,map,filename,'gif','LoopCount',Inf,...
            'DelayTime',gift1);%loopcountֻ����i==1��ʱ�������
    else
        imwrite(I,map,filename,'gif','WriteMode','append',...
            'DelayTime',gift);%DelaylayTime��������gif�ļ��Ĳ��ſ���
    end
end
close all