function [S] = ANNOTATION_GRAPHIC_GENERALIZED(x_word,fs,phrase_info,phrase_info_num,n3sgram,folder, id)


%% n3sgram to LSF
% P = 40; % order of LSF P=2D+1 D is the number of formant
c = [];
for kk=1:size(n3sgram,2),
    s = n3sgram(:,kk);
    ceps = abs(ifft(log(abs((s)))));
    c(:,kk) = ceps(1:11);
end       
         
disp('Nn3sgram -> mcc done!')
n0 = 20;
indx = buffer(1:size(n3sgram,2),2*n0+1,2*n0,'nodelay');
clear D ma;
for kk =1:size(indx,2),
    
    a = sum(c(:,indx(:,kk)) * (-n0:n0).')  ./sum((-n0:n0).^2);
    D(kk) = sum(a.^2);
    
    
end
         
xy = [];
%% Locations of event target

S = D;

[~, localmaxima] = findpeaks(S);
[~, localminima] = findpeaks(-S);
D_D = diff(D);
SS = D_D.^2;
[~,localmaximaSS]= findpeaks(SS);
[~,localminimaSS]= findpeaks(-SS);

optima= unique(sort([1  localmaxima localminima localminimaSS+1  localmaximaSS+1  ]));
eliminateSet = [];
% for ii=2:length(phrase_info_num(1,:))
%     eliminateSet = [eliminateSet fix(phrase_info_num(1,ii))-4:fix(phrase_info_num(1,ii))+4];
% end
optima = setdiff(optima,[1 eliminateSet]);

optima= unique(sort([optima length(S)]));




hf = figure;
ha = axes('position',[0.1 0.3 0.8 0.6]);

plot((0:length(x_word)-1)/fs*1000,x_word);
hold on;
S_norm =S/max(S)*max(x_word);
[~,locs]=findpeaks(S_norm);
            hold on;
plot(S_norm);
hp=plot(locs,S_norm(locs),'v','markersize',5,'markerfacecolor','r');

hold on;
set(hp,'color',[0.860000 0.440000 0.580000])
%     plot(localmaximaSS+1, S(localmaximaSS+1),'dr');
% plot(localminimaSS+1, S(localminimaSS+1),'sb');

plot([phrase_info_num(1,:) size(n3sgram,2)],0 ,'ok');
for kk=1:length(phrase_info_num(1,:)),
    
    
    h1 = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,-0.01 ,phrase_info(1,kk));
    set(h1,'fontsize',16);
    set(h1,'color',[0.680000 0.090000 0.130000])
end
set(hp,'hittest','off')
set(gca,'fontsize',16);





hstart = uicontrol('style','pushbutton','string','Start',...
    'units','normalized','position',[0.2 0.1 0.2 0.1],...
    'callback',@startgin);
hstop = uicontrol('style','pushbutton','string','Done',...
    'units','normalized','position',[0.6 0.1 0.2 0.1],...
    'callback',@stopgin,'enable','off');
    function startgin(hObj,handles,eventdat)
        set(hObj,'Enable','off')
        set(hstop,'enable','on')
        set(hf,'WindowButtonMotionFcn',@changepointer)
        set(ha,'ButtonDownFcn',@getpoints)
    end
    function stopgin(hObj,handles,eventdat)
        set(hObj,'Enable','off')
        set(hstart,'enable','on')
        set(hf,'Pointer','arrow')
        set(hf,'WindowButtonMotionFcn',[])
        set(ha,'ButtonDownFcn',@getpoints)
        xy = getappdata(hf,'xypoints');
        line(xy(:,1),xy(:,2))
        
        load([folder, '/' id '.mat'])

        phrase_info_num_part(1,:)=[0;sort(unique(fix(xy(:,1))))];
        phrase_info_num_part(2,:)= [phrase_info_num_part(1,2:end)-phrase_info_num_part(1,1:end-1) size(n3sgram,2) - phrase_info_num_part(1,end)];

        save([folder, '/' id '.mat'],'minima1','phrase_info_num_part')

    end
    function changepointer(hObj,handles,eventdat)
        axlim = get(ha,'Position');
        fglim = get(hf,'Position');
        x1 = axlim(1)*fglim(3) + fglim(1);
        x2 = (axlim(1)+axlim(3))*fglim(3) + fglim(1);
        y1 = axlim(2)*fglim(4) + fglim(2);
        y2 = (axlim(2)+axlim(4))*fglim(4) + fglim(2);
        pntr = get(0,'PointerLocation');
        if pntr(1)>x1 && pntr(1)<x2 && pntr(2)>y1 && pntr(2)<y2
            set(hf,'Pointer','crosshair')
        else
            set(hf,'Pointer','arrow')
        end
    end
    function getpoints(hObj,~,~)
        cp = get(hObj,'CurrentPoint');
        line(cp(1,1),cp(1,2),'linestyle','none',...
            'marker','o','color','r')
        xy = getappdata(hf,'xypoints');
        xy = [xy;cp(1,1:2)];
        setappdata(hf,'xypoints',xy);
    end
end