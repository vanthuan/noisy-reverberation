function TD_ANALYSIS_MCC_GRAPHIC_ANNOTATE_MINIMA(ndBs,id,jj,phrase_info,phrase_info_num,n3sgram, f0raw)

modified_events= [];
chosenEvents=[];
%% n3sgram to LSF
% P = 40; % order of LSF P=2D+1 D is the number of formant
c = [];
for kk=1:size(n3sgram,2),
    s = n3sgram(:,kk);
    ceps = abs(ifft(log(abs((s)))));
    c(:,kk) = ceps(1:13);
end       
         
disp('Nn3sgram -> mcc done!')
n0 = 20;
indx = buffer(1:size(n3sgram,2),2*n0+1,2*n0);
clear D ma;
for kk =1:size(indx,2),
    indices = indx(:,kk);
    c_comp= zeros(13,length(indices));
    c_comp(:,find(indices > 0)) = c(:,indices(indices > 0));
    a = sum(c_comp * (-n0:n0).')  ./sum((-n0:n0).^2);
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

[fBin,timebin] = size(n3sgram);
% waterfall(db(n3sgram))
% view(ha,[0.5 63.6]);
imagesc((1:timebin),(0:fBin-1)/1024*16000, db(n3sgram));
axis tight xy; 
colormap('jet')
ylim([0 5000]);

S = -200*log10(S) +500;
hold on; hp= plot(S,'b', 'linewidth',1.7);
hold on; plot((1:timebin),100*log10(f0raw) ,'k','linewidth',1.5);
plot(fix([phrase_info_num(1,:)+1 size(n3sgram,2)]), S(fix([phrase_info_num(1,:) size(n3sgram,2)-1])+1),'ro','markersize',5)
plot([optima size(n3sgram,2)],[ S(optima) 0]+500,'m.','markersize',10);

plot([phrase_info_num(1,:)+1 size(n3sgram,2)],0 ,'ok');
for kk=1:length(phrase_info_num(1,:)),
       
    h1 = text(phrase_info_num(1,kk) +phrase_info_num(2,kk)/2 ,min(min(db(n3sgram(1,:)))) ,phrase_info(1,kk));
    set(h1,'fontsize',16);
    set(h1,'color',[0.680000 0.130000 0.130000])
end
title([ndBs{jj} 'dB'])                                            
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
        load([id '.mat'])
        minima = sort(unique([fix(xy(:,1)); size(n3sgram,2)]));
        modified_events = unique(modified_events);
        save([id '.mat'],'minima','phrase_info_num' ,'modified_events')

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
       
        
        switch get(hf,'SelectionType')
            case 'normal'
                disp('single click')
                line(cp(1,1),cp(1,2),'linestyle','none',...
                    'marker','o','color','k','markersize',5)
                xy = getappdata(hf,'xypoints');
                xy = [xy;cp(1,1:2)];
            case 'open'            
                disp('double click')

                
            case 'alt'
                disp('right click')
                chosenEvent = fix(cp(1,1));
                modified_events =[ modified_events chosenEvent];
                line(cp(1,1),cp(1,2),'linestyle','none',...
                    'marker','o','color','k','markersize',7, 'markerfacecolor','k');
                save([id '_chosen_event.mat'],'chosenEvent');
                chosenEvents =[chosenEvents chosenEvent];
                save([id '_chosen_events.mat'],'chosenEvents')

                xy = getappdata(hf,'xypoints');
                xy = [xy;cp(1,1:2)];

        end
        
        
        setappdata(hf,'xypoints',xy);
    end
end