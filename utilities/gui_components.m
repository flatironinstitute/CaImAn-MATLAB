function gui_components(rval_space,rval_time,max_pr,sizeA,A,Cn,options)
    

        figure;CC = plot_contours(A,Cn,options,0); title('Selected components','fontweight','bold','fontsize',14); close all;

        plt_cont = @(x) plot_contours(A,Cn,options,0,[],CC,[],x);

        d1 = options.d1;
        d2 = options.d2;

        figure;
            screensize = get(0,'Screensize' );
            fac = min(min((screensize(3:4)-200)./[2*d2,1.2*d1]),10);
            fig_size = round([100 100 fac*2*d2 fac*1.2*d1]);
            set(gcf, 'Position', fig_size);

            imagesc(Cn); axis equal; axis tight; axis off;
                set(gca,'Position',[0.01,0.1,0.5,0.8]);

        thr_space = @(x,y,z,w,v) find((rval_space >= x) & (rval_time >= y) & (sizeA >= z) & (sizeA < w) & (log(1-max_pr) < v));

        hsl1 = uicontrol('Style','slider','Min',0,'Max',1,...
                        'SliderStep',[0.01,0.01],'Value',options.space_thresh,'Tag','space_corr',...
                        'Position',[round(0.55*fig_size(3)) 50 20 round(0.85*fig_size(4))]);
        txt1 = uicontrol('Style','text','Position',[round(0.55*fig_size(3))-50 50+round(0.85*fig_size(4)) 100 50],'String','space threshold' );
        %ved1 = uicontrol('Style','edit','Position',[round(0.55*fig_size(3))-50 50+round(0.85*fig_size(4)) 100 50],'String',['val=',num2str(x1)] );

        hsl2 = uicontrol('Style','slider','Min',0,'Max',1,...
                        'SliderStep',[0.01,0.01],'Value',options.time_thresh,'Tag','time_corr',...
                        'Position',[round(0.6*fig_size(3)) 50 20 round(0.85*fig_size(4))]);            
        txt2 = uicontrol('Style','text','Position',[round(0.6*fig_size(3))-50 50+round(0.85*fig_size(4)) 100 50],'String','time threshold' );


        hsl3 = uicontrol('Style','slider','Min',0,'Max',120,...
                        'SliderStep',[0.01,0.01],'Value',options.min_size_thr,'Tag','min_size',...
                        'Position',[round(0.65*fig_size(3)) 50 20 round(0.85*fig_size(4))]);
        txt3 = uicontrol('Style','text','Position',[round(0.65*fig_size(3))-50 50+round(0.85*fig_size(4)) 100 50],'String','minimum size' );


        hsl4 = uicontrol('Style','slider','Min',0,'Max',500,...
                        'SliderStep',[.01,.05],'Value',options.max_size_thr,'Tag','max_size',...
                        'Position',[round(0.7*fig_size(3)) 50 20 round(0.85*fig_size(4))]);              

        txt4 = uicontrol('Style','text','Position',[round(0.7*fig_size(3))-50 50+round(0.85*fig_size(4)) 100 50],'String','maximum size' );

        hsl5 = uicontrol('Style','slider','Min',-10,'Max',0,...
                        'SliderStep',[0.01,0.01],'Value',log10(1-options.max_pr_thr),'Tag','max_pr',...
                        'Position',[round(0.75*fig_size(3)) 50 20 round(0.85*fig_size(4))]);              

        txt5 = uicontrol('Style','text','Position',[round(0.75*fig_size(3))-50 50+round(0.85*fig_size(4)) 100 50],'String','log(1-max_pr)' );

        pb = uicontrol('Style','pushbutton','Callback',@close_Callback);

        set(hsl1,'Callback',@(hObject,eventdata) plt_cont(thr_space(get(hObject,'Value'),hsl2.Value,hsl3.Value,hsl4.Value,hsl5.Value)))            
        set(hsl2,'Callback',@(hObject,eventdata) plt_cont(thr_space(hsl1.Value,get(hObject,'Value'),hsl3.Value,hsl4.Value,hsl5.Value)))   
        set(hsl3,'Callback',@(hObject,eventdata) plt_cont(thr_space(hsl1.Value,hsl2.Value,get(hObject,'Value'),hsl4.Value,hsl5.Value)))            
        set(hsl4,'Callback',@(hObject,eventdata) plt_cont(thr_space(hsl1.Value,hsl2.Value,hsl3.Value,get(hObject,'Value'),hsl5.Value)))   
        set(hsl5,'Callback',@(hObject,eventdata) plt_cont(thr_space(hsl1.Value,hsl2.Value,hsl3.Value,hsl4.Value,get(hObject,'Value'))))      

%     function close_Callback(hObject, eventdata)
% 
%             close all; 
%     end

end