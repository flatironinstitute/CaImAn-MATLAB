function varargout = ROI_GUI(varargin)
% ROI_GUI MATLAB code for ROI_GUI.template_fig
%      ROI_GUI, by itself, creates a new ROI_GUI or raises the existing
%      singleton*.
%
%      H = ROI_GUI returns the handle to a new ROI_GUI or the handle to
%      the existing singleton*.
%
%      ROI_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROI_GUI.M with the given input arguments.
%
%      ROI_GUI('Property','Value',...) creates a new ROI_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ROI_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ROI_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%             
%      varargin : 
%           ROIvars : for each ROI/neuron (see CNMFsetparams for further
%           documentation)
%               size
%               fitness
%               fitness delta
%               spatial corr
%               temporal corr
%               
%
%
%           A :
%           options : initial threshold values on the informtation of
%           ROIvars
%           template : 
%           CC : 
%           keep : 
%        
%          
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ROI_GUI

% Last Modified by GUIDE v2.5 31-Jul-2017 11:38:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ROI_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @ROI_GUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before ROI_GUI is made visible.
function ROI_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ROI_GUI (see VARARGIN)

% Choose default command line output for ROI_GUI
handles.output = hObject;

%% INITIATE VARIABLES

handles.sizetrain=15;
handles.maxcompsize=500;
handles.idx = false;
handles.selection = false;
handles.isimple =false;
hanldes.iscomputed =false;
handles.simplekeep = zeros(size(handles.sizetrain));
handles.genval = 0; 
handles.isdiscarded = false;
%from varargin
handles.A = varargin{1};
handles.options = varargin{2};
handles.template = varargin{3};
handles.CC = varargin{4};
handles.keep = varargin{5};
handles.ROIvars = varargin{6};

handles.cellnum = sum(handles.keep);

%we will define C graph and info of the neuron once a click has been made
set(handles.min_fit_delta,'String',num2str(handles.options.min_fitness_delta));
set(handles.min_fit,'String',num2str(handles.options.min_fitness));
set(handles.mincellsize,'String',num2str(handles.options.min_size_thr));
set(handles.maxcellsize,'String',num2str(handles.options.max_size_thr));
set(handles.rval_sp,'String',num2str(handles.options.space_thresh));
set(handles.rval_t,'String',num2str(handles.options.time_thresh));
set(handles.numcells,'String',num2str(handles.cellnum));
set(handles.General,'String',num2str(0));


set(handles.slider_min_fit_delta,'Value',handles.options.min_fitness_delta);
set(handles.slider_min_fit,'Value',handles.options.min_fitness);
set(handles.slider_min,'Value',handles.options.min_size_thr);
set(handles.slider_max,'Value',handles.options.max_size_thr);
set(handles.slider_rval_sp,'Value',handles.options.space_thresh);
set(handles.slider_rval_t,'Value',handles.options.time_thresh);
set(handles.SliderGeneral,'Value',0);

handles.keep = (handles.ROIvars.rval_space > handles.options.space_thresh) &...
(handles.ROIvars.rval_time > handles.options.time_thresh) & ... 
(handles.ROIvars.sizeA >= handles.options.min_size_thr) &...
(handles.ROIvars.sizeA <= handles.options.max_size_thr) & ...
(handles.ROIvars.fitness <= handles.options.min_fitness) & ...
(handles.ROIvars.fitness_delta <= handles.options.min_fitness_delta);
handles.A_keep = handles.A(:,handles.keep); 
handles.disc = zeros(size(handles.keep));
axes(handles.template_fig);
colormap gray

[~,~,im] = plot_contours(handles.A_keep,handles.template,handles.options,0,[],handles.CC,[],find(handles.keep));
set(im,'ButtonDownFcn',@(hObject,eventdata)ROI_GUI('template_fig_ButtonDownFcn',hObject,eventdata,guidata(hObject)));
%the image funtion will not fire until hit test is turned on
set(im,'HitTest','on');
axes(handles.traceplot);
plot(handles.ROIvars.C(1,:));
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = ROI_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
OUT = {handles.A_keep; handles.options; handles.keep};
handles.output = OUT;
varargout{1} = handles.output;

%% SPACE CORRELATION THRESHOLD
function slider_rval_sp_Callback(hObject, eventdata, handles)
if ~ handles.isimple
    handles.options.space_thresh = get(hObject,'Value');
    set(handles.rval_sp,'String',num2str(handles.options.space_thresh));
    set(handles.computing_text,'Visible','on');
    set(handles.numcells,'String',num2str(handles.cellnum));
    set(handles.ispace,'Value',0);
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end


%% TIME CORRELATION THRESHOLD
function slider_rval_t_Callback(hObject, eventdata, handles)
if ~ handles.isimple
    handles.options.time_thresh = get(hObject,'Value');
    set(handles.rval_t,'String',num2str(handles.options.time_thresh));
    set(handles.computing_text,'Visible','on');
    set(handles.numcells,'String',num2str(handles.cellnum));
    set(handles.istimecorr,'Value',0);
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end



%% MAXIMUM CELL SIZE
function slider_max_Callback(hObject, eventdata, handles)
if ~ handles.isimple
    handles.options.max_size_thr = round(get(hObject,'Value'));
    set(handles.maxcellsize,'String',num2str(handles.options.max_size_thr));
    set(handles.computing_text,'Visible','on');
    set(handles.numcells,'String',num2str(handles.cellnum));
    set(handles.ismax,'Value',0);
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end


%% MINIMUM CELL SIZE
function slider_min_Callback(hObject, eventdata, handles)
if ~ handles.isimple
    handles.options.min_size_thr = round(get(hObject,'Value'));
    set(handles.mincellsize,'String',num2str(handles.options.min_size_thr));
    set(handles.computing_text,'Visible','on');
    set(handles.numcells,'String',num2str(handles.cellnum));
    set(handles.ismin,'Value',0);
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end



%% MAX EVENT EXCEPTIONALITY
function slider_min_fit_Callback(hObject, eventdata, handles)
if ~ handles.isimple
    handles.options.min_fitness = get(hObject,'Value');
    set(handles.min_fit,'String',num2str(handles.options.min_fitness));

    set(handles.computing_text,'Visible','on');
    set(handles.numcells,'String',num2str(handles.cellnum));
    set(handles.isfit,'Value',0);
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end


%% MAX EVENT EXCEPTIONALITY delta
function slider_min_fit_delta_Callback(hObject, eventdata, handles)
if ~ handles.isimple
    handles.options.min_fitness_delta = get(hObject,'Value');
    set(handles.min_fit_delta,'String',num2str(handles.options.min_fitness_delta));

    set(handles.computing_text,'Visible','on');
    set(handles.numcells,'String',num2str(handles.cellnum));
    set(handles.isfitdelta,'Value',0);
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end

% --- Executes on slider movement.
function SliderGeneral_Callback(hObject, eventdata, handles)
if handles.isimple
    handles.genval = get(hObject,'Value');
    set(handles.min_fit,'String',num2str(handles.genval));

    set(handles.computing_text,'Visible','on'); 
    set(handles.numcells,'String',num2str(handles.cellnum));
    handles = replot(handles,eventdata); 
    guidata(hObject, handles);
end

%% check boxes
function ispace_Callback(hObject, eventdata, handles)
if get(handles.ispace,'Value')&& ~ handles.isimple
    handles.options.space_thresh = 0;
    set(handles.rval_sp,'String',num2str(handles.options.space_thresh));
    
    set(handles.computing_text,'Visible','on');
    
    set(handles.numcells,'String',num2str(handles.cellnum));
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end


function istimecorr_Callback(hObject, eventdata, handles)
if get(handles.istimecorr,'Value') && ~ handles.isimple
    handles.options.time_thresh = 0;
    set(handles.rval_t,'String',num2str(handles.options.time_thresh));
    
    set(handles.computing_text,'Visible','on');
    set(handles.numcells,'String',num2str(handles.cellnum));
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end


function ismax_Callback(hObject, eventdata, handles)
if get(handles.ismax,'Value')&& ~ handles.isimple
    handles.options.max_size_thr = 500;
    set(handles.maxcellsize,'String',num2str(handles.options.max_size_thr));
    set(handles.computing_text,'Visible','on');
    set(handles.numcells,'String',num2str(handles.cellnum));
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end


function ismin_Callback(hObject, eventdata, handles)
if get(handles.ismin,'Value')&& ~ handles.isimple
    handles.options.min_fitness_delta = 5;
    set(handles.mincellsize,'String',num2str(handles.options.min_size_thr));

    set(handles.computing_text,'Visible','on');
    set(handles.numcells,'String',num2str(handles.cellnum));
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end

function isfit_Callback(hObject, eventdata, handles)
if get(handles.isfit,'Value')&& ~ handles.isimple
    handles.options.min_fitness = 0;
    set(handles.min_fit,'String',num2str(handles.options.min_fitness));

    set(handles.computing_text,'Visible','on');
    set(handles.numcells,'String',num2str(handles.cellnum));
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
    
end

function isfitdelta_Callback(hObject, eventdata, handles)
if get(handles.isfitdelta,'Value') && ~ handles.isimple
    handles.options.min_fitness_delta = 0;
    set(handles.min_fit_delta,'String',num2str(handles.options.min_fitness_delta));

    set(handles.computing_text,'Visible','on');
    handles = replot(handles,eventdata);
    set(handles.numcells,'String',num2str(handles.cellnum));
    guidata(hObject, handles);
end


%% infer the weights and stuff
function simplemode_Callback(hObject, eventdata, handles)
if handles.iscomputed && ~ handles.isimple % we go direclty into the simple mode
    handles.isimple =true;

    set(handles.SliderGeneral,'Visible','on');
    set(handles.textGeneral,'Visible','on');
    set(handles.General,'Visible','on');
    set(handles.text53,'Visible','on');
    set(handles.text52,'Visible','on');

    set(handles.text23,'Visible','off');
    set(handles.text24,'Visible','off');
    set(handles.text25,'Visible','off');
    set(handles.text26,'Visible','off');
    set(handles.text36,'Visible','off');
    set(handles.test35,'Visible','off');
    set(handles.text37,'Visible','off');
    set(handles.text38,'Visible','off');
    set(handles.text39,'Visible','off');
    set(handles.text40,'Visible','off');
    set(handles.text15,'Visible','off');
    set(handles.text14,'Visible','off');
    set(handles.text16,'Visible','off');
    set(handles.text17,'Visible','off');
    set(handles.text19,'Visible','off');
    set(handles.text20,'Visible','off');

    set(handles.mincellsize,'Visible','off');
    set(handles.maxcellsize,'Visible','off');
    set(handles.min_fit,'Visible','off');
    set(handles.min_fit_delta,'Visible','off');
    set(handles.rval_t,'Visible','off');
    set(handles.rval_sp,'Visible','off');

    set(handles.Slider_min,'Visible','off');
    set(handles.Slider_max,'Visible','off');
    set(handles.Slider_min_fit,'Visible','off');
    set(handles.Slider_min_fit_delta,'Visible','off');
    set(handles.Slider_rval_t,'Visible','off');
    set(handles.Slider_rval_sp,'Visible','off');

    set(handles.ismin,'Visible','off');
    set(handles.ismax,'Visible','off');
    set(handles.isfit,'Visible','off');
    set(handles.isfitdelta,'Visible','off');
    set(handles.ispace,'Visible','off');
    set(handles.istimecorr,'Visible','off');
    
    set(handles.simplemode,'String','Regular Mode');
    handles.selection =false;   
    set(handles.badcomp_simple,'Visible','off');
    set(handles.goodcomp_simple,'Visible','off');
    set(handles.switchmode,'Visible','on');
    set(handles.find,'Visible','on');
    set(handles.remove_accept,'Visible','on');
else
    if ~ handles.iscomputed
        handles.selection=true;
        set(handles.badcomp_simple,'Visible','on');
        set(handles.goodcomp_simple,'Visible','on');
        set(handles.switchmode,'Visible','off');
        set(handles.find,'Visible','off');
        set(handles.remove_accept,'Visible','off');
    else 
        if handles.isimple
                handles.isimple =false;

                set(handles.SliderGeneral,'Visible','off');
                set(handles.textGeneral,'Visible','off');
                set(handles.General,'Visible','off');
                set(handles.text53,'Visible','off');
                set(handles.text52,'Visible','off');

                set(handles.text23,'Visible','on');
                set(handles.text24,'Visible','on');
                set(handles.text25,'Visible','on');
                set(handles.text26,'Visible','on');
                set(handles.text36,'Visible','on');
                set(handles.test35,'Visible','on');
                set(handles.text37,'Visible','on');
                set(handles.text38,'Visible','on');
                set(handles.text39,'Visible','on');
                set(handles.text40,'Visible','on');
                set(handles.text15,'Visible','on');
                set(handles.text14,'Visible','on');
                set(handles.text16,'Visible','on');
                set(handles.text17,'Visible','on');
                set(handles.text19,'Visible','on');
                set(handles.text20,'Visible','on');

                set(handles.mincellsize,'Visible','on');
                set(handles.maxcellsize,'Visible','on');
                set(handles.min_fit,'Visible','on');
                set(handles.min_fit_delta,'Visible','on');
                set(handles.rval_t,'Visible','on');
                set(handles.rval_sp,'Visible','on');

                set(handles.Slider_min,'Visible','on');
                set(handles.Slider_max,'Visible','on');
                set(handles.Slider_min_fit,'Visible','on');
                set(handles.Slider_min_fit_delta,'Visible','on');
                set(handles.Slider_rval_t,'Visible','on');
                set(handles.Slider_rval_sp,'Visible','on');

                set(handles.ismin,'Visible','on');
                set(handles.ismax,'Visible','on');
                set(handles.isfit,'Visible','on');
                set(handles.isfitdelta,'Visible','on');
                set(handles.ispace,'Visible','on');
                set(handles.istimecorr,'Visible','on');
                set(handles.simplemode,'String','Simple Mode');
        end 
    end
end
guidata(hObject, handles);


function badcomp_simple_Callback(hObject, eventdata, handles)
if handles.selection
    handles.keep(handles.i)=false;
    handles.i =+1;
    if(handles.i==handles.sizetrain)
        handles = compute(hObject);
    end
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end


function goodcomp_simple_Callback(hObject, eventdata, handles)
if handles.selection
    handles.keep(handles.i)=true;
    handles.i =+1;
    if(handles.i==handles.sizetrain)
        handles = compute(hObject,handles);
    end
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end


function Switchmode_Callback(hObject, eventdata, handles)
if ~ handles.selection
    if handles.isdiscarded 
        handles.isdiscarded =false;
        set(handles.Switchmode,'String','Accepted');
        set(handles.Switchmode,'ForegroundColor',[.0,1.0,.0]);
        set(handles.remove_accept,'String','BAD');
        set(handles.remove_accept,'ForegroundColor',[1.0,.0,.0]);
    else
        handles.isdiscarded =true;
        set(handles.Switchmode,'String','discarded');
        set(handles.Switchmode,'ForegroundColor',[1.0,.0,.0]);
        set(handles.remove_accept,'String','GOOD');
        set(handles.remove_accept,'ForegroundColor',[.0,1.0,.0]);
    end
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end


function find_Callback(hObject, eventdata, handles)
if ~ handles.selection
    disp('we will soon find it ')
%     axes(handles.template_fig);
%     [x,y,button]=ginput(1)
%     pixel=round([x y]);
%     if button==1
%         disp(['Adding pixel at:' num2str(fliplr(pixel))])
%         newcenters=[newcenters; fliplr(pixel)];
%         int_x = round(newcenters(end,1)) + (-sx:sx);
%         if int_x(1)<1
%             int_x = int_x + 1 - int_x(1);
%         end
%         if int_x(end)>options.d1
%             int_x = int_x - (int_x(end)-options.d1);
%         end
%         int_y = round(newcenters(end,2)) + (-sx:sx);
%         if int_y(1)<1
%             int_y = int_y + 1 - int_y(1);
%         end
%         if int_y(end)>options.d2
%             int_y = int_y - (int_y(end)-options.d2);
%         end
%         [INT_x,INT_y] = meshgrid(int_x,int_y);
%         coor = sub2ind([options.d1,options.d2],INT_x(:),INT_y(:));
%         Ypatch = reshape(Y(int_x,int_y,:),(2*sx+1)^2,T);                        
%         Y_res = Ypatch - A(coor,:)*C;
%         Y_res = bsxfun(@minus, Y_res, median(Y_res,2));
%         [atemp, ctemp, ~, ~, newcenter, ~] = greedyROI(reshape(Y_res,2*sx+1,2*sx+1,T), 1, options);
%         %[atemp, ctemp] = initialize_components(reshape(Y_res,2*sx+1,2*sx+1,T), 1,sx,options);  % initialize
%         % find contour
%         a_srt = sort(atemp,'descend');
%         ff = find(cumsum(a_srt.^2) >= cont_threshold*sum(a_srt.^2),1,'first');
%         K = K + 1;
%         A(coor,K) = atemp/norm(atemp);
%         C(K,:) = ctemp*norm(atemp);
%         new_center = com(A(:,end),options.d1,options.d2);
%         newcenters(end,:) = new_center;
%         scatter(new_center(2),new_center(1),'mo'); hold on; 
%         CC{K} = contour(reshape(A(:,end),options.d1,options.d2),[0,0]+a_srt(ff),'Linecolor',[1,0,1]/2);
%         CC{K}(CC{K}<1) = NaN;
%         CC{K}(:,CC{K}(1,:)>options.d2) = NaN;
%         CC{K}(:,CC{K}(2,:)>options.d1) = NaN;
%         hold on;
%         colormap(fig,cmap);
%         
%         drawnow;
%     end
end



%% click on the fig

function template_fig_ButtonDownFcn(hObject,eventdata,handles)
figc = handles.traceplot;
set(handles.computing_text,'Visible','on');
%coordinates = ginput(1)~
set(hObject.Parent,'Units','pixels');
%pos = get(hObject.Parent,'Position');
coor = get(hObject.Parent,'CurrentPoint');
coorx = coor(1,1);
coory = coor(1,2);
%width = pos(3);
%height = pos(4);
% we divide width and h with the one of the matrix we get the ratio and we
% use it to compare against the pos of the component
%hratio = height / handles.options.d1;
%wratio = width / handles.options.d2;
idx = false;
K = size(handles.A);
if not(handles.selection)
    for i = 1:K(2)
        if (handles.keep(i)==1 && not(handles.isdiscarded)) || (handles.keep(i)==0 && not(handles.isdiscarded))
            x = handles.ROIvars.cm(i,2);
            y = handles.ROIvars.cm(i,1);
            distx = abs(x - coorx);
            disty = abs(y - coory);
            if and(distx < handles.options.gSig(1)+1,disty < handles.options.gSig(1)+1)
               display('click on a comp')
               dist = disty*distx;
               if idx
                   if bdist > dist
                       idx = i;
                       bdist = dist;
                   end
               else 
                  idx = i;
                  bdist = dist; 
               end
            end
        end
    end
end
display(idx)
if idx
    set(handles.choosetext,'Visible','off');
    handles.traceplot = handles.ROIvars.C(idx,:);
    set(handles.num_disp,'String',num2str(idx));
    %set info neur    
    set(handles.timecorr,'String',num2str(handles.ROIvars.rval_space(idx)));
    set(handles.cellsize,'String',num2str(handles.ROIvars.sizeA(idx)));
    set(handles.spacecorr,'String',num2str(handles.ROIvars.rval_time(idx)));
    set(handles.fitness,'String',num2str(handles.ROIvars.fitness(idx)));
    axes(figc);
    plot(handles.ROIvars.C(idx,:));
    
else 
    set(handles.choosetext,'Visible','on');
    set(handles.num_disp,'String','none');
    set(handles.timecorr,'String','none');
    set(handles.cellsize,'String','none');
    set(handles.spacecorr,'String','none');
    set(handles.fitness,'String','none');
end
handles.idx = idx;
handles.traceplot = figc;
handles = replot(handles,eventdata);
set(handles.numcells,'String',num2str(handles.cellnum));
set(handles.computing_text,'Visible','off');
set(hObject.Parent,'Units','normalized');
guidata(hObject, handles);

%% discarding and accepting

function remove_accept_Callback(hObject, eventdata, handles)
if handles.idx && ~ handles.selection
    if handles.isdiscarded
        display('accepting');
        handles.disc(handles.idx)=0;
        handles = replot(handles,eventdata);
    else
        display('deletion');
        handles.disc(handles.idx)=1;
        handles = replot(handles,eventdata);
    end
    guidata(hObject, handles);
end

function [handles] = replot(handles,eventdata)
if not(handles.selection)
    if not(handles.isimple) % computing by comparing everything
        handles.keep = (handles.ROIvars.rval_space > handles.options.space_thresh) & (handles.ROIvars.rval_time > handles.options.time_thresh) & ... 
            (handles.ROIvars.sizeA >= handles.options.min_size_thr) & (handles.ROIvars.sizeA <= handles.options.max_size_thr) & ...
            (handles.ROIvars.fitness <= handles.options.min_fitness) & (handles.ROIvars.fitness_delta <= handles.options.min_fitness_delta) ...
            & not(handles.disc);
    else % how we compute the keep for the simple mode, we discard the one that the user discarded himself
        handles.keep = handles.simplekeep & not(handles.disc);
    end
    if handles.isdiscarded
       handles.A_keep = handles.A(:,not(handles.keep));
    else
       handles.A_keep = handles.A(:,handles.keep);
    end
    axes(handles.template_fig);
    [~,~,im] = plot_contours(handles.A_keep,handles.template,handles.options,0,[],handles.CC,[],find(handles.keep)); 
    disp('replot')
    set(im,'ButtonDownFcn',@(hObject,eventdata)ROI_GUI('template_fig_ButtonDownFcn',hObject,eventdata,guidata(hObject)));
    handles.cellnum = sum(handles.keep);
    set(handles.computing_text,'Visible','off');
else % the procedure of defining a test for the 
    handles.A_keep = handles.A(handles.i); 
    axes(handles.template_fig);
    plot_contours(handles.A_keep,handles.template,handles.options,0,[],handles.CC,[],handles.i); 
    handles.cellnum = 1;
    set(handles.numdisp,'String',num2str(handles.i));
    set(handles.numdisp,'String',num2str(handles.i));
    set(handles.computing_text,'Visible','off');
end


%% FINISHING PIPELINE FUNCTIONS

% --- Executes on button press in refine.
function refine_Callback(hObject, eventdata, handles)
if ~ handles.selection
end 

function extract_Callback(hObject, eventdata, handles)
if ~ handles.selection
end

function deconvolve_Callback(hObject, eventdata, handles)
if ~ handles.selection
end

function save_Callback(hObject, eventdata, handles)
if ~ handles.selection
    
end


%% PLOT





% %% EXPORT BUTTON
% % --- Executes on button press in save_button.
% function save_button_Callback(hObject, eventdata, handles)
% % hObject    handle to save_button (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% OUT = {handles.A_keep; handles.options; handles.keep};
% handles.output = OUT;
% guidata(hObject, handles);


%% CLOSE
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
display('exit')
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

%% CREATE FUNCTIONS
%%TODO : show the selected neuron

function template_fig_CreateFcn(~, ~, handles)

function ID_CreateFcn(~, ~, handles)

function session_CreateFcn(~, ~, handles)

function numcells_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor','white');

%space correlation
function slider_rval_sp_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function rval_sp_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor','white');

%time correlation
function slider_rval_t_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);


function rval_t_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor','white');

%max cell size
function slider_max_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function maxcellsize_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor','white');

%min cell size
function slider_min_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);


function mincellsize_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor','white');

%min fit
function slider_min_fit_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function min_fit_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor','white');

%min fit delta
function slider_min_fit_delta_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function min_fit_delta_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor','white');

%infos about neuron
function timecorr_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor','white');

function cellsize_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor','white');

function spacecorr_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor','white');

function fitness_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor','white');

function traceplot_CreateFcn(hObject, ~, handles)

%others
function SliderGeneral_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);


function figure1_ButtonDownFcn(hObject,~,handles)


function handles = compute_simple(handles)
%normalizing
MAXcellsize = 500; %% to be defined somwhere else
                                                           %REVIEW
set(handles.computing_text,'Visible','on');

val = handles.ROIvars;
K = size(handles.A);
for i = 1:K
    val(i).sizeA = val(i).sizeA/MAXcellsize;
    val(i).fitness = exp(val(i).fitness);
    val(i).fitness_delta = exp(val(i).fitness_delta);
end
infocomp = [val.rval_space,val.rval_time,val.sizeA,val.fitness,val.fitness_delta];
test = infocomp(handles.sizetrain+1:end);
train = infocomp(1:handles.sizetrain);
group = handles.keep(1:handles.sizetrain);
[C,err,~,~,coeff] = classify(test,train,group,'Quadratic');
% instanciating the weights.
handles.K = coeff(1,2).const;
handles.L = coeff(1,2).linear;
handles.Q = coeff(1,2).quadratic;

display(' error is : '+err)
handles.keep(handles.sizetrain+1:end)=C; % new keep values

% Function to compute K + L*v + v'*Q*v for multiple vectors
handles.f = @(v) K + v*L + sum((v*Q) .* v, 2);
%feval(f,v)
%displaying the simple mode view
handles.isimple =true;

set(handles.SliderGeneral,'Visible','on');
set(handles.textGeneral,'Visible','on');
set(handles.General,'Visible','on');
set(handles.text53,'Visible','on');
set(handles.text52,'Visible','on');

set(handles.text23,'Visible','off');
set(handles.text24,'Visible','off');
set(handles.text25,'Visible','off');
set(handles.text26,'Visible','off');
set(handles.text36,'Visible','off');
set(handles.test35,'Visible','off');
set(handles.text37,'Visible','off');
set(handles.text38,'Visible','off');
set(handles.text39,'Visible','off');
set(handles.text40,'Visible','off');
set(handles.text15,'Visible','off');
set(handles.text14,'Visible','off');
set(handles.text16,'Visible','off');
set(handles.text17,'Visible','off');
set(handles.text19,'Visible','off');
set(handles.text20,'Visible','off');

set(handles.mincellsize,'Visible','off');
set(handles.maxcellsize,'Visible','off');
set(handles.min_fit,'Visible','off');
set(handles.min_fit_delta,'Visible','off');
set(handles.rval_t,'Visible','off');
set(handles.rval_sp,'Visible','off');

set(handles.Slider_min,'Visible','off');
set(handles.Slider_max,'Visible','off');
set(handles.Slider_min_fit,'Visible','off');
set(handles.Slider_min_fit_delta,'Visible','off');
set(handles.Slider_rval_t,'Visible','off');
set(handles.Slider_rval_sp,'Visible','off');

set(handles.ismin,'Visible','off');
set(handles.ismax,'Visible','off');
set(handles.isfit,'Visible','off');
set(handles.isfitdelta,'Visible','off');
set(handles.ispace,'Visible','off');
set(handles.istimecorr,'Visible','off');

set(handles.simplemode,'String','Regular Mode');
handles.i = 1; %reset iteration for next time 
handles.selection =false;   
set(handles.badcomp_simple,'Visible','off');
set(handles.goodcomp_simple,'Visible','off');
set(handles.switchmode,'Visible','on');
set(handles.find,'Visible','on');
set(handles.remove_accept,'Visible','on');

set(handles.computing_text,'Visible','off');
guidata(handles.figure1, handles);


% function [model, llh] = logitBin(X, t, lambda)
% % Logistic regression for binary classification optimized by Newton-Raphson method.
% % Input:
% %   X: d x n data matrix
% %   t: 1 x n label (0/1)
% %   lambda: regularization parameter
% % Output:
% %   model: trained model structure
% %   llh: loglikelihood
% % Written by Mo Chen (sth4nth@gmail.com).
% if nargin < 3
%     lambda = 1e-2;
% end
% X = [X; ones(1,size(X,2))];
% [d,n] = size(X);
% tol = 1e-4;
% maxiter = 100;
% llh = -inf(1,maxiter);

% %     end
%     if incr < tol; break; end
% end
% llh = llh(2:iter);
% model.w = w;


