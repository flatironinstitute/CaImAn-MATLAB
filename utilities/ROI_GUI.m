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
%       
%       DOCUMENTATION :
%           this will create a User interface for the rest of the pipeline described in demoscript and 
%           demo pipeline. LOOK at demo pipeline to understand how to use this GUI and what to pass it
%           
%           It is a mix of everyhting found in those pipelines. In the future it will be also called by 
%           the initialization GUI 
%
%           run GUI will let you, refine the components manually, change the parameters and see the results,
%           analyse the trace of each components, use a classification algorithm to
%           find the components instead, add and remove components, save the ROIS
%           and finish the pipeline with button clicks. 
%           It is still an optional method.
%
% WARNING : do not use simplemode for now - minimum MATLAB VERSION 2014
%                   save is only compatible with the JSON matlabpackage
%      varargin : 
%           A : intial spatial components
%           options : initial threshold values on the informtation of
%           template : Cn
%           C : initial calcium traces
%           b : bakcground
%           f : background activity
%           Y :video matrix
%           P : parmaters found during the preprocessing
%           (OPTIONAL)
%           Coor coordingates of the contours
%           S : infered spikes
%
%
%   varargout :
%       TODO edit here
%       ROIvars.fitness
%              .fitness_delta
%                   fitness values of each components
%              .cm  : center of mass
%              .sizeA
%                   size of each compoents
%              .rval_space
%                   spacial correlation 
%              .rval_time 
%                   time correlation
%               
%       A_keep
%           the kept spatial components
%       handles.options 
%           updated options
%       handles.Cdec
%           deconvolved calcium traces
%       handles.YrA
%           updated residuals
%       handles.CC
%           updated contours of each components
%       handles.b
%           updated spatial background matrix
%       handles.f
%           updated temporal background matrix
%       handles.S 
%           found spikes after deconvolution
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ROI_GUI

% Last Modified by GUIDE v2.5 02-Aug-2017 12:27:12

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

%User defined values (look at the demos to knwo more about those values)
handles.sizetrain=25; %number of training data that the user needs to check
handles.MAXCELLSIZE=500; % maximum size of a component
handles.df_percentile = 30; % the percentile filtering along 
handles.window = 1000; % windiwd to compute noise
handles.min_spikes =3; % find spikes resulting in transients above min_sp x noise level


%Values for the program
handles.idx = false;
handles.selection = false;
handles.isimple =false;
handles.iscomputed =false;
handles.isdiscarded = false;
handles.isextracted = false;
handles.isdeconvolved = false;
handles.isrefined = false;
handles.saving = false;
handles.hasnew =false; 
handles.genval = 0; 
handles.shownum  =0;
handles.simplekeep = zeros(size(handles.sizetrain));


%from varargin

handles.Y = varargin{1};
handles.A = varargin{2};
handles.P = varargin{3};
handles.options = varargin{4};
handles.template = varargin{5};
handles.C = varargin{6};
handles.b = varargin{7};
handles.f = varargin{8};

Yr = reshape(handles.Y,handles.options.d1*handles.options.d2,size(handles.C,2));

[handles.A,handles.b,handles.C] = update_spatial_components(...
    Yr,handles.C,handles.f,[handles.A,handles.b],handles.P,handles.options);

handles.P.p = 0;    % set AR temporarily to zero for speed
[handles.C,handles.f,handles.P,handles.S,handles.YrA] = update_temporal_components(...
    Yr,handles.A,handles.b,handles.C,handles.f,handles.P,handles.options);


% compute center of mass
[handles.ROIvars.rval_space,handles.ROIvars.rval_time,handles.ROIvars.max_pr,...
    handles.ROIvars.sizeA,handles.keep] = classify_components(handles.Y,handles.A,handles.C,handles.b,handles.f,handles.YrA,handles.options);
handles.ROIvars.C =handles.C;
handles.ROIvars.cm = com(handles.A,handles.options.d1,handles.options.d2);
[d1,d2,handles.T] = size(handles.Y);                                % dimensions of dataset
handles.d = d1*d2;
handles.Yr = reshape(handles.Y,handles.d,handles.T);
handles.cellnum = sum(handles.keep);
handles.newneur = zeros(handles.cellnum);


% compute time variation of the components ( needed in the GUI )
traces = prctfilt(handles.ROIvars.C+handles.YrA,8,1000,100);
handles.ROIvars.fitness = compute_event_exceptionality(traces,0);
handles.ROIvars.fitness_delta = compute_event_exceptionality(diff(traces,[],2),0);

% We set everything into the GUI from the values that have sent from the User and computed before
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
set(handles.unrefinedtext,'Visible','off');

%computing a first batch of kept component to display something (could have called replot instead)
handles.keep = (handles.ROIvars.rval_space > handles.options.space_thresh) &...
(handles.ROIvars.rval_time > handles.options.time_thresh) & ... 
(handles.ROIvars.sizeA >= handles.options.min_size_thr) &...
(handles.ROIvars.sizeA <= handles.options.max_size_thr) & ...
(handles.ROIvars.fitness <= handles.options.min_fitness) & ...
(handles.ROIvars.fitness_delta <= handles.options.min_fitness_delta);
handles.A_keep = handles.A(:,handles.keep); 
%the values  for what the User choose to define as discarded and accepted component (not computed)
handles.disc = zeros(size(handles.keep));
handles.accp = zeros(size(handles.keep));
%first plottings
axes(handles.template_fig);
colormap gray
[handles.CC,~,im] = plot_contours(handles.A_keep,handles.template,handles.options,0,[],[],[],find(handles.keep));
set(im,'ButtonDownFcn',@(hObject,eventdata)ROI_GUI('template_fig_ButtonDownFcn',hObject,eventdata,guidata(hObject)));
%the image funtion will not fire until hit test is turned on
set(im,'HitTest','on');
axes(handles.traceplot);
plot(handles.ROIvars.C(1,:));
guidata(hObject, handles);
uiwait()

% --- Outputs from this function are returned to the command line.
%TODO need more output from this pipeline
function varargout = ROI_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
OUT = {handles.A_keep; handles.options; handles.Cdec;handles.ROIvars;handles.YrA;handles.CC;handles.b; handles.f;handles.S };
handles.output = OUT;
varargout{1} = handles.output;


%%SLIDERS -- 
%only usable if not being in simple mode always need to send the handles to the main 
%hObject of the GUI

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


%% check boxes for each sliders. puts the value of the slider to 0
% works only of not clicked on already
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


%% Will transform the UI in a simple one using a 
%classification method to define the right components
function simplemode_Callback(hObject, eventdata, handles) %%todebug
if ~ handles.hasnew && ~ handles.isdeconvolved % scope of when it should be clickable
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
        set(handles.text35,'Visible','off');
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
        set(handles.Switchmode,'Visible','on');
        set(handles.find,'Visible','on');
        set(handles.remove_accept,'Visible','on');
    else
        if handles.selection
            handles.selection=false;
            set(handles.badcomp_simple,'Visible','off');
            set(handles.goodcomp_simple,'Visible','off');
            set(handles.Switchmode,'Visible','on');
            set(handles.find,'Visible','on');
            set(handles.remove_accept,'Visible','on');
            handles.i = 1;
        else
            if ~ handles.iscomputed % we are in selection mode for the training dataset
                handles.selection=true;
                set(handles.badcomp_simple,'Visible','on');
                set(handles.goodcomp_simple,'Visible','on');
                set(handles.Switchmode,'Visible','off');
                set(handles.find,'Visible','off');
                set(handles.remove_accept,'Visible','off');
                handles.ranset = randperm(handles.cellnum,handles.sizetrain);
                handles.i =1;
            else  
                if handles.isimple %% we go back to regular mode
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
                        set(handles.text35,'Visible','on');
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
    end
else 
        set(handles.unrefinedtext,'Visible','on');
end
handles = replot(handles,eventdata);
guidata(hObject, handles);


function badcomp_simple_Callback(hObject, eventdata, handles)
if handles.selection && ~ handles.isrefined
    handles.keep(handles.i)=false;
    handles.i =handles.i +1;
    if(handles.i==handles.sizetrain)
        handles = compute_simple(handles);
    end
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end


function goodcomp_simple_Callback(hObject, eventdata, handles)
if handles.selection && ~ handles.isrefined
    handles.keep(handles.i)=true;
    handles.i =handles.i +1;
    if(handles.i==handles.sizetrain)
        handles = compute_simple(handles);
    end
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end



function ordercomp_Callback(hObject, eventdata, handles) %%todebug
if ~ handles.selection 
    if ~ handles.hasnew
        [handles.A,handles.ROIvars.C,~,~] = order_ROIs(handles.A,handles.ROIvars.C);
        handles.shownum = 1;
        handles = replot(handles,eventdata);
        guidata(hObject,handles);
    else 
        set(handles.unrefinedtext,'Visible','on');
        guidata(hObject,handles);
    end
end

%%reusing the function from manually refine components.m with 
%littles tweaks to make it work here. 
function find_Callback(hObject, eventdata, handles) %%todebug
if ~ handles.selection && ~ handles.isrefined
    disp('we will soon find it ')
    axes(handles.template_fig);
    [x,y,button]=ginput(1)
    pixel=round([x y]);
    if button==1
        disp(['Adding pixel at:' num2str(fliplr(pixel))])
        handles.ROIvars.cm=[handles.ROIvars.cm; fliplr(pixel)];
        int_x = round(handles.ROIvars.cm(end,1)) + (-handles.options.gSig:handles.options.gSig);
        if int_x(1)<1
            int_x = int_x + 1 - int_x(1);
        end
        if int_x(end)>handles.options.d1
            int_x = int_x - (int_x(end)-handles.options.d1);
        end
        int_y = round(handles.ROIvars.cm(end,2)) + (-handles.options.gSig:handles.options.gSig);
        if int_y(1)<1
            int_y = int_y + 1 - int_y(1);
        end
        if int_y(end)>handles.options.d2
            int_y = int_y - (int_y(end)-handles.options.d2);
        end
        [INT_x,INT_y] = meshgrid(int_x,int_y);
        coor = sub2ind([handles.options.d1,handles.options.d2],INT_x(:),INT_y(:));
        Ypatch = reshape(handles.Y(int_x,int_y,:),(2*handles.options.gSig+1)^2,size(handles.ROIvars.C,2));                        
        Yres = Ypatch - handles.A(coor,:)*handles.ROIvars.C;
        Yres = bsxfun(@minus, Yres, median(Yres,2));
        [atemp, ctemp, ~, ~, ~, ~] = greedyROI(reshape(Yres,2*handles.options.gSig+1,...
            2*handles.options.gSig+1,handles.T), 1, handles.options);
        %[atemp, ctemp] = initialize_components(reshape(Y_res,2*sx+1,2*sx+1,T), 1,sx,options);  % initialize
        % find contour
        a_srt = sort(atemp,'descend');
        ff = find(cumsum(a_srt.^2) >= handles.options.cont_threshold*sum(a_srt.^2),1,'first');
        K = size(handles.A,2) + 1;
        %maj of the general values of the GUI
        handles.accp(K) = 0;
        handles.newneur(K) = 1;
        handles.keep(K)=1;
        handles.disc(K)=0;
        %creating the new contours
        handles.A(coor,K) = atemp/norm(atemp);
        handles.CC{K} = contour(reshape(handles.A(:,end),handles.options.d1,handles.options.d2),[0,0]+a_srt(ff),'Linecolor',[1,0,1]/2);
        handles.CC{K}(handles.CC{K}<1) = NaN;
        handles.CC{K}(:,handles.CC{K}(1,:)>handles.options.d2) = NaN;
        handles.CC{K}(:,handles.CC{K}(2,:)>handles.options.d1) = NaN;
        handles.ROIvars.C(K,:) = ctemp*norm(atemp);
        new_center = com(handles.A(:,end),handles.options.d1,handles.options.d2);
        handles.ROIvars.cm(end,:) = new_center;
        handles.hasnew =true;
        handles = replot(handles,eventdata);
        guidata(hObject,handles);
    end
end


%% when a click on the figure has been made 
% will look for a close components and display all its informations
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
K = size(handles.A,2);
if ~ handles.selection
    for i = 1:K
        if (handles.keep(i)==1 && not(handles.isdiscarded)) || (handles.keep(i)==0 && handles.isdiscarded)
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
    if handles.newneur(idx) == 1
        set(handles.unrefinedtext,'Visible','on');
        set(handles.choosetext,'Visible','on');
        axes(figc);
        plot(handles.ROIvars.C(idx,:));
    else
        set(handles.unrefinedtext,'Visible','off');
        set(handles.choosetext,'Visible','off');
        handles.traceplot = handles.ROIvars.C(idx,:);
        set(handles.num_disp,'String',num2str(idx));
        %set info neur    
        set(handles.timecorr,'String',num2str(handles.ROIvars.rval_space(idx)));
        set(handles.cellsize,'String',num2str(handles.ROIvars.sizeA(idx)));
        set(handles.spacecorr,'String',num2str(handles.ROIvars.rval_time(idx)));
        set(handles.fitness,'String',num2str(handles.ROIvars.fitness(idx)));
        axes(figc);
        if handles.isdeconvolved 
            T = size(handles.ROIvars.C,2);
            hold off;
            plot(1:T,handles.ROIvars.C(handles.idx,:),'--k'); hold on; plot(1:T,handles.ROIvars.Cdec(handles.idx,:),'r','linewidth',2);
            spt = find(handles.S(handles.idx,:));
            if spt(1) == 1; spt(1) = []; end
            hold on; scatter(spt,repmat(-0.25,1,length(spt)),'m*')
            title(['Component ',num2str(handles.idx)]);
            legend('Fluorescence DF/F','Deconvolved','Spikes')
        else
            plot(handles.ROIvars.C(idx,:));
        end
    end
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


%will switch between accepted and discaded components. and will switch the button type in accordance with it
function Switchmode_Callback(hObject, eventdata, handles)
if ~ handles.selection
    if handles.isdiscarded 
        handles.isdiscarded =false;
        set(handles.Switchmode,'String','discarded');
        set(handles.Switchmode,'ForegroundColor',[1.0,.0,.0]);
        set(handles.remove_accept,'String','BAD');
        set(handles.remove_accept,'ForegroundColor',[1.0,.0,.0]);
    else
        handles.isdiscarded =true;
        set(handles.Switchmode,'String','Accepted');
        set(handles.Switchmode,'ForegroundColor',[.0,1.0,.0]);
        set(handles.remove_accept,'String','GOOD');
        set(handles.remove_accept,'ForegroundColor',[.0,1.0,.0]);
    end
    handles = replot(handles,eventdata);
    guidata(hObject, handles);
end

%% discarding and accepting PART
%depending on what type of componenent is displayed, 
%will either discard or accept a component when it has been selected
function remove_accept_Callback(hObject, eventdata, handles)
if handles.idx && ~ handles.selection && ~ handles.isrefined
    if handles.isdiscarded
        display('accepting');
        handles.accp(handles.idx)=1;
        handles.disc(handles.idx)=0;
        handles = replot(handles,eventdata);
    else
        display('deletion');
        handles.disc(handles.idx)=1;
        handles.accp(handles.idx)=0;
        handles = replot(handles,eventdata);
    end
    guidata(hObject, handles);
end


%% The main function to show the components and do some computations that makes more sense here.
function [handles] = replot(handles,eventdata)
if ~ handles.selection
    if not(handles.isimple) % computing by comparing everything
        handles.keep = ((handles.ROIvars.rval_space > handles.options.space_thresh) & ...
            (handles.ROIvars.rval_time > handles.options.time_thresh) & ... 
            (handles.ROIvars.sizeA >= handles.options.min_size_thr) & ...
            (handles.ROIvars.sizeA <= handles.options.max_size_thr) & ...
            (handles.ROIvars.fitness <= handles.options.min_fitness) & ...
            (handles.ROIvars.fitness_delta <= handles.options.min_fitness_delta));    
        original = size(find(handles.newneur),1);
        if original % if new components have been added by the user manually and not refined yet (not comparable still)
            s = size(handles.keep);
            handles.keep(s:s+original) = handles.newneur(s:s+original);
        end
    end 
    handles.keep = or(and(handles.keep,not(handles.disc)),handles.accp);    %adding the hard accepting and discarding 
    if handles.isdiscarded %% show the discarding components or the accepted ones
       handles.A_keep = handles.A(:,not(handles.keep));
       axes(handles.template_fig);
       [~,json_file,im] = plot_contours(handles.A_keep,handles.template,...
     handles.options,handles.shownum,[],[],[],find(not(handles.keep)),handles.ROIvars.cm(not(handles.keep),:)); 
    else 
       handles.A_keep = handles.A(:,handles.keep);
       axes(handles.template_fig);
       [~,json_file,im] = plot_contours(handles.A_keep,handles.template,...
            handles.options,handles.shownum,[],[],[],find(handles.keep),handles.ROIvars.cm(handles.keep,:)); 
    end
    if handles.saving %% we will have a box to ask the user where to save the datas
        [filename, pathname] = uiputfile({'*.json'},'Save as');
        if isequal(filename,0) || isequal(pathname,0)
            disp('User selected Cancel')
        else
           savejson('jmesh',json_file,fullfile(pathname,filename)); 
           disp(['User selected ',fullfile(pathname,filename)])
        end
    end
    disp('replot') % this part is to add to the plotted image the function back.
    % else it will not work on the image but on the axes
    set(im,'ButtonDownFcn',@(hObject,eventdata)ROI_GUI('template_fig_ButtonDownFcn',hObject,eventdata,guidata(hObject)));
    handles.cellnum = sum(handles.keep);
    set(handles.computing_text,'Visible','off');
else % the procedure of defining a training dataset for the simple mode
    handles.A_keep = handles.A(:,handles.ranset(handles.i)); 
    axes(handles.template_fig);
    plot_contours(handles.A_keep,handles.template,handles.options,0,[],[],[],handles.ranset(handles.i),handles.ROIvars.cm(handles.ranset(handles.i),:)); 
    handles.cellnum = 1;
    set(handles.timecorr,'String',num2str(handles.ROIvars.rval_space(handles.ranset(handles.i))));
    set(handles.cellsize,'String',num2str(handles.ROIvars.sizeA(handles.ranset(handles.i))));
    set(handles.spacecorr,'String',num2str(handles.ROIvars.rval_time(handles.ranset(handles.i))));
    set(handles.fitness,'String',num2str(handles.ROIvars.fitness(handles.ranset(handles.i))));
    axes(handles.traceplot);
    plot(handles.ROIvars.C(handles.ranset(handles.i),:));
    set(handles.num_disp,'String',num2str(handles.ranset(handles.i)));
    set(handles.computing_text,'Visible','off');
end


%% FINISHING PIPELINE FUNCTIONS

% --- Executes on button press in refine.
% will either do a first cnmf for the all the components if some of them were just added
%or refine the accepted components to produce the final dataset to be treated. 
function refine_Callback(hObject, eventdata, handles) 
if ~ handles.selection && ~ handles.isrefined && ~ handles.isextracted && ~ handles.isdeconvolved 
    handles.shownum = 0;
    if handles.hasnew 
        %redo the pipeline
        set(handles.computing_text,'Visible','on');
        handles.Yr = reshape(handles.Y,handles.d,handles.T);
        [handles.A,handles.b,handles.ROIvars.C] = update_spatial_components(...
            handles.Yr,handles.ROIvars.C,handles.f,[handles.A,handles.b],handles.P,handles.options);
        % update temporal components
        handles.P.p = 0;    % set AR temporarily to zero for speed
        [handles.ROIvars.C,handles.f,handles.P,handles.S,handles.YrA] = update_temporal_components(...
            handles.Yr,handles.A,handles.b,handles.ROIvars.C,handles.f,handles.P,handles.options);

        [handles.ROIvars.rval_space,handles.ROIvars.rval_time,~,handles.ROIvars.sizeA,~] = classify_components(...
            handles.Y,handles.A,handles.ROIvars.C,handles.b,handles.f,handles.YrA,handles.options);
        
        traces = prctfilt(handles.ROIvars.C+handles.YrA,8,1000,100);
        handles.ROIvars.fitness = compute_event_exceptionality(traces,0);
        handles.ROIvars.fitness_delta = compute_event_exceptionality(diff(traces,[],2),0);
        handles.newneur = zeros(size(handles.cellnum));
        handles.hasnew =false;
        set(handles.computing_text,'Visible','off');
    else
        set(handles.computing_text,'Visible','on');
        handles.P.p = 2;    % restore AR value
        % A first merge 
        %check that there is an updated keep here
        [handles.A,handles.ROIvars.C,~,~,handles.P,handles.S] = merge_components(handles.Yr,handles.A(:,handles.keep),...
            handles.b,handles.ROIvars.C(handles.keep,:),handles.f,handles.P,handles.S,handles.options);

        [handles.A,handles.b,handles.ROIvars.C] = update_spatial_components(handles.Yr,handles.ROIvars.C,...
            handles.f,[handles.A,handles.b],handles.P,handles.options);

        [handles.ROIvars.C,handles.f,handles.P,handles.S,handles.YrA] = update_temporal_components(...
            handles.Yr,handles.A,handles.b,handles.ROIvars.C,handles.f,handles.P,handles.options);

        [handles.ROIvars.rval_space,handles.ROIvars.rval_time,~,handles.ROIvars.sizeA,handles.keep] = classify_components(...
            handles.Y,handles.A,handles.ROIvars.C,handles.b,handles.f,handles.YrA,handles.options);

        handles.cellnum = sum(handles.keep);
        handles.newneur = zeros(handles.cellnum);

        traces = prctfilt(handles.ROIvars.C+handles.YrA,8,1000,100);
        handles.ROIvars.fitness = compute_event_exceptionality(traces,0);
        handles.ROIvars.fitness_delta = compute_event_exceptionality(diff(traces,[],2),0);

        handles.disc = zeros(size(handles.keep));
        handles.accp = zeros(size(handles.keep));

        handles.isrefined = true;
        set(handles.refine,'ForegroundColor',[1.0,.0,.0]);
        set(handles.computing_text,'Visible','off');
    end
    display('refined')
    handles = replot(handles,eventdata);
    guidata(hObject,handles);
end 

%% will do the deconvolution part (first will extract the dff values)
% once this has been made it will get plotted inside the traces axes. 
function deconvolve_Callback(hObject, eventdata, handles) 
%deconvolve and extract
if ~ handles.selection && handles.isrefined && ~ handles.isdeconvolved 
    set(handles.computing_text,'Visible','on');
    handles.shownum = 0;

    F = diag(sum(handles.A_keep.^2)) * ...
        (handles.ROIvars.C(handles.keep,:) + handles.YrA(handles.keep,:));                      % fluorescence
    Fd = prctfilt(F,handles.df_percentile,handles.window);                      % detrended fluorescence
    Bc = prctfilt((handles.A_keep'*handles.b)*handles.f,30,1000,300,0) + (F-Fd);       % background + baseline for each component
    F_dff = Fd./Bc;
    K = size(handles.A_keep,2);
    handles.ROIvars.Cdec = zeros(size(handles.ROIvars.C(handles.keep,:)));
    handles.kernels = cell(K,1);
    for it = 1:K
        disp(it)
        [handles.ROIvars.Cdec(it,:),handles.S(it,:),handles.kernels{it}] = deconvCa(...
            F_dff(it,:), [], handles.min_spikes, true, false, [], 20, [], 0);
    end
    handles.isdeconvolved = true;
    handles = replot(handles,eventdata);
    display('deconvolved and extracted')
    set(handles.computing_text,'Visible','off');
    set(handles.deconvolve,'ForegroundColor',[1.0,.0,.0]);
    guidata(hObject,handles);
end


function save_Callback(hObject, eventdata, handles) %totest
if ~ handles.selection
    handles.saving = true;
    handles = replot(handles,eventdata);
    guidata(hObject,handles);
end


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



%% an additionnal function to free the simplemode one 
% will do the computations and save a function that allows to play with 
%the classification and define where new components needs to go. 
function handles = compute_simple(handles)
set(handles.computing_text,'Visible','on');%REVIEW
%normalizing the values
val = handles.ROIvars;
K = size(handles.A,2);
val.sizeA = val.sizeA/handles.MAXCELLSIZE;
val.fitness = exp(val.fitness);
val.fitness_delta = exp(val.fitness_delta);
val.rval_space(val.rval_space<0)=0;
val.rval_time(val.rval_time<0)=0;
infocomp = [val.rval_space,val.rval_time,val.sizeA,val.fitness,val.fitness_delta];
test = infocomp(setdiff(1:K,handles.ranset),:);
train = infocomp(handles.ranset,:);
group = handles.keep(handles.ranset);

%this should work with larger training set than the one of demoMovie with
% a sizetrain of (30-40) else there is not enough values ( and especially 
% with the fitness values we have here)
%else try  (but only if matlab is the latest)
% Mdl = fitcdiscr(train,group,'OptimizeHyperparameters','auto');
% handles.keep(handles.sizetrain+1:end) = predict(Mdl,test);
% it should produce better results
[C,err,~,~,coeff] = classify(test,train,group,'Quadratic'); 
handles.K = coeff(1,2).const;
handles.L = coeff(1,2).linear;
handles.Q = coeff(1,2).quadratic;
display(' error is : '+err)
%______

handles.keep(setdiff(1:K,handles.ranset))=C; % new keep values
handles.simplekeep =handles.keep;
% Function to compute K + L*v + v'*Q*v for multiple vectors
handles.F = @(v) handles.K + v*handles.L + sum((v*handles.Q) .* v, 2);
%feval(f,v)
%displaying the simple mode view
handles.isimple =true;
handles.iscomputed =true;
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
set(handles.text35,'Visible','off');
set(handles.text37,'Visible','off');
set(handles.text38,'Visible','off');
set(handles.text39,'Visible','off');
set(handles.text40,'Visible','off');
set(handles.text15,'Visible','off');
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

set(handles.slider_min,'Visible','off');
set(handles.slider_max,'Visible','off');
set(handles.slider_min_fit,'Visible','off');
set(handles.slider_min_fit_delta,'Visible','off');
set(handles.slider_rval_t,'Visible','off');
set(handles.slider_rval_sp,'Visible','off');

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
set(handles.Switchmode,'Visible','on');
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



