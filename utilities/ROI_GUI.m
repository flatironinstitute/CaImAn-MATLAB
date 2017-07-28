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
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ROI_GUI

% Last Modified by GUIDE v2.5 26-Jul-2017 16:46:21

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
handles.A = varargin{1};
handles.options = varargin{2};
handles.template = varargin{3};
handles.CC = varargin{4};
handles.keep = varargin{5};
handles.ROIvars = varargin{6};
handles.cellnum = sum(handles.keep);
handles.idx = false;
%we will define C graph and info of the neuron once a click has been made
set(handles.min_fit_delta,'String',num2str(handles.options.min_fitness_delta));
set(handles.min_fit,'String',num2str(handles.options.min_fitness));
set(handles.mincellsize,'String',num2str(handles.options.min_size_thr));
set(handles.maxcellsize,'String',num2str(handles.options.max_size_thr));
set(handles.rval_sp,'String',num2str(handles.options.space_thresh));
set(handles.rval_t,'String',num2str(handles.options.time_thresh));
set(handles.numcells,'String',num2str(handles.cellnum));

set(handles.slider_min_fit_delta,'Value',handles.options.min_fitness_delta);
set(handles.slider_min_fit,'Value',handles.options.min_fitness);
set(handles.slider_min,'Value',handles.options.min_size_thr);
set(handles.slider_max,'Value',handles.options.max_size_thr);
set(handles.slider_rval_sp,'Value',handles.options.space_thresh);
set(handles.slider_rval_t,'Value',handles.options.time_thresh);

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
set(im,'ButtonDownFcn',@(hObject,eventdata,figc)ROI_GUI('template_fig_ButtonDownFcn',hObject,eventdata,handles.traceplot,guidata(hObject)));
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
handles.options.space_thresh = get(hObject,'Value');
set(handles.rval_sp,'String',num2str(handles.options.space_thresh));
set(handles.computing_text,'Visible','on');
[handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
set(handles.numcells,'String',num2str(handles.cellnum));
set(handles.ispace,'Value',0);
set(handles.computing_text,'Visible','off');
guidata(hObject, handles);


%% TIME CORRELATION THRESHOLD
function slider_rval_t_Callback(hObject, eventdata, handles)
handles.options.time_thresh = get(hObject,'Value');
set(handles.rval_t,'String',num2str(handles.options.time_thresh));
set(handles.computing_text,'Visible','on');
[handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
set(handles.numcells,'String',num2str(handles.cellnum));
set(handles.istimecorr,'Value',0);
set(handles.computing_text,'Visible','off');
guidata(hObject, handles);



%% MAXIMUM CELL SIZE
function slider_max_Callback(hObject, eventdata, handles)
handles.options.max_size_thr = round(get(hObject,'Value'));
set(handles.maxcellsize,'String',num2str(handles.options.max_size_thr));
set(handles.computing_text,'Visible','on');
[handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
set(handles.numcells,'String',num2str(handles.cellnum));
set(handles.ismax,'Value',0);
set(handles.computing_text,'Visible','off');
guidata(hObject, handles);


%% MINIMUM CELL SIZE
function slider_min_Callback(hObject, eventdata, handles)
handles.options.min_size_thr = round(get(hObject,'Value'));
set(handles.mincellsize,'String',num2str(handles.options.min_size_thr));
set(handles.computing_text,'Visible','on');
[handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
set(handles.numcells,'String',num2str(handles.cellnum));
set(handles.ismin,'Value',0);
set(handles.computing_text,'Visible','off');
guidata(hObject, handles);



%% MAX EVENT EXCEPTIONALITY
function slider_min_fit_Callback(hObject, eventdata, handles)

handles.options.min_fitness = get(hObject,'Value');
set(handles.min_fit,'String',num2str(handles.options.min_fitness));

set(handles.computing_text,'Visible','on');
[handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
set(handles.numcells,'String',num2str(handles.cellnum));
set(handles.isfit,'Value',0);
set(handles.computing_text,'Visible','off');
guidata(hObject, handles);

%% MAX EVENT EXCEPTIONALITY delta
% --- Executes on slider movement.
function slider_min_fit_delta_Callback(hObject, eventdata, handles)
% hObject    handle to slider_min_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.options.min_fitness_delta = get(hObject,'Value');
set(handles.min_fit_delta,'String',num2str(handles.options.min_fitness_delta));

set(handles.computing_text,'Visible','on');
[handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
set(handles.numcells,'String',num2str(handles.cellnum));
set(handles.isfitdelta,'Value',0);
set(handles.computing_text,'Visible','off');
guidata(hObject, handles);


%% check boxes
function ispace_Callback(hObject, eventdata, handles)
if get(handles.ispace,'Value')
    handles.options.space_thresh = 0;
    set(handles.rval_sp,'String',num2str(handles.options.space_thresh));
    
    set(handles.computing_text,'Visible','on');
    [handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
    set(handles.numcells,'String',num2str(handles.cellnum));
    set(handles.computing_text,'Visible','off');
    guidata(hObject, handles);
end

function istimecorr_Callback(hObject, eventdata, handles)
if get(handles.istimecorr,'Value')
    handles.options.time_thresh = 0;
    set(handles.rval_t,'String',num2str(handles.options.time_thresh));
    
    set(handles.computing_text,'Visible','on');
    [handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
    set(handles.numcells,'String',num2str(handles.cellnum));
    set(handles.computing_text,'Visible','off');
    guidata(hObject, handles);
end

function ismax_Callback(hObject, eventdata, handles)
if get(handles.ismax,'Value')
    handles.options.max_size_thr = 500;
    set(handles.maxcellsize,'String',num2str(handles.options.max_size_thr));
    set(handles.computing_text,'Visible','on');
    [handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
    set(handles.numcells,'String',num2str(handles.cellnum));
    set(handles.computing_text,'Visible','off');
    guidata(hObject, handles);
end

function ismin_Callback(hObject, eventdata, handles)
if get(handles.ismin,'Value')
    handles.options.min_fitness_delta = 5;
    set(handles.mincellsize,'String',num2str(handles.options.min_size_thr));

    set(handles.computing_text,'Visible','on');
    [handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
    set(handles.numcells,'String',num2str(handles.cellnum));
    set(handles.computing_text,'Visible','off');
    guidata(hObject, handles);
end

function isfit_Callback(hObject, eventdata, handles)
if get(handles.isfit,'Value')
    handles.options.min_fitness = 0;
    set(handles.min_fit,'String',num2str(handles.options.min_fitness));

    set(handles.computing_text,'Visible','on');
    [handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
    set(handles.numcells,'String',num2str(handles.cellnum));
    set(handles.computing_text,'Visible','off');
    guidata(hObject, handles);
end

function isfitdelta_Callback(hObject, eventdata, handles)
if get(handles.isfitdelta,'Value')
    handles.options.min_fitness_delta = 0;
    set(handles.min_fit_delta,'String',num2str(handles.options.min_fitness_delta));

    set(handles.computing_text,'Visible','on');
    [handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
    set(handles.numcells,'String',num2str(handles.cellnum));
    set(handles.computing_text,'Visible','off');
    guidata(hObject, handles);
end
%% click on the fig

function template_fig_ButtonDownFcn(hObject,eventdata, figc, handles)
set(handles.computing_text,'Visible','on');
%coordinates = ginput(1)
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

for i = 1:K(2)
    if handles.keep(i)==1
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
[handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
set(handles.numcells,'String',num2str(handles.cellnum));
set(handles.computing_text,'Visible','off');
set(hObject.Parent,'Units','normalized');
guidata(hObject, handles);


%% discarding

function pushbutton14_Callback(hObject, eventdata, handles)
if handles.idx
    display('deletion')
    handles.disc(handles.idx)=1;
    [handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars,handles.traceplot,handles);
end
guidata(hObject, handles);

%% PLOT
function [A_keep,keep,cellnum] = replot(A,template,options,fig,CC,ROIvars,figcx,handles)
keep = (ROIvars.rval_space > options.space_thresh) & (ROIvars.rval_time > options.time_thresh) & ... 
    (ROIvars.sizeA >= options.min_size_thr) & (ROIvars.sizeA <= options.max_size_thr) & ...
    (ROIvars.fitness <= options.min_fitness) & (ROIvars.fitness_delta <= options.min_fitness_delta) ...
    & not(handles.disc);
A_keep = A(:,keep);    

axes(fig);
[~,~,im] = plot_contours(A_keep,template,options,0,[],CC,[],find(keep)); 
set(im,'ButtonDownFcn',@(hObject,eventdata,figc)ROI_GUI('template_fig_ButtonDownFcn',hObject,eventdata,figcx,guidata(hObject)));
cellnum = sum(keep);



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
% eventdata  reserved - to be defined in a future version of MATLAB
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

function template_fig_CreateFcn(hObject, eventdata, handles)

function ID_CreateFcn(hObject, eventdata, handles)

function session_CreateFcn(hObject, eventdata, handles)

function numcells_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

%space correlation
function slider_rval_sp_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function rval_sp_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

%time correlation
function slider_rval_t_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);


function rval_t_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

%max cell size
function slider_max_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function maxcellsize_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

%min cell size
function slider_min_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);


function mincellsize_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

%min fit
function slider_min_fit_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function min_fit_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

%min fit delta
function slider_min_fit_delta_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function min_fit_delta_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

%infos about neuron
function timecorr_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function cellsize_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function spacecorr_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function fitness_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function traceplot_CreateFcn(hObject, eventdata, handles)

