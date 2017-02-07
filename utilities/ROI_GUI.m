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

% Last Modified by GUIDE v2.5 06-Feb-2017 13:07:59

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
%varargin = varargin{1};
handles.A = varargin{1};
handles.options = varargin{2};
handles.template = varargin{3};
handles.CC = varargin{4};
handles.keep = varargin{5};
handles.ROIvars = varargin{6};
handles.cellnum = sum(handles.keep);
handles.A_keep = handles.A(:,handles.keep);

set(handles.numcells,'string',num2str(handles.cellnum));
set(handles.slider_rvalt,'Value',handles.options.time_thresh);
set(handles.rval_t,'string',num2str(handles.options.time_thresh));
set(handles.slider_rvalsp,'Value',handles.options.space_thresh);
set(handles.rval_sp,'string',num2str(handles.options.space_thresh));
set(handles.slider_max,'Value',handles.options.max_size_thr);
set(handles.maxcellsize,'string',num2str(handles.options.max_size_thr));
set(handles.slider_min,'Value',handles.options.min_size_thr);
set(handles.mincellsize,'string',num2str(handles.options.min_size_thr));
axes(handles.template_fig);
colormap gray
plot_contours(handles.A,handles.template,handles.options,0,[],handles.CC,[],find(handles.keep)); 
guidata(hObject, handles);

uiwait()

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
function slider_rvalsp_Callback(hObject, eventdata, handles)
handles.options.space_thresh = get(hObject,'Value');
set(handles.rval_sp,'string',num2str(handles.options.space_thresh));
set(handles.computing_text,'Visible','on');
[handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars);
set(handles.numcells,'string',num2str(handles.cellnum));
set(handles.computing_text,'Visible','off');
guidata(hObject, handles);


%% TIME CORRELATION THRESHOLD
function slider_rvalt_Callback(hObject, eventdata, handles)
handles.options.time_thresh = get(hObject,'Value');
set(handles.rval_t,'string',num2str(handles.options.time_thresh));
set(handles.computing_text,'Visible','on');
[handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars);
set(handles.numcells,'string',num2str(handles.cellnum));
set(handles.computing_text,'Visible','off');
guidata(hObject, handles);



%% MAXIMUM CELL SIZE
function slider_max_Callback(hObject, eventdata, handles)
handles.options.max_size_thr = round(get(hObject,'Value'));
set(handles.maxcellsize,'string',num2str(handles.options.max_size_thr));
set(handles.computing_text,'Visible','on');
[handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars);
set(handles.numcells,'string',num2str(handles.cellnum));
set(handles.computing_text,'Visible','off');
guidata(hObject, handles);


%% MINIMUM CELL SIZE
function slider_min_Callback(hObject, eventdata, handles)
handles.options.min_size_thr = round(get(hObject,'Value'));
set(handles.mincellsize,'string',num2str(handles.options.min_size_thr));
set(handles.computing_text,'Visible','on');
[handles.A_keep,handles.keep,handles.cellnum] = replot(handles.A,handles.template,handles.options,handles.template_fig,handles.CC,handles.ROIvars);
set(handles.numcells,'string',num2str(handles.cellnum));
set(handles.computing_text,'Visible','off');
guidata(hObject, handles);

%% PLOT
function [A_keep,keep,cellnum] = replot(A,template,options,fig,CC,ROIvars)

keep = (ROIvars.rval_space > options.space_thresh) & (ROIvars.rval_time > options.time_thresh) & (ROIvars.max_pr > options.max_pr_thr) & (ROIvars.sizeA >= options.min_size_thr) & (ROIvars.sizeA <= options.max_size_thr);
A_keep = A(:,keep);    
axes(fig);
plot_contours(A_keep,template,options,0,[],CC,[],find(keep)); 
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

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

%% CREATE FUNCTIONS

function rval_sp_Callback(hObject, eventdata, handles)

function rval_t_Callback(hObject, eventdata, handles)

function mincellsize_Callback(hObject, eventdata, handles)

function maxcellsize_Callback(hObject, eventdata, handles)

function template_fig_CreateFcn(hObject, eventdata, handles)

function ID_CreateFcn(hObject, eventdata, handles)

function session_CreateFcn(hObject, eventdata, handles)

function numcells_CreateFcn(hObject, eventdata, handles)

function slider_rvalsp_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function rval_sp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function slider_rvalt_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function rval_t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slider_max_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function maxcellsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slider_min_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function mincellsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
