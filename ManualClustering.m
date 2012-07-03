function varargout = ManualClustering(varargin)
% MANUALCLUSTER M-file for ManualClustering.fig
%      MANUALCLUSTER, by itself, creates a new MANUALCLUSTER or raises the existing
%      singleton*.
%
%      H = MANUALCLUSTER returns the handle to a new MANUALCLUSTER or the handle to
%      the existing singleton*.
%
%      MANUALCLUSTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOSCLUSTER.M with the given input arguments.
%
%      MANUALCLUSTER('Property','Value',...) creates a new MOSCLUSTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ManualClustering_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ManualClustering_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ManualClustering

% Last Modified by GUIDE v2.5 01-Jul-2012 21:14:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ManualClustering_OpeningFcn, ...
                   'gui_OutputFcn',  @ManualClustering_OutputFcn, ...
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


% --- Executes just before ManualClustering is made visible.
function ManualClustering_OpeningFcn(hObject, eventdata, handles, model, filename)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ManualClustering (see VARARGIN)

% Choose default command line output for ManualClustering
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

assert(isa(model,'ClusteringHelper'), 'This must be passed with one or more clustering helpers');

if nargin > 4
    set(gcf,'Name',filename);
end

% display this data
handles.modelData = model;
mosSetFileButtons(hObject,handles,'off');
% Update handles structure
guidata(hObject, handles);
NewModel(hObject,handles);    

% UIWAIT makes ManualClustering wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ManualClustering_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.modelData;
delete(handles.figure1);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in opMerge.
function opMerge_Callback(hObject, eventdata, handles)
% hObject    handle to opMerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.modelData = merge(handles.modelData,GetSelectedIds(hObject,handles));
guidata(hObject,handles);
NewModel(hObject,handles);


% --- Executes on button press in opSplit.
function opSplit_Callback(hObject, eventdata, handles)
% hObject    handle to opSplit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.modelData = split(handles.modelData,GetSelectedIds(hObject,handles));
guidata(hObject,handles);
NewModel(hObject,handles);


% --- Executes on button press in opDelete.
function opDelete_Callback(hObject, eventdata, handles)
% hObject    handle to opDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.modelData = delete(handles.modelData,GetSelectedIds(hObject,handles));
guidata(hObject,handles);
NewModel(hObject,handles);


% --- Executes on button press in opReproject.
function opReproject_Callback(hObject, eventdata, handles)
% hObject    handle to opReproject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.modelData = reproject(handles.modelData,GetSelectedIds(hObject,handles));
guidata(hObject,handles);
NewModel(hObject,handles);


% --- Executes on button press in opStrip.
function opStrip_Callback(hObject, eventdata, handles)
% hObject    handle to opStrip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.modelData = strip(handles.modelData, GetSelectedIds(hObject,handles));
guidata(hObject,handles);
NewModel(hObject,handles);


% --- Executes on button press in opRefit.
function opRefit_Callback(hObject, eventdata, handles)
% hObject    handle to opRefit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.modelData = refit(handles.modelData); ..., GetSelectedIds(hObject,handles));
guidata(hObject,handles);
NewModel(hObject,handles);


% --- Executes on button press in opGroup.
function opGroup_Callback(hObject, eventdata, handles)
% hObject    handle to opGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.modelData = group(handles.modelData, GetSelectedIds(hObject,handles));
guidata(hObject,handles);
NewModel(hObject,handles);


% --- Executes on button press in opSingle.
function opSingle_Callback(hObject, eventdata, handles)
% hObject    handle to opSingle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.modelData = singleUnit(handles.modelData, GetSelectedIds(hObject,handles));
guidata(hObject,handles);
[clusIds groups] = getClusterIds(handles.modelData);
[fp fn snr frac] = getStats(handles.modelData);
su = hasTag(handles.modelData,'SingleUnit');
set(handles.stats,'Data',num2cell([clusIds' groups' fp' fn' snr' frac' su']));


% --- Executes on button press in opLDA.
function opLDA_Callback(hObject, eventdata, handles)
% hObject    handle to opLDA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure
plotLDAs(handles.modelData,'clusIds',GetSelectedIds(hObject,handles));


% --- Executes on button press in opTime.
function opTime_Callback(hObject, eventdata, handles)
% hObject    handle to opTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure
plotTimeFeatures(handles.modelData,'clusIds',GetSelectedIds(hObject,handles));


% --- Executes on button press in opPrev.
function opPrev_Callback(hObject, eventdata, handles)
% hObject    handle to opPrev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.fileNum > 1
    handles.fileNum = handles.fileNum - 1;
    mosLoadFileData(hObject, handles);
end


% --- Executes on button press in opSave.
function opSave_Callback(hObject, eventdata, handles)
% hObject    handle to opSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fn = handles.fileNames{handles.fileNum};
if strfind(fn,'clustered') == 1
    backup = strrep(fn,'clustered','backup');
    handles.fileNames{handles.fileNum} = strrep(fn,'clustered','finalized');
    if exist(fn,'file'), movefile(fn, backup); end
elseif strfind(fn,'finalized') == 1;
else
    error('Do not know how to save for this file name %s',fn);
end
data = handles.modelData;
data.Model = compress(data.Model);
save(handles.fileNames{handles.fileNum},'-v7.3','-struct','data');
set(handles.lblFilename,'String',handles.fileNames{handles.fileNum});
guidata(hObject,handles);


% --- Executes on button press in opNext.
function opNext_Callback(hObject, eventdata, handles)
% hObject    handle to opNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.fileNum < length(handles.fileNames)
    handles.fileNum = handles.fileNum + 1;
    mosLoadFileData(hObject, handles);
end


% --- Executes on selection change in opSelection.
function opSelection_Callback(hObject, eventdata, handles)
% hObject    handle to opSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns opSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from opSelection


% --- Executes during object creation, after setting all properties.
function opSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to opSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lbSelection.
function lbSelection_Callback(hObject, eventdata, handles)
% hObject    handle to lbSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lbSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lbSelection

UpdateDisplay(hObject,handles)


% --- Executes during object creation, after setting all properties.
function lbSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function mosSetFileButtons(hObjects,handles,status)
set(handles.opPrev,'Enable',status);
set(handles.opNext,'Enable',status);
set(handles.opSave,'Enable',status);


function clusIds = GetSelectedIds(hObject,handles)
clusters = get(handles.lbSelection,'String');
clusIds = cellfun(@str2num,clusters(get(handles.lbSelection,'Value')));


function mosLoadFileData(hObject,handles)
handles.modelData = load(handles.fileNames{handles.fileNum});
handles.modelData.Model = uncompress(handles.modelData.Model,handles.modelData);
if isfield(handles.modelData,'ClusterTags')
    handles.modelData = rmfield(handles.modelData,'ClusterTags');
end
set(handles.lblFilename,'String',handles.fileNames{handles.fileNum});
NewModel(hObject,handles)


function NewModel(hObject,handles)
% This function is called whenever the model changed.

handles.modelData = updateInformation(handles.modelData);
[handles.cc handles.cctime] = getCrossCorrs(handles.modelData);

% update cluster list and table with information
guidata(hObject, handles);
[clusIds groups] = getClusterIds(handles.modelData);
set(handles.lbSelection, 'String', num2cell(clusIds));
set(handles.lbSelection, 'Value', 1:length(clusIds));
[fp fn snr frac] = getStats(handles.modelData);
su = hasTag(handles.modelData, 'SingleUnit');
set(handles.stats, 'Data', num2cell([clusIds; groups; 100 * fp; 100 * fn; snr; 100 * frac; su]'));

% create CCG and waveform plots
if isfield(handles, 'ccg')
    delete(handles.ccg)
    delete(handles.wave)
    delete(handles.projAxes)
end
handles.ccg = plotCrossCorrs(handles.modelData, 'figure', hObject);
handles.wave = plotWaveforms(handles.modelData, 'figure', hObject);
[handles.projAxes, handles.projPlots] = plotProjections(handles.modelData, 'position', get(handles.projection, 'Position'));
guidata(hObject, handles);

UpdateDisplay(hObject,handles);


function UpdateDisplay(hObject,handles)
% update the display with all the selected information

% simple not-quite-safe lock to prevent race conditions
persistent locked
if isempty(locked), locked = false; end
if locked, return, end
locked = true; %#ok

% find out currently selected clusters
clusIds = GetSelectedIds(hObject, handles);
k = numel(clusIds);
[n, chans] = size(handles.wave);

% plot contamination
axes(handles.contamination);
[~, hdl] = plotContaminations(handles.modelData);
set(hdl, 'ButtonDownFcn', @ContamMatrixClicked);

% toggle projection plots
state = {'off', 'on'};
for i = 1 : n
    set(handles.projPlots(i, :), 'Visible', state{any(i == clusIds) + 1})
end

% switch on and resize waveform plots for selected units
pos = get(handles.wavepanel, 'Position');
for i = 1 : k
    for j = 1 : chans
        p = [pos(1) + (i - 1) * pos(3) / k, pos(2) + (j - 1) * pos(4) / chans, pos(3) / k, pos(4) / chans];
        set(handles.wave(clusIds(i), j), 'Position', p)
    end
end
for i = 1 : n
    for j = 1 : chans
        ShowPlot(handles.wave(i, j), ismember(i, clusIds));
    end
end

% switch on and resize CCG plots for selected units
pos = get(handles.ccgpanel, 'Position');
for i = 1 : k
    for j = 1 : k
        p = [pos(1) pos(2) 0 0] + [(i - 1) * pos(3), (k - j) * pos(4), pos(3), pos(4)] / k;
        set(handles.ccg(clusIds(i), clusIds(j)), 'Position', p)
    end
end
for i = 1 : n
    for j = 1 : n
        ShowPlot(handles.ccg(i, j), all(ismember([i j], clusIds)));
    end
end

% release lock
locked = false;


function ShowPlot(hdl, on)
% Show or hide plot

if nargin < 2 || on, state = 'on'; else state = 'off'; end
ch = get(hdl, 'Children');
if ~strcmp(state, get(ch(1), 'Visible'))
    set(ch, 'Visible', state)
end


function ContamMatrixClicked(hObject, eventdata)

handles = guidata(hObject);
cp = get(handles.contamination, 'CurrentPoint');
cp = round(cp(1, 1:2));
if strcmp(get(handles.figure1, 'SelectionType'), 'normal')
    set(handles.lbSelection, 'Value', unique(cp))
else
    sel = get(handles.lbSelection, 'Value');
    if all(ismember(cp, sel))
        set(handles.lbSelection, 'Value', setdiff(sel, cp))
    else
        set(handles.lbSelection, 'Value', union(sel, cp))
    end
end
UpdateDisplay(hObject,handles);


% --- Executes on button press in accepbutton.
function accepbutton_Callback(hObject, eventdata, handles)
% hObject    handle to accepbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over accepbutton.
function accepbutton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to accepbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

model = handles.modelData;
