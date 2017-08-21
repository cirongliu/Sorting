function varargout = Modulation(varargin)
% MODULATION MATLAB code for Modulation.fig
%      MODULATION, by itself, creates a new MODULATION or raises the existing
%      singleton*.
%
%      H = MODULATION returns the handle to a new MODULATION or the handle to
%      the existing singleton*.
%
%      MODULATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODULATION.M with the given input arguments.
%
%      MODULATION('Property','Value',...) creates a new MODULATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Modulation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Modulation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Modulation

% Last Modified by GUIDE v2.5 11-Apr-2016 17:36:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Modulation_OpeningFcn, ...
                   'gui_OutputFcn',  @Modulation_OutputFcn, ...
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


% --- Executes just before Modulation is made visible.
function Modulation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Modulation (see VARARGIN)

%% Initial show
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.text, 'String', 'Running, the manual modulation is ready');
set(handles.text8, 'String', varargin{end});
handles.Input = varargin;%
handles.Backup = varargin;
handles.Clusters = zeros(1,7);%Label of clusters
handles.manual_sign = 0;%Label of drawing
handles.manual_label = [];%Points of drawing
if varargin{1}~=7
    handles.Clusters(varargin{1}+1:end) = NaN;%Label of clusters
    switch varargin{1}
        case 1
            set(handles.checkbox1,'Enable','On');
            set(handles.checkbox2,'Enable','Off');
            set(handles.checkbox3,'Enable','Off');
            set(handles.checkbox4,'Enable','Off');
            set(handles.checkbox5,'Enable','Off');
            set(handles.checkbox6,'Enable','Off');
            set(handles.checkbox7,'Enable','Off');
        case 2
            set(handles.checkbox1,'Enable','On');
            set(handles.checkbox2,'Enable','On');
            set(handles.checkbox3,'Enable','Off');
            set(handles.checkbox4,'Enable','Off');
            set(handles.checkbox5,'Enable','Off');
            set(handles.checkbox6,'Enable','Off');
            set(handles.checkbox7,'Enable','Off');
        case 3
            set(handles.checkbox1,'Enable','On');
            set(handles.checkbox2,'Enable','On');
            set(handles.checkbox3,'Enable','On');
            set(handles.checkbox4,'Enable','Off');
            set(handles.checkbox5,'Enable','Off');
            set(handles.checkbox6,'Enable','Off');
            set(handles.checkbox7,'Enable','Off');
        case 4
            set(handles.checkbox1,'Enable','On');
            set(handles.checkbox2,'Enable','On');
            set(handles.checkbox3,'Enable','On');
            set(handles.checkbox4,'Enable','On');
            set(handles.checkbox5,'Enable','Off');
            set(handles.checkbox6,'Enable','Off');
            set(handles.checkbox7,'Enable','Off');
        case 5
            set(handles.checkbox1,'Enable','On');
            set(handles.checkbox2,'Enable','On');
            set(handles.checkbox3,'Enable','On');
            set(handles.checkbox4,'Enable','On');
            set(handles.checkbox5,'Enable','On');
            set(handles.checkbox6,'Enable','Off');
            set(handles.checkbox7,'Enable','Off');
        case 6
            set(handles.checkbox1,'Enable','On');
            set(handles.checkbox2,'Enable','On');
            set(handles.checkbox3,'Enable','On');
            set(handles.checkbox4,'Enable','On');
            set(handles.checkbox5,'Enable','On');
            set(handles.checkbox6,'Enable','On');
            set(handles.checkbox7,'Enable','Off');
        case 7
            set(handles.checkbox1,'Enable','On');
            set(handles.checkbox2,'Enable','On');
            set(handles.checkbox3,'Enable','On');
            set(handles.checkbox4,'Enable','On');
            set(handles.checkbox5,'Enable','On');
            set(handles.checkbox6,'Enable','On');
            set(handles.checkbox7,'Enable','On');
    end
end
handles.selection = 0;
handles.collection = [];
BOX_shows(handles)
%radio
set(handles.radio1,'Enable','Off');
set(handles.radio2,'Enable','Off');
set(handles.radio3,'Enable','Off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Choose default command line output for Modulation
% handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Modulation wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Modulation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.temp_cluster(:,1);
varargout{2} = handles.temp_cluster_label;
% assignin('base','Final_cluster_label',handles.temp_cluster(:,1));
% assignin('base','Final_cluster_type',temp_cluster_label);
delete(handles.figure1);

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
Index = get(hObject,'Value');%Value
if Index == 1
    handles.Clusters(1) = 1;
    set(handles.text, 'String', 'Showing the Cluster 1');
    guidata(hObject, handles)
    BOX_shows(handles);
    spike_properties(handles)
else
    handles.Clusters(1) = 0;
    set(handles.text, 'String', []);
    guidata(hObject, handles)
    BOX_shows(handles);   
end
guidata(hObject, handles)


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
Index = get(hObject,'Value');%Value
if Index == 1
    handles.Clusters(2) = 1;
    set(handles.text, 'String', 'Showing the Cluster 2');
    guidata(hObject, handles)
    BOX_shows(handles);  
    spike_properties(handles)
else
    handles.Clusters(2) = 0;
    set(handles.text, 'String', []);
    guidata(hObject, handles)
    BOX_shows(handles); 
end
guidata(hObject, handles)

% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
Index = get(hObject,'Value');%Value
if Index == 1
    handles.Clusters(3) = 1; 
    set(handles.text, 'String', 'Showing the Cluster 3');
    guidata(hObject, handles)
    BOX_shows(handles);
    spike_properties(handles)
else
    handles.Clusters(3) = 0;  
    set(handles.text, 'String', []);
    guidata(hObject, handles)
    BOX_shows(handles);
end
guidata(hObject, handles)

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
Index = get(hObject,'Value');%Value
if Index == 1
    handles.Clusters(4) = 1;
    set(handles.text, 'String', 'Showing the Cluster 4');
    guidata(hObject, handles)
    BOX_shows(handles);   
    spike_properties(handles)
else
    handles.Clusters(4) = 0;
    set(handles.text, 'String', []);
    guidata(hObject, handles)
    BOX_shows(handles);   
end
guidata(hObject, handles)

% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
Index = get(hObject,'Value');%Value
if Index == 1
    handles.Clusters(5) = 1;
    set(handles.text, 'String', 'Showing the Cluster 5');
    guidata(hObject, handles)
    BOX_shows(handles); 
    spike_properties(handles);
else
    handles.Clusters(5) = 0;
    set(handles.text, 'String', []);
    guidata(hObject, handles)
    BOX_shows(handles);   
end

% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
Index = get(hObject,'Value');%Value
if Index == 1
    handles.Clusters(6) = 1;
    set(handles.text, 'String', 'Showing the Cluster 6');
    guidata(hObject, handles)
    BOX_shows(handles); 
    spike_properties(handles)
else
    handles.Clusters(6) = 0;
    set(handles.text, 'String', []);
    guidata(hObject, handles)
    BOX_shows(handles); 
end
guidata(hObject, handles)

% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7
Index = get(hObject,'Value');%Value
if Index == 1
    handles.Clusters(7) = 1;
    set(handles.text, 'String', 'Showing the Cluster 7');
    guidata(hObject, handles)
    BOX_shows(handles);
    spike_properties(handles);
else
    handles.Clusters(7) = 0;
    set(handles.text, 'String', []);
    guidata(hObject, handles)
    BOX_shows(handles); 
end
guidata(hObject, handles)

% --- Executes on button press in pushbutton4.
% --- Mouse gets the circle
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Index = get(hObject,'Value');%Value
if Index == 1    
    handles.manual_sign = handles.manual_sign+1;
    if mod(handles.manual_sign,2) == 1%First time
        %radio
        set(handles.radio1,'Enable','On');
        set(handles.radio2,'Enable','On');
        set(handles.radio3,'Enable','On');
        set(handles.text, 'String', 'Select a space');
        guidata(hObject, handles);
    else
        %radio
        set(handles.radio1,'Enable','On');
        set(handles.radio2,'Enable','On');
        set(handles.radio3,'Enable','On');
        set(handles.radio1,'Value',0);
        set(handles.radio2,'Value',0);
        set(handles.radio3,'Value',0);
        set(handles.text, 'String', 'Select a space again');
        set(handles.axes1,'ButtonDownFcn','');
        set(handles.axes2,'ButtonDownFcn','');
        set(handles.axes3,'ButtonDownFcn','');
        handles.manual_label = [];
        set(gcf,'pointer','arrow')
        guidata(hObject, handles);
        BOX_shows(handles);
    end
    guidata(hObject, handles);
end


% --- Executes on button press in pushbutton1.
% --- Selection function
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%Using for the cluster selection
Index = get(hObject,'Value');%Value
if Index == 1
    if isempty( find( handles.Clusters == 1 ) ) == 1
        set(handles.text, 'String', 'The selection of cluster is empty and try again');
    else
        if handles.selection == 0
            handles.temp_cluster = [zeros(length(handles.Input{2}),1),handles.Input{2}];
            handles.temp_cluster_pattern = handles.Input{3};
            set(handles.text, 'String', 'The cluster is selected');
            temp = find( handles.Clusters == 1 );%find the which cluster is selected
            New_label = [];
            for i = 1:length(temp)
                label_temp = find(handles.temp_cluster(:,2) == temp(i));
                New_label = [New_label;label_temp];
            end
            clear label_temp
            handles.selection = handles.selection+1;% Type + 1
            handles.temp_cluster(New_label,1) = handles.selection;
            %Feature
            handles.Input{2}(New_label,:) = NaN;
            handles.Input{3}(New_label,:) = NaN;
            for i = 1: length(temp)
                if temp(i) == 1
                    set(handles.checkbox1,'Enable','Off');
                end
                if temp(i) == 2
                    set(handles.checkbox2,'Enable','Off');
                end
                if temp(i) == 3
                    set(handles.checkbox3,'Enable','Off');
                end
                if temp(i) == 4
                    set(handles.checkbox4,'Enable','Off');
                end
                if temp(i) == 5
                    set(handles.checkbox5,'Enable','Off');
                end
                if temp(i) == 6
                    set(handles.checkbox6,'Enable','Off');
                end
                if temp(i) == 7
                    set(handles.checkbox7,'Enable','Off');
                end
            end
            button=questdlg('PLS select type','Type Selection','Single unit','Multi-unit','default');
            if strcmp(button,'Single unit')
                if handles.selection < 10
                    handles.collection = [handles.collection;['Single 0',num2str(handles.selection)]];
                else
                    handles.collection = [handles.collection;['Single ',num2str(handles.selection)]];
                end
                set(handles.listbox,'String',handles.collection);
            else
                if handles.selection < 10
                    handles.collection = [handles.collection;['Multip 0',num2str(handles.selection)]];
                else
                    handles.collection = [handles.collection;['Multip ',num2str(handles.selection)]];
                end
                set(handles.listbox,'String',handles.collection);
            end
            handles.Clusters(temp) = NaN;
            guidata(hObject, handles)
            BOX_shows(handles)
        else
            set(handles.text, 'String', 'The cluster is selected');
            temp = find( handles.Clusters == 1 );%find the which cluster is selected
            New_label = [];
            for i = 1:length(temp)
                label_temp = find(handles.temp_cluster(:,2) == temp(i));
                New_label = [New_label;label_temp];
            end
            clear label_temp
            handles.selection = handles.selection+1;% Type + 1
            handles.temp_cluster(New_label,1) = handles.selection;
            handles.Input{2}(New_label,:) = NaN;
            handles.Input{3}(New_label,:) = NaN;
            for i = 1: length(temp)
                if temp(i) == 1
                    set(handles.checkbox1,'Enable','Off');
                end
                if temp(i) == 2
                    set(handles.checkbox2,'Enable','Off');
                end
                if temp(i) == 3
                    set(handles.checkbox3,'Enable','Off');
                end
                if temp(i) == 4
                    set(handles.checkbox4,'Enable','Off');
                end
                if temp(i) == 5
                    set(handles.checkbox5,'Enable','Off');
                end
                if temp(i) == 6
                    set(handles.checkbox6,'Enable','Off');
                end
                if temp(i) == 7
                    set(handles.checkbox7,'Enable','Off');
                end
            end
            button=questdlg('PLS select type','Type Selection','Single unit','Multi-unit','default');
            if strcmp(button,'Single unit')
                if handles.selection < 10
                    handles.collection = [handles.collection;['Single 0',num2str(handles.selection)]];
                else
                    handles.collection = [handles.collection;['Single ',num2str(handles.selection)]];
                end
                set(handles.listbox,'String',handles.collection);
            else
                if handles.selection < 10
                    handles.collection = [handles.collection;['Multip 0',num2str(handles.selection)]];
                else
                    handles.collection = [handles.collection;['Multip ',num2str(handles.selection)]];
                end
                set(handles.listbox,'String',handles.collection);
            end
%             handles.collection = [handles.collection;['Neuron ',num2str(handles.selection)]];
%             set(handles.listbox,'String',handles.collection);
            handles.Clusters(temp) = NaN; 
            guidata(hObject, handles)
            BOX_shows(handles)
            set(handles.text, 'String', ['Neuron ',num2str(handles.selection),' is selected']);
        end      
        guidata(hObject, handles)
    end
end


% --- Executes on button press in pushbutton2.
% --- Re-set function
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Index = get(hObject,'Value');%Value
if Index == 1
    set(handles.checkbox1,'Enable','On');
    set(handles.checkbox2,'Enable','On');
    set(handles.checkbox3,'Enable','On');
    set(handles.checkbox4,'Enable','On');
    set(handles.checkbox5,'Enable','On');
    set(handles.checkbox6,'Enable','On');
    set(handles.checkbox7,'Enable','On');
    
    set(handles.checkbox1,'Value',0);
    set(handles.checkbox2,'Value',0);
    set(handles.checkbox3,'Value',0);
    set(handles.checkbox4,'Value',0);
    set(handles.checkbox5,'Value',0);
    set(handles.checkbox6,'Value',0);
    set(handles.checkbox7,'Value',0);
    
    handles.Input = handles.Backup;%
    handles.Clusters = zeros(1,7);%Label of clusters
    if handles.Backup{1}~=7
        handles.Clusters(handles.Backup{1}+1:end) = NaN;%Label of clusters
        for i = (handles.Backup{1}+1) : 7
            if i == 1
                set(handles.checkbox1,'Enable','Off');
            end
            if i == 2
                set(handles.checkbox2,'Enable','Off');
            end
            if i == 3
                set(handles.checkbox3,'Enable','Off');
            end
            if i == 4
                set(handles.checkbox4,'Enable','Off');
            end
            if i == 5
                set(handles.checkbox5,'Enable','Off');
            end
            if i == 6
                set(handles.checkbox6,'Enable','Off');
            end
            if i == 7
                set(handles.checkbox7,'Enable','Off');
            end
        end
    end
    handles.selection = 0;
    handles.collection = [];
    clear handles.temp_cluster temp_cluster    
    set(handles.listbox,'String',[]);
    set(handles.text, 'String', 'Restarting again');
    guidata(hObject, handles)
    
    BOX_shows(handles)
end


% --- Executes on button press in pushbutton3.
% --- Re-cluster algorithm
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Index = get(hObject,'Value');%Value
if Index == 1
    %Pattern recluster
    if isfield(handles,'temp_cluster') == 1
        temp = find(handles.temp_cluster(:,1) == 0);
        NewPattern = handles.temp_cluster_pattern(temp,:);%Recluster
    else
        NewPattern = handles.Input{3};%Recluster
    end
    set(handles.text, 'String', 'Cluster is starting, PLS wait for a moment');
    
    
    
    if size(NewPattern,1)<1000
        opts.p = size(NewPattern,1);
    else
        opts.p = 1000;
    end
    NumCluster = [1:7];
    %OLPP
    Value_SPE = [];
    parfor k = 1:length(NumCluster)
        temptemp = LSC(NewPattern,k,opts);
        Value_SPE = [Value_SPE,temptemp];
    end
    clear temptemp k    
    %%%%%% CalinskiHarabasz Index
    eva = evalclusters(NewPattern,Value_SPE,'CalinskiHarabasz');%Spectral-cluster
    Max_gap = eva.OptimalK;
    label_cluster = Value_SPE(:,Max_gap);
      
%     NumCluster = [1:7];
%     Value = [];
%     parfor k = 1:length(NumCluster)
%         [New_Fea] = Landmark_selection(NewPattern,k);%Land-mark selection
%         eva = evalclusters(New_Fea,'kmeans','gap','KList',[k]);%Spectral-cluster
%         Value = [Value;eva.CriterionValues];
%     end
%     Max_gap = find( Value==max(Value));
%     label_cluster = LSC(NewPattern,Max_gap);%Land-spectrum cluster
    
    set(handles.text, 'String', 'Clustering is finished');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Re-organize
    for i = 1:Max_gap
        label = find(label_cluster == i);
        handles.Input{2}(temp(label),1) = i;
    end
    handles.Input{1} = Max_gap;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.checkbox1,'Enable','On');
    set(handles.checkbox2,'Enable','On');
    set(handles.checkbox3,'Enable','On');
    set(handles.checkbox4,'Enable','On');
    set(handles.checkbox5,'Enable','On');
    set(handles.checkbox6,'Enable','On');
    set(handles.checkbox7,'Enable','On');
    
    set(handles.checkbox1,'Value',0);
    set(handles.checkbox2,'Value',0);
    set(handles.checkbox3,'Value',0);
    set(handles.checkbox4,'Value',0);
    set(handles.checkbox5,'Value',0);
    set(handles.checkbox6,'Value',0);
    set(handles.checkbox7,'Value',0);
    
    handles.Clusters = zeros(1,7);%Label of clusters
    if Max_gap~=7
        handles.Clusters(Max_gap+1:end) = NaN;%Label of clusters
        for i = (Max_gap+1) : 7
            if i == 1
                set(handles.checkbox1,'Enable','Off');
            end
            if i == 2
                set(handles.checkbox2,'Enable','Off');
            end
            if i == 3
                set(handles.checkbox3,'Enable','Off');
            end
            if i == 4
                set(handles.checkbox4,'Enable','Off');
            end
            if i == 5
                set(handles.checkbox5,'Enable','Off');
            end
            if i == 6
                set(handles.checkbox6,'Enable','Off');
            end
            if i == 7
                set(handles.checkbox7,'Enable','Off');
            end
        end
    end
    handles.temp_cluster(:,2) = handles.Input{2};
    guidata(hObject, handles)
    BOX_shows(handles)
    guidata(hObject, handles)
end



% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Index = get(hObject,'Value');%Value
if Index == 1
    set(handles.text, 'String', 'New cluster is generated');
    set(handles.axes1,'ButtonDownFcn','');
    set(handles.axes2,'ButtonDownFcn','');
    set(handles.axes3,'ButtonDownFcn','');
    set(handles.radio1,'Enable','Off');
    set(handles.radio2,'Enable','Off');
    set(handles.radio3,'Enable','Off'); 
    set(handles.radio1,'Value',0);
    set(handles.radio2,'Value',0);
    set(handles.radio3,'Value',0);
    set(gcf,'pointer','arrow');
    if handles.axes_Num == 1
        temp_point_X = handles.Input{3}(:,1);temp_point_Y = handles.Input{3}(:,2);
        X = handles.manual_label(:,1);Y = handles.manual_label(:,2);
        in = inpolygon(temp_point_X,temp_point_Y,X,Y);
    end
    if handles.axes_Num == 2
        temp_point_X = handles.Input{3}(:,1);temp_point_Y = handles.Input{3}(:,3);
        X = handles.manual_label(:,1);Y = handles.manual_label(:,2);
        in = inpolygon(temp_point_X,temp_point_Y,X,Y);
    end
    if handles.axes_Num == 3
        temp_point_X = handles.Input{3}(:,2);temp_point_Y = handles.Input{3}(:,3);
        X = handles.manual_label(:,1);Y = handles.manual_label(:,2);
        in = inpolygon(temp_point_X,temp_point_Y,X,Y);
    end
    %Re-orgnize
    handles.Input{2}(in,:) = 8;    
    temp = unique(handles.Input{2});
    if isempty( isnan(temp) )~=1
        temp(isnan(temp)) = [];
    end
    for i = 1:length(temp)
        row = find(handles.Input{2}==temp(i));
        handles.Input{2}(row,:) = i;
        handles.Clusters(i) = 0; 
    end
    if handles.selection ~= 0%Update cluster
        handles.temp_cluster(:,2) = handles.Input{2};
    end
    if length(temp) == 1      
        set(handles.checkbox1,'Enable','On');
        set(handles.checkbox2,'Enable','Off');  
        set(handles.checkbox3,'Enable','Off');
        set(handles.checkbox4,'Enable','Off');
        set(handles.checkbox5,'Enable','Off');
        set(handles.checkbox6,'Enable','Off');
        set(handles.checkbox7,'Enable','Off');
        set(handles.checkbox1,'Value',0);
        set(handles.checkbox2,'Value',0);
        set(handles.checkbox3,'Value',0);
        set(handles.checkbox4,'Value',0);
        set(handles.checkbox5,'Value',0);
        set(handles.checkbox6,'Value',0);
        set(handles.checkbox7,'Value',0);
        handles.Clusters(2) = NaN;
        handles.Clusters(3) = NaN;
        handles.Clusters(4) = NaN;
        handles.Clusters(5) = NaN;
        handles.Clusters(6) = NaN;
        handles.Clusters(7) = NaN;
        handles.manual_label = [];
    end
    if length(temp) == 2
        set(handles.checkbox1,'Enable','On');
        set(handles.checkbox2,'Enable','On');
        set(handles.checkbox3,'Enable','Off');
        set(handles.checkbox4,'Enable','Off');
        set(handles.checkbox5,'Enable','Off');
        set(handles.checkbox6,'Enable','Off');
        set(handles.checkbox7,'Enable','Off');
        set(handles.checkbox1,'Value',0);
        set(handles.checkbox2,'Value',0);
        set(handles.checkbox3,'Value',0);
        set(handles.checkbox4,'Value',0);
        set(handles.checkbox5,'Value',0);
        set(handles.checkbox6,'Value',0);
        set(handles.checkbox7,'Value',0);
        handles.Clusters(3) = NaN;
        handles.Clusters(4) = NaN;
        handles.Clusters(5) = NaN;
        handles.Clusters(6) = NaN;
        handles.Clusters(7) = NaN;
        handles.manual_label = [];
    end
    if length(temp) == 3
        set(handles.checkbox1,'Enable','On');
        set(handles.checkbox2,'Enable','On');
        set(handles.checkbox3,'Enable','On');
        set(handles.checkbox4,'Enable','Off');
        set(handles.checkbox5,'Enable','Off');
        set(handles.checkbox6,'Enable','Off');
        set(handles.checkbox7,'Enable','Off');
        set(handles.checkbox1,'Value',0);
        set(handles.checkbox2,'Value',0);
        set(handles.checkbox3,'Value',0);
        set(handles.checkbox4,'Value',0);
        set(handles.checkbox5,'Value',0);
        set(handles.checkbox6,'Value',0);
        set(handles.checkbox7,'Value',0);
        handles.Clusters(4) = NaN;
        handles.Clusters(5) = NaN;
        handles.Clusters(6) = NaN;
        handles.Clusters(7) = NaN;
        handles.manual_label = [];
    end
    if length(temp) == 4
        set(handles.checkbox1,'Enable','On');
        set(handles.checkbox2,'Enable','On');
        set(handles.checkbox3,'Enable','On');
        set(handles.checkbox4,'Enable','On');
        set(handles.checkbox5,'Enable','Off');
        set(handles.checkbox6,'Enable','Off');
        set(handles.checkbox7,'Enable','Off');
        set(handles.checkbox1,'Value',0);
        set(handles.checkbox2,'Value',0);
        set(handles.checkbox3,'Value',0);
        set(handles.checkbox4,'Value',0);
        set(handles.checkbox5,'Value',0);
        set(handles.checkbox6,'Value',0);
        set(handles.checkbox7,'Value',0);
        handles.Clusters(5) = NaN;
        handles.Clusters(6) = NaN;
        handles.Clusters(7) = NaN;
        handles.manual_label = [];
    end
    if length(temp) == 5
        set(handles.checkbox1,'Enable','On');
        set(handles.checkbox2,'Enable','On');
        set(handles.checkbox3,'Enable','On');
        set(handles.checkbox4,'Enable','On');
        set(handles.checkbox5,'Enable','On');
        set(handles.checkbox6,'Enable','Off');
        set(handles.checkbox7,'Enable','Off');
        set(handles.checkbox1,'Value',0);
        set(handles.checkbox2,'Value',0);
        set(handles.checkbox3,'Value',0);
        set(handles.checkbox4,'Value',0);
        set(handles.checkbox5,'Value',0);
        set(handles.checkbox6,'Value',0);
        set(handles.checkbox7,'Value',0);
        handles.Clusters(6) = NaN;
        handles.Clusters(7) = NaN;
        handles.manual_label = [];
    end
    if length(temp) == 6
        set(handles.checkbox1,'Enable','On');
        set(handles.checkbox2,'Enable','On');
        set(handles.checkbox3,'Enable','On');
        set(handles.checkbox4,'Enable','On');
        set(handles.checkbox5,'Enable','On');
        set(handles.checkbox6,'Enable','On');
        set(handles.checkbox7,'Enable','Off');
        set(handles.checkbox1,'Value',0);
        set(handles.checkbox2,'Value',0);
        set(handles.checkbox3,'Value',0);
        set(handles.checkbox4,'Value',0);
        set(handles.checkbox5,'Value',0);
        set(handles.checkbox6,'Value',0);
        set(handles.checkbox7,'Value',0);
        handles.Clusters(7) = NaN;
        handles.manual_label = [];
    end 
    if length(temp) == 7
        set(handles.checkbox1,'Enable','On');
        set(handles.checkbox2,'Enable','On');
        set(handles.checkbox3,'Enable','On');
        set(handles.checkbox4,'Enable','On');
        set(handles.checkbox5,'Enable','On');
        set(handles.checkbox6,'Enable','On');
        set(handles.checkbox7,'Enable','On');
        set(handles.checkbox1,'Value',0);
        set(handles.checkbox2,'Value',0);
        set(handles.checkbox3,'Value',0);
        set(handles.checkbox4,'Value',0);
        set(handles.checkbox5,'Value',0);
        set(handles.checkbox6,'Value',0);
        set(handles.checkbox7,'Value',0);
        handles.manual_label = [];
    end
    BOX_shows(handles);   
    guidata(hObject, handles);
end



% --- Executes on button press in pushbutton4.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function BOX_shows(handles)
Clusters = handles.Clusters;
label_cluster = handles.Input{2};
Features = handles.Input{3};
%% Orginal clustering result
color = [[255 0 0];[255 165 0];[255 255 0];[0 255 0];[0 255 255];[0 0 255];[191 0 191]];
color = color./255;
color = single(color);
%%%%% X-Y axis
axes(handles.axes1)
cla reset
for i = 1:length(Clusters)
    if Clusters(i) == 0
        temp = find(label_cluster == i);
        hold on
        plot(Features(temp,1),Features(temp,2),'.','Color',color(i,:));
        hold off
    end
    if Clusters(i) == 1 % Selection
        temp = find(label_cluster == i);
        hold on
        plot(Features(temp,1),Features(temp,2),'.','Color',[0 0 0]);
        hold off
    end
    if isnan(Clusters(i)) == 1
    end
end
set(gca,'xtick',[])
set(gca,'ytick',[])
%%%%% X-Z axis
axes(handles.axes2)
cla reset
for i = 1:length(Clusters)
    if Clusters(i) == 0
        temp = find(label_cluster == i);
        hold on
        plot(Features(temp,1),Features(temp,3),'.','Color',color(i,:));
        hold off
    end
    if Clusters(i) == 1 % Selection
        temp = find(label_cluster == i);
        hold on
        plot(Features(temp,1),Features(temp,3),'.','Color',[0 0 0]);
        hold off
    end
    if isnan(Clusters(i)) == 1
    end
end
set(gca,'xtick',[])
set(gca,'ytick',[])
%%%%% Y-Z axis
axes(handles.axes3)
cla reset
for i = 1:length(Clusters)
    if Clusters(i) == 0
        temp = find(label_cluster == i);
        hold on
        plot(Features(temp,2),Features(temp,3),'.','Color',color(i,:));
        hold off
    end
    if Clusters(i) == 1 % Selection
        temp = find(label_cluster == i);
        hold on
        plot(Features(temp,2),Features(temp,3),'.','Color',[0 0 0]);
        hold off
    end
    if isnan(Clusters(i)) == 1
    end
end
set(gca,'xtick',[])
set(gca,'ytick',[])

% --- Executes on button press in Listbox.
function listbox_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radio1.
function radio1_Callback(hObject, eventdata, handles)
% hObject    handle to radio1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio1
Index = get(hObject,'Value');%Value
if Index == 1
    set(handles.radio1,'Enable','Off');
    set(handles.radio2,'Enable','Off');
    set(handles.radio3,'Enable','Off');
    axes(handles.axes1);
    handles.axes_Num = 1;
    set(gcf,'Pointer','cross');
    set(handles.text, 'String', 'Select the area where you want');
    set(handles.axes1,'ButtonDownFcn',{@axes_ButtonDownFcn,handles});
    guidata(hObject, handles)
end

% --- Executes on button press in radio2.
function radio2_Callback(hObject, eventdata, handles)
% hObject    handle to radio2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio2
Index = get(hObject,'Value');%Value
if Index == 1
    set(handles.radio1,'Enable','Off');
    set(handles.radio2,'Enable','Off');
    set(handles.radio3,'Enable','Off');    
    axes(handles.axes2)
    handles.axes_Num = 2;
    set(gcf,'Pointer','cross');
    set(handles.text, 'String', 'Select the area where you want');
    set(handles.axes2,'ButtonDownFcn',{@axes_ButtonDownFcn,handles});
    guidata(hObject, handles);
end

% --- Executes on button press in radio3.
function radio3_Callback(hObject, eventdata, handles)
% hObject    handle to radio3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio3
Index = get(hObject,'Value');%Value
if Index == 1
    set(handles.radio1,'Enable','Off');
    set(handles.radio2,'Enable','Off');    
    set(handles.radio3,'Enable','Off');
    axes(handles.axes3)
    handles.axes_Num = 3;
    set(gcf,'Pointer','cross');
    set(handles.text, 'String', 'Select the area where you want');
    set(handles.axes3,'ButtonDownFcn',{@axes_ButtonDownFcn,handles});
    guidata(hObject, handles);
end

% --- Executes on mouse press over axes background.
function axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.axes_Num == 1
    handles = guidata(handles.axes1);
end
if handles.axes_Num == 2
    handles = guidata(handles.axes1);
end 
if handles.axes_Num == 3
    handles = guidata(handles.axes1);
end

pt = get(gca,'CurrentPoint');
pt = [pt(1,1),pt(1,2)];
handles.manual_label = [handles.manual_label;pt];
hold on
plot(handles.manual_label(:,1),handles.manual_label(:,2),'-k*');
hold off
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Index = get(hObject,'Value');%Value
if Index == 1
    if isfield(handles,'temp_cluster') == 1
        if isempty( find(handles.temp_cluster==0) ) == 1
            guidata(hObject, handles);
            for i = 1:size(handles.collection,1)
                temp = find( handles.temp_cluster(:,1) ==  i); 
                if strcmp( handles.collection(i,1:6), 'Single') == 1
                    temp_cluster_label(temp,1) = 'S';
                else
                    temp_cluster_label(temp,1) = 'M';
                end
            end
            assignin('base','Final_cluster_label',handles.temp_cluster(:,1));
            assignin('base','Final_cluster_type',temp_cluster_label);
            handles.temp_cluster_label = temp_cluster_label;%Label
            guidata(hObject, handles);
            uiresume(handles.figure1);
%             close(handles.figure1);
        else
            set(handles.text, 'String', 'Attention!! Not complete yet');
            guidata(hObject, handles);
        end
    else
        set(handles.text, 'String', 'Attention!! Not complete yet');
        guidata(hObject, handles);
    end
end


% --- Executes during object deletion, before destroying properties.
function text_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to text7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function [t, c] = sc_acorr(a)
% function [l,c] =  sc_acorr(a,lag,N);
% trying it the easy way... may be too slow to work out, but the code was
% already written...

% defaults - need to figure out how to make specifiable
binwidth = 2;
win = [-200 200];
t = win(1):binwidth:win(2);

crosscorr_result = zeros(size(t));
for s = 1:length(a)
    x = a(a(:) >= a(s) + win(1) & a(:) <= a(s) + win(end));
    y = x - a(s);
    if ~isempty(x)
        n = histc(y(y~=0),t);
        if ~isempty(n)
            crosscorr_result = crosscorr_result + reshape(n,1,length(crosscorr_result));
        end
    end
end
t = t + binwidth/2;
c = crosscorr_result ./ (sum(crosscorr_result));

%Calculate the ISI interval, autocorrelation and waveform
function spike_properties(handles)
Clusters = handles.Clusters;
label_cluster = handles.Input{2};
candidate = find(Clusters == 1) ;
New_label = [];
Name = [];
for i = 1:length(candidate)
    New_label = [New_label;find(label_cluster == candidate(i))];
    if i == 1
        Name = [Name,num2str(candidate(i))];
    else
        Name = [Name,' & ',num2str(candidate(i))];
    end
end
if length(candidate) > 1
    figure,
    for i = 1:length(candidate)
        if i == 1
            temp_label{i} = find(label_cluster == candidate(i));
        else
            temp_label{i} = [temp_label{i-1};find(label_cluster == candidate(i))];
        end
        temp = handles.Input{5}(temp_label{i},1);
        temp = sort(temp);
        diff_temp = diff(temp);
        subplot(length(candidate),1,i),hist(diff_temp,0:1:1000);
        ax = gca;
        hold on
        plot([2 2],[0 ax.YLim(2)],'-r');
        plot([10 10],[0 ax.YLim(2)],'-g');
        hold off
        xlim([0 50]);
        if i == 1
            title(['ISI in Cluster ',num2str(candidate(1))]);
        else
            name = [];
            for j = 1:i
                name = [name,num2str(candidate(j)),' & '];       
            end
            title(['ISI in Cluster ',name]);
        end
    end
 
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     % Observed data
    %     n1 = length(find(diff_temp1<=2)); N1 = length(diff_temp1);
    %     n2 = length(find(diff_temp2<=2)); N2 = length(diff_temp2);
    %     % Pooled estimate of proportion
    %     p0 = (n1+n2) / (N1+N2);
    %     % Expected counts under H0 (null hypothesis)
    %     n10 = N1 * p0;
    %     n20 = N2 * p0;
    %     % Chi-square test, by hand
    %     observed = [n1 N1-n1 n2 N2-n2];
    %     expected = [n10 N1-n10 n20 N2-n20];
    %     chi2stat = sum((observed-expected).^2 ./ expected);
    %     p = 1 - chi2cdf(chi2stat,1);
    %     xlabel(['Within 10ms in Chi-square test and P value is ',num2str(p)]);
end
%Waveform
figure,
subplot(2,1,1)
template(1,:) = mean(handles.Input{4}(New_label,:));
template(2,:) = std(handles.Input{4}(New_label,:));
plotshaded(1:length(template(1,:)),[template(1,:)+template(2,:);template(1,:)-template(2,:)],'k');
hold on
plot(template(1,:),'-k');
% plot([0.3*40*2+1,0.3*40*2+1],[min(template(1,:)) max(template(1,:))],'--b');%Spike time
plot([0.2*40*2+1,0.2*40*2+1],[min(template(1,:)) max(template(1,:))],'--b');%Spike time, peak time
hold off
axis tight
set(gca,'xtick',[])
title(['Waveform in Cluster ',Name]);
% Autocorrelation
[t, c] = sc_acorr(handles.Input{5}(New_label,1));%Spike time
c(end) = [];t(end) = [];
subplot(2,1,2)
bar(t,c);
xlabel('Time lag (ms)');
title(['Autocorrelation in Cluster ',Name])
box off
set(gca,'ytick',[])
