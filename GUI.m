function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 12-Jun-2017 12:51:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;
handles.FileName = -1;
handles.FilePath = -1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in singleLineCBox.
function singleLineCBox_Callback(hObject, eventdata, handles)
% hObject    handle to singleLineCBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleLineCBox
handles.drawCurveCBox.Value = 1;
handles.drawTrajectoryCBox.Value = 1;
handles.exportCBox.Value = 1;
handles.clearCBox.Value = 1;


% --- Executes on button press in drawCurveCBox.
function drawCurveCBox_Callback(hObject, eventdata, handles)
% hObject    handle to drawCurveCBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of drawCurveCBox
if handles.drawCurveCBox.Value == 0
    handles.singleLineCBox.Value = 0;
end


% --- Executes on button press in drawTrajectoryCBox.
function drawTrajectoryCBox_Callback(hObject, eventdata, handles)
% hObject    handle to drawTrajectoryCBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of drawTrajectoryCBox
if handles.drawTrajectoryCBox.Value == 0
    handles.singleLineCBox.Value = 0;
end


% --- Executes on button press in exportCBox.
function exportCBox_Callback(hObject, eventdata, handles)
% hObject    handle to exportCBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exportCBox
if handles.exportCBox.Value == 0
    handles.clearCBox.Value = 0;
    handles.singleLineCBox.Value = 0;
end


% --- Executes on button press in clearCBox.
function clearCBox_Callback(hObject, eventdata, handles)
% hObject    handle to clearCBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of clearCBox
if handles.clearCBox.Value == 0
    handles.singleLineCBox.Value = 0;
end
if handles.clearCBox.Value == 1
    handles.exportCBox.Value = 1;
end


function frameRateEdit_Callback(hObject, eventdata, handles)
% hObject    handle to frameRateEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frameRateEdit as text
%        str2double(get(hObject,'String')) returns contents of frameRateEdit as a double
handles.timeFactorStatic.String = ['Zeitfaktor: ' num2str(round(str2double(handles.frameRateEdit.String)/str2double(handles.dataResolutionEdit.String),2))];


function dataResolutionEdit_Callback(hObject, eventdata, handles)
% hObject    handle to dataResolutionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataResolutionEdit as text
%        str2double(get(hObject,'String')) returns contents of dataResolutionEdit as a double
handles.timeFactorStatic.String = ['Zeitfaktor: ' num2str(round(str2double(handles.frameRateEdit.String)/str2double(handles.dataResolutionEdit.String),2))];


% --- Executes on button press in fileChooserBtn.
function fileChooserBtn_Callback(hObject, eventdata, handles)
% hObject    handle to fileChooserBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, filepath] = uigetfile('*.mat');
if filename ~= 0
    handles.FileName = filename;
    handles.FilePath = filepath;
    handles.fileChooserStatic.String = filename;
    guidata(hObject, handles);
end


% --- Executes on button press in renderBtn.
function renderBtn_Callback(hObject, eventdata, handles)
% hObject    handle to renderBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
render(handles);


%Function to render a single image or an animation showing differet aspects
%of the actors simulated behaviour.
function render(handles)
startTime = tic;

singleLineTrajectory = handles.singleLineCBox.Value;
drawBodyCurve = handles.drawCurveCBox.Value;
drawTrajectory = handles.drawTrajectoryCBox.Value;
doExport =  handles.exportCBox.Value;
clearFigure = handles.clearCBox.Value;
dataResolution = str2double(handles.dataResolutionEdit.String);
framerate = str2double(handles.frameRateEdit.String);

if ~drawBodyCurve && ~drawTrajectory
    msgbox('''Kurve zeichnen'' oder ''Bahnverlauf zeichnen'' muss gesetzt sein.','Error','error')
    return
end

if isequal(handles.FilePath, -1)
    msgbox('Keine Datei gewählt.','Error','error')
    return
end

config = load(strcat(handles.FilePath, handles.FileName));
config = config.config;

maxI = size(config, 1);
framenumber = round(maxI/1000*min([1000 dataResolution]));
indeces = round(linspace(1, maxI, framenumber));
f1 = figure('InnerPosition', [0 0 1280 720]);
colormap([0.3 0.3 0.3])
lightangle(40,15)
xlabel('x');ylabel('y');zlabel('z');
daspect([1 1 1]);
xlim([-60 60]);ylim([-60 60]);zlim([0 200]);
view(40,25);
hold on
if  drawTrajectory
    if ~drawBodyCurve
        xlim([-100 100]);ylim([-100 100]);zlim([140 180]);
        view(60,20)
    end
    headData = zeros(3,framenumber);
end
if drawBodyCurve
    bodyData = cell(1,framenumber);
end
i2 = 1;
tic
for i = indeces
    x = config(i,:);
    if drawTrajectory && ~drawBodyCurve
        headData(:,i2) = line2(x(1), x(2), x(3), x(4), x(5),0);
    elseif ~drawTrajectory && drawBodyCurve
        [~,bodyData{:,i2}] = line2(x(1), x(2), x(3), x(4), x(5),1);
    elseif drawBodyCurve && drawTrajectory
        [headData(:,i2),bodyData{:,i2}] = line2(x(1), x(2), x(3), x(4), x(5),1);
    end
    i2 = i2 + 1;
end
toc
if doExport
    %     v = VideoWriter('out3', 'MPEG-4'); %only use when using win7+ or MacOSX 10.7+
    v = VideoWriter('out3', 'Motion JPEG AVI'); % works on linux
    v.FrameRate = framerate;
    v.Quality = 100;
    open(v);
    meanLoopTime = 0;
    for i = 1:framenumber
        loopStartTime = tic;
        if clearFigure
            if drawTrajectory
                if singleLineTrajectory
                    plot3(headData(1,1:i),headData(2,1:i),headData(3,1:i),'-r');
                else
                    i2 = max(1, i-1);
                    plot3(headData(1,i2:i),headData(2,i2:i),headData(3,i2:i),'-r');
                end
            end
            if drawBodyCurve
                x = bodyData{:,i};
                h = streamtube({[x(1,:);x(2,:);x(3,:)]'}, 3);
                h.FaceLighting = 'gouraud';
                h.AmbientStrength = 0.8;
                h.DiffuseStrength = 0.7;
                h.SpecularStrength = 0.8;
                h.SpecularExponent = 15;
                h.BackFaceLighting = 'unlit';
                lightangle(40,15)
                shading interp
            end
            writeVideo(v, getframe(f1));
            cla(f1);%replace with delete() twice
        else
            if drawTrajectory
                j = plot3(headData(1,1:i),headData(2,1:i),headData(3,1:i),'-r');
            end
            if drawBodyCurve
                x = bodyData{:,i};
                h = streamtube({[x(1,:);x(2,:);x(3,:)]'}, 3);
                h.FaceLighting = 'gouraud';
                h.AmbientStrength = 0.8;
                h.DiffuseStrength = 0.7;
                h.SpecularStrength = 0.8;
                h.SpecularExponent = 15;
                h.BackFaceLighting = 'unlit';
                shading interp
            end
            writeVideo(v, getframe(f1));
            if drawTrajectory
                delete(j)
            end
        end
        loopEndTime = toc(loopStartTime);
        meanLoopTime = (meanLoopTime * (i-1) + loopEndTime)/(i);
        rTime = round(meanLoopTime*(framenumber-i),2);
        fprintf('%d of %d Frames (%.2f%%); estimated time remaining: %.0fh %.0fm %.0fs\n',...
            i, framenumber, round(i/framenumber*100,2), floor(rTime/3600),floor(mod(rTime/60,60)), floor(mod(rTime,60)));
    end
else
    tic
    if drawTrajectory
        plot3(headData(1,:),headData(2,:),headData(3,:),'-r');
    end
    if drawBodyCurve
        for i = 1:framenumber
            x = bodyData{:,i};
            h = streamtube({[x(1,:);x(2,:);x(3,:)]'}, 3);
            h.FaceLighting = 'gouraud';
            h.AmbientStrength = 0.8;
            h.DiffuseStrength = 0.7;
            h.SpecularStrength = 0.8;
            h.SpecularExponent = 15;
            h.BackFaceLighting = 'unlit';
        end
        shading interp
    end
    toc
end
% savefig(f1, 'fig_script_compact.fig', 'compact');
rTime = toc(startTime);
fprintf('Time elapsed: %.0fh %.0fm %.0fs %.0fms\n', ...
    floor(rTime/3600),floor(mod(rTime/60,60)), floor(mod(rTime,60)), (rTime-floor(rTime))*1000);

%Function to compute an amount of points on the circle segment described by
%the values obtained by the simulation.
function [p3,h] = line2(k, sinPhi, cosPhi, s, d, plotCurve)
r = 1/k;
pM = [r * cosPhi; r * sinPhi; d];
p1 = [0;0;d];
a = s/r;
cosa = cos(a);
sina = sin(a);
% p2 per rotation matrix
% x = p1-pM;
% n = cross(pM, [0;0;1]);
% n = -n/norm(n);
% n1 = n(1); n2 = n(2); n3 = n(3);
%R = [n1*n1*(1-cosa)+cosa n1*n2*(1-cosa)-n3*sina n1*n3*(1-cosa)+n2*sina;
%     n2*n1*(1-cosa)+n3*sina n2*n2*(1-cosa)+cosa n2*n3*(1-cosa)-n1*sina;
%      n3*n1*(1-cosa)-n2*sina n3*n2*(1-cosa)+n1*sina n3*n3*(1-cosa)+cosa];

%p2 = R*x+pM;
%or shorter
% p2 = n*dot(n,x)+cosa*cross(cross(n,x),n)+sina*cross(n,x)+pM;

p2 = p1 + (pM-p1)*(1-cosa) + [0;0;r*sina];
v1 = p1 - pM;
v2 = p2 - pM;
p3 = -cross(v2, cross(v1,v2));
p3 = p3/norm(p3)*d + p2;
h = [];
if plotCurve
    segments = 100;
    v3 = cross(cross(v1,v2),v1);
    v3 = r*v3/norm(v3);
    t = linspace(0,atan2(norm(cross(v1,v2)),dot(v1,v2)),segments);
    v = v1*cos(t)+v3*sin(t)+pM;
    h = horzcat([0; 0; 0], v, p3);
end
