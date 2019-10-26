function varargout = Face_recognition(varargin)
% FACE_RECOGNITION MATLAB code for Face_recognition.fig
%      FACE_RECOGNITION, by itself, creates a new FACE_RECOGNITION or raises the existing
%      singleton*.
%
%      H = FACE_RECOGNITION returns the handle to a new FACE_RECOGNITION or the handle to
%      the existing singleton*.
%
%      FACE_RECOGNITION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FACE_RECOGNITION.M with the given input arguments.
%
%      FACE_RECOGNITION('Property','Value',...) creates a new FACE_RECOGNITION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Face_recognition_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Face_recognition_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Face_recognition

% Last Modified by GUIDE v2.5 11-Jan-2019 23:17:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Face_recognition_OpeningFcn, ...
                   'gui_OutputFcn',  @Face_recognition_OutputFcn, ...
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


% --- Executes just before Face_recognition is made visible.
function Face_recognition_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Face_recognition (see VARARGIN)

% Choose default command line output for Face_recognition
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.text_state,'String',' '); 
set(handles.textaccu,'String',''); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UIWAIT makes Face_recognition wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Face_recognition_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index = get(hObject,'Value');
test_file = handles.test_file;
set(handles.test_name,'String',test_file(index).name);
imagetest = imread(cat(2,test_file(index).folder,'\',test_file(index).name));
imshow(imagetest,'Parent',handles.axes1);

handles.index = index;
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k = handles.k;  
    set(handles.text_state,'String','Loading.....');
    %files_path = uigetdir(' '); %Getting the directory where the images are.
    if k==0 | not(isnumeric(k))
        k=1;
    end
    
    Train_images_path = '.\Normalized images\train_images';
    Train_images = fullfile(Train_images_path,'*.jpg');
    Load_train_images = dir(Train_images);

    D_matrix = [];     %empty d matrix
    Label_matrix = [];  %Label matrix for assiging the person name
    Dimension_vector = 4096;    % 64*64 vector of dimension 4096
    for i = 1:length(Load_train_images)
        image = imread(cat(2,Load_train_images(i).folder,'\',Load_train_images(i).name));
        image_vector = reshape(image,1,Dimension_vector);  % Resized the images in 64*64=4096 vector 
        %computing eigenfaces from training set
        D_matrix = [D_matrix;image_vector];
        Label_matrix = [Label_matrix;Load_train_images(i).name(1:3)];   
    end
    %% perform PCA of D 
    % Calculating mean of D_matrix
    mean_D = mean(D_matrix);
    mean_D = floor(mean_D);
    D_matrix = double(D_matrix);
   
    % get covariance matrix from removing mean
    Cov_matrix = D_matrix-mean_D;
    Covariance = 1/(size(Cov_matrix,1)-1)*Cov_matrix*Cov_matrix'; % Covariance matrix computation
    [eigvect,eigval] = eig(Covariance); % Eigenvectors and eigenvalues of the convariance matrix

    %% Finding PCA of D
    PHI = Cov_matrix'*eigvect; %transformation matrix PHI
    Ft = []; % Feature vector matrix
    for j = 1:size(D_matrix,1)
        Ft(j,:) = D_matrix(j,:)*PHI; % get Feature vector of each image and store in Ft
    end
    %% Test data set part 
    test_image_path = '.\Normalized images\test_images';
    Test_images = fullfile(test_image_path,'*.jpg');
    Load_Test_images = dir(Test_images);
    % Loading an image from the test set
    Test_label_matrix = []; %%assigh the test label matrix for person name 
    Error = 0;
    for t = 1:size(Load_Test_images,1)
    image_test = imread(cat(2,Load_Test_images(t).folder,'\',Load_Test_images(t).name));
    image_test = double(image_test);
    Test_label_matrix = [Test_label_matrix;Load_Test_images(t).name(1:8)];

    % Reshaping the test image in 1*4096 
    test_img_vector = reshape(image_test,1,Dimension_vector);
    phi_test = test_img_vector*PHI; % Projecting the test image into PCA space
    Distance = pdist2(Ft,phi_test); % find a distance between phi_test(phi_q) and Ft(feature vector)
    [Min index] = sort(Distance);   % get sort distance between ft and phi_test
    
    for number_faces = 1:k % loop for 1 to k faces
    Iz = D_matrix(index(number_faces),:); % Most similar k face
    if strcmp(Label_matrix(index(number_faces),:),Load_Test_images(t).name(1:3))
        break;
    end
    if number_faces==k
        Error = Error+1;
        number_faces=1;
        break;
    end
    end
    
    end
    %% correct recognition accuracy measurement
    accuracy = (1-(Error/size(Load_Test_images,1)))*100;
    accuracy = floor(accuracy);
    percentage = num2str(accuracy);
    accuracy2 = cat(2,percentage,'%');
    
    handles.test_file = Load_Test_images;
    handles.Phi = PHI;
    handles.Ft = Ft;
    handles.D = D_matrix;
    handles.train_file = Load_train_images;
    handles.Lt = Label_matrix;
    guidata(hObject,handles);
    
   
    set(handles.listbox1,'String',Test_label_matrix);
    set(handles.text_state,'String','Done!');
    set(handles.textaccu,'String',num2str(accuracy2));

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles) %%%%%%%%%%%%%%
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index = handles.index;
test_file = handles.test_file; 
PHI = handles.Phi;
Ft = handles.Ft;
train_file = handles.train_file;
D_matrix = handles.D;
k = handles.k;
Label_matrix = handles.Lt;

    if k==0 | not(isnumeric(k))
        k=1;
    end

    itest = imread(cat(2,test_file(index).folder,'\',test_file(index).name));
    itest = double(itest);
    
    % Reshaping the test image
    test_img_vector = reshape(itest,1,64*64);
    phi_test = test_img_vector*PHI; % Projecting the test image into PCA space

    Distance = pdist2(Ft,phi_test);
    [Min position] = sort(Distance); % finding the k's minimum distances

    for number_faces = 1:k 
    Iz = D_matrix(position(number_faces),:); % Most similar k face
    

    if strcmp(Label_matrix(position(number_faces),:),test_file(index).name(1:3))
        break;
    end
    if number_faces==k
        number_faces=1;
        break;
    end
    end
    
    
%     [Min position] = min(Distance); % finding the minimum distance
    Iz = D_matrix(position(number_faces),:); % Most similar face
    Iz = reshape(Iz,64,64);
    
    set (handles.text4,'String',train_file(position(number_faces)).name);
    imshow(Iz,[],'Parent', handles.axes2);


function k_faces_Callback(hObject, eventdata, handles)
% hObject    handle to k_faces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k = uint8(str2double(get(hObject,'String')));

handles.k = k;
guidata(hObject,handles);

% Hints: get(hObject,'String') returns contents of k_faces as text
%        str2double(get(hObject,'String')) returns contents of k_faces as a double


% --- Executes during object creation, after setting all properties.
function k_faces_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_faces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
k = uint8(str2double(get(hObject,'String')));

handles.k = k;
guidata(hObject,handles);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
images = dir('Main_Image\*.jpg'); %read image files 
images_features = dir('text_file\*.txt'); %read text files 
P_features = [13; 20; 50; 20; 34; 34; 16; 50; 48; 50]; %predefiend feature points 
Error  = true; 
iterations = 0;
N = 64; % define N for 64 * 64 images 
feature_matrix = []; %create empty matrix 
%% for single image
text_path = fullfile(images_features(1).folder,images_features(1).name); %load text path for the first image
F1 = load(text_path);
F1_bar = F1'; %get f_bar
F_vector = F1_bar(1:end); %get f_vector 


while Error 
    %P_features = P_features';
    [A,B] = FindTransformation(P_features,F_vector); % get AB transformation from defind features and f vector 
    F_bar_transpose = ApplyTransformation(A,B,F_vector); % to get f bar transpose 
    
   %% for all images
   
    for a = 1 : length(images_features)
        Features_path = fullfile(images_features(a).folder,images_features(a).name); % load text path for all images 
        Fp = load(Features_path); 
        Fp = Fp';  %    transpose of all features 
        Fp = Fp(1:end); 
        F2_bar = F_bar_transpose';
        [A1,B1] = FindTransformation(F2_bar,Fp); %get A1 and B1 
        f1_bar_prime = ApplyTransformation(A1,B1,Fp); % apply a1 and b1 trasformation to all features transpose 
        feature_matrix = [feature_matrix ; f1_bar_prime]; % create the matrix of feature matrix and f tanspose i
    end
    F_mean = mean(feature_matrix); % get avg value of feature matrix 
    
    if immse(F_mean,F_vector) <= 0.00000001 % compare means feature matrix with f vector 
        Error = false;
        break;
    end
    F_vector = F_mean;
    iterations = iterations + 1;
end

for a = 1:length(images_features) 
    Features_path = fullfile(images_features(a).folder,images_features(a).name);
    Fp = load(Features_path);
    Fp = Fp';
    Fp = Fp(1:end);
    F2_bar = F_vector';
    [A,B] = FindTransformation(F2_bar,Fp);
    
    normalization = zeros(64,64); % conver in 64 *64 image 
    normalization_i = imread(cat(2,images(a).folder,'\',images(a).name));
    image_irgb = rgb2gray(normalization_i); % convert rgb to gray image 
    mkdir('Normalized images'); % create the new directory 
    
    for i = 1:N 
        for j = 1:N
            xy_origin = inv(A)*([i;j]-B);
            xy_origin = ceil(xy_origin);
            if xy_origin(1,1)<=0
                xy_origin(1,1)=1;
            end
            if xy_origin(1,1)>240
                xy_origin(1,1)=240;
            end
            if xy_origin(2,1)<=0
                xy_origin(2,1)=1;
            end
            if xy_origin(2,1)>320
                xy_origin(2,1)=320;
            end
            normalization(j,i) = image_irgb(xy_origin(2),xy_origin(1));
        end
    end
    normalization_point = uint8(normalization);
    imwrite(normalization_point,fullfile('normalized images',images(a).name));
end
