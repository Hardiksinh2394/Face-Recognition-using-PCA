    clear all; close all; clc;
%% Training data set part 
Train_images_path = '.\Normalized images\train_images';
Train_images = fullfile(Train_images_path,'*.jpg');
Load_train_images = dir(Train_images);

D_matrix = []; %empty d matrix 
Label_matrix = []; %Lt matrix for assiging the person name 
Dimension_vector = 4096; % 64*64 vector of dimension 4096
for i = 1:length(Load_train_images)
    image = imread(cat(2,Load_train_images(i).folder,'\',Load_train_images(i).name));
    image_vector = reshape(image,1,Dimension_vector); % Resized the images in 64*64=4096 vector 
    %computing eigenfaces from training set
    D_matrix = [D_matrix;image_vector]; % we get D_matrix(p*d) where each row of d corspond to training image 
    Label_matrix = [Label_matrix;Load_train_images(i).name(1:3)]; %assige label and give a name of the person in Lt
end
%% perform PCA of D 
% Calculating mean of D_matrix 
mean_D = mean(D_matrix); %get mean of D_matrix 
mean_D = floor(mean_D); %get round matrix 
D_matrix = double(D_matrix); 


% get covariance matrix from removing mean
Cov_matrix = D_matrix-mean_D;
Covariance = 1/(size(Cov_matrix,1)-1)*Cov_matrix*Cov_matrix'; % Covariance matrix computation
[eigvect,eigval] = eig(Covariance); % Eigenvectors and eigenvalues of the convariance matrix

%% Finding PCA of D
PHI = Cov_matrix'*eigvect; %transformation matrix PHI
Ft = []; % create empty Feature vector matrix
for j = 1:size(D_matrix,1) 
    Ft(j,:) = D_matrix(j,:)*PHI; % get Feature vector of each image and store in Ft
end

%% Test data set part 
test_image_path = '.\Normalized images\test_images';
Test_images = fullfile(test_image_path,'*.jpg');
Load_Test_images = dir(Test_images);
% Loading an image from the test set
Test_label_matrix = []; %assigh the test label matrix for person name 
Error = 0;

for t = 1:size(Load_Test_images,1) 
image_test = imread(cat(2,Load_Test_images(t).folder,'\',Load_Test_images(t).name));
image_test = double(image_test);
figure, imshow(image_test,[])
%assign the label and give name of person 
Test_label_matrix = [Test_label_matrix;Load_Test_images(t).name(1:3)]; 

% Reshaping the test image in 1*4096 
test_img_vector = reshape(image_test,1,Dimension_vector); 
phi_test = test_img_vector*PHI; % Projecting the test image into PCA space
Distance = pdist2(Ft,phi_test); % find a distance between phi_test(phi_q) and Ft(feature vector)
[Min index] = sort(Distance); % get sort distance between ft and phi_test

k=10; %number of faces in which you want to look
for number_faces = 1:k %% loop for 1 to k faces 
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

close all;
end
%% correct recognition accuracy measurement
accuracy = (1-(Error/size(Load_Test_images,1)))*100;