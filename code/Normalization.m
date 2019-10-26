clear all;
close all;
clc;
%%
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
    
