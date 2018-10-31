% -----------------------------------------------------------------------
% SIRT reconstruction algorithm using the projection and backprojection <
% operations defined in ASTRA toolbox
% -----------------------------------------------------------------------

close all;
clear;
clc;

tic;

%Add folders to MATLAB path
% addpath(genpath('C:/Users/banjak/Desktop/Post-doc/1.Programs/Reconstruction/Reconstruction-SART-TV-FISTA/Data'));
% addpath(genpath('C:/Users/banjak/Desktop/Post-doc/1.Programs/Reconstruction/Reconstruction-SART-TV-FISTA/ASTRA-toolbox'));
% addpath(genpath('C:/Users/banjak/Desktop/Post-doc/1.Programs/Reconstruction/Reconstruction-SART-TV-FISTA/denoising_FISTA-TV'));
% addpath(genpath('C:/Users/banjak/Desktop/Post-doc/1.Programs/Reconstruction/Reconstruction-SART-TV-FISTA/Matlab-toolbox'));
% addpath(genpath('C:/Users/banjak/Desktop/Post-doc/1.Programs/Reconstruction/Reconstruction-SART-TV-FISTA/Spot-toolbox'));

% Create 3D volume geometry
nb_rows = 256;
nb_columns = 256;
nb_slices = 256;

vol_geom = astra_create_vol_geom(nb_rows, nb_columns, nb_slices);

nb_projs = 69;

%% Create projection geometry
N = 256; % nb. of columns (Width in ImageJ)
pix_size_x = 1; %(value=(rec resolution)/(projection resolution))

M = 256; % nb. of rows (Height in imageJ)
pix_size_y = 1;

fid_angle =fopen('angles_PdSiO2.rawtlt','rt');

for i=1:nb_projs
angles(1,i) = str2double(fgetl(fid_angle));
end

angles = angles*pi/180;

proj_geom = astra_create_proj_geom('parallel3d', pix_size_x, pix_size_y, M, N, angles);

Voxel = [nb_rows, nb_columns, nb_slices];

%%
%-------------------------------------------------------------------
%---------------Reconstruction algorithm and parameters-------------
%-------------------------------------------------------------------
regularization = 1;
FISTA = 1;
lambda = 1; % SIRT relaxation factor
Niter = 150; % nb. of iterations
Algo_type = 'SIRT'; % options: SART, SIRT, OS-SART
Nsubsets = 9; % Number of subsets (considered only if OS-SART was chosen)
alpha_TV = 1; % relaxation parameter for TV minimization
Ntv = 10; % nb. of iterations for TV minimization
t_acc = 1; % parameter for FISTA acceleration technique

%% Generate projection data
% Create the Spot operator for ASTRA using the GPU.
W = opTomo('cuda', proj_geom, vol_geom);
nb_pixels_sino = nb_projs*M*N;
%--------------------------------------------------------------------
%--------------------------------------------------------------------

% Load saved sinogram
file = fopen('PdSiO2_aligned_tilt_series_256_256_69_16int.raw');
sino_import = fread(file,'int16','ieee-be');

% file = fopen('PdSiO2_aligned_tilt_series_1024_1024_69_16int.raw');
% sino_import = fread(file,'int16');

sino_import = reshape(sino_import,[M N nb_projs]);
p = zeros(M,nb_projs,N);

for i=1:nb_projs    
slice_import = sino_import(:,:,i);
p(:,i,:) = slice_import;
end
name_sino = sprintf('./Reconstruction_results/Sino_%d_%d_%d.raw',M,N,nb_projs);
file_sino = fopen(name_sino,'w');
fwrite(file_sino,sino_import,'float');
fclose(file_sino);
p = reshape(p,[nb_pixels_sino 1]);

%% Reconstruction 
nb_voxels = nb_rows*nb_columns*nb_slices;

if(strcmp(Algo_type,'SART'))
    Nsubsets = nb_projs;
elseif(strcmp(Algo_type,'SIRT'))
    Nsubsets = 1;
else
    SubdivisionCheck(nb_projs, Nsubsets);
end

% partition of the projection indices
Ns = nb_projs/Nsubsets; % nb. of projections for each subset
NraysSub = Ns*N*M; % nb. of rays for each subset
Wsubsets = cell(1,Nsubsets);
for k=1:Nsubsets
    Wsubsets{k} = ((k-1)*NraysSub+1):(k*NraysSub);
end

% Initialize image by zeros
rec = zeros(nb_rows,nb_columns,nb_slices);

% Transform into vector form
y = rec(:); 
y_temp2 = rec(:);

Rsubs = cell(1,Nsubsets);
Csubs = cell(1,Nsubsets);

fprintf('Reconstruction using %s algorithm \n',Algo_type)
if(regularization && FISTA)
    disp('TV regularization and FISTA acceleration are applied');
elseif(regularization)
    disp('TV regularization is applied');
end


%% Iterate
for i =1:Niter
    
    fprintf('Iteration: %d \n',i);
    
   y_prev = y;
    
    for k=1:Nsubsets
        
        % Create weighting matrices R and C during 1st iteration only
        if(i==1)
            Rsubs{k} = 1 ./ (W(Wsubsets{k},:)*ones(size(y)));
            Csubs{k} = 1 ./ (W(Wsubsets{k},:)'*ones(NraysSub,1));
            Rsubs{k}(Rsubs{k}==Inf) = 0;
            Csubs{k}(Csubs{k}==Inf) = 0;
        end
        
        %compute projection difference
        u = p(Wsubsets{k}) - W(Wsubsets{k},:)*y;
        
        %update image y
        y = y + lambda*Csubs{k}.*(W(Wsubsets{k},:)'*(Rsubs{k}.*u));
            
        % non-negative regularization
        y(y<0) = 0; 
    end
        
       %=================================================
       %-------------- TV regularization ----------------
       %=================================================
      if(regularization==1)
          disp('TV regularization');
%           % reshape y into a 3D matrix
%           y = reshape(y, W.vol_size);
%           y_prev = reshape(y_prev, W.vol_size);
%           d=im3Dnorm(y-y_prev,'L2');
%           y = TVkernel(y, alpha_TV, d, Ntv);
%           % reshape y into a 1D array
%           y = reshape(y, [nb_voxels 1]);
           
       %=================================================
       %---------- CUDA Implementation  -----------------
       %=================================================
       y = reshape(y, W.vol_size);  
       y=minimizeTV(y,alpha_TV,Ntv);
       y = reshape(y, [nb_voxels 1]); 
      end
       %=================================================
       
       %=================================================
       %------------- FISTA acceleration ----------------
       %=================================================
       if(regularization && FISTA)
            disp('FISTA');
			t_acc_inc = (1 + sqrt(1 + 4 *t_acc.^2))/2;
            y_temp1 = y;
            y = y + ((t_acc-1)/t_acc_inc)*(y-y_temp2);
            t_acc = t_acc_inc;
            y_temp2 = y_temp1;
       end
     
    %======================================================    
    %-----Save reconstructed volume at curent iteration----
    %======================================================
    if(~mod(i,5))
        reconstruction_sirt = reshape(y, W.vol_size); % Transform into 3D matrix form
        name_img = sprintf('./Reconstruction_results/SIRT_rec_%d_%d_%d_Niter=%d.raw',nb_columns,nb_rows,nb_slices,i);
        file_img = fopen(name_img,'w');
        fwrite(file_img,reconstruction_sirt,'double');
        fclose(file_img);    
    end
    
end
%end of iteration
%%

% Show 3D reconstructed image
% figure, imshow3D(reconstruction_sirt);
% title('SIRT Reconstruction');

toc;

%%
% Function: check if the given number of susbsets can divide the total number of projections
function SubdivisionCheck(nb_projs, Nsubsets)
    if(mod(nb_projs,Nsubsets) ~= 0)
        error('Nsubsets must divide the total nb. of projections');
    end
end

%%
