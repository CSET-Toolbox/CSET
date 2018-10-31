%--------------------------------------------------------------------------
% This file is part of the ASTRA Toolbox
%
% Copyright: 2010-2014, iMinds-Vision Lab, University of Antwerp
%                 2014, CWI, Amsterdam
% License: Open Source under GPLv3
% Contact: astra@uantwerpen.be
% Website: http://sf.net/projects/astra-toolbox
%--------------------------------------------------------------------------
close all;

addpath('../');

%% Example 1: 2D parallel beam, cuda

% configuration
proj_count = 180;

% load phantom
I = phantom(256);

% create projection and volume geometries
det_count = size(I, 1);
angles = linspace2(0, pi, proj_count);
proj_geom = astra_create_proj_geom('parallel', 1, det_count, angles);
vol_geom = astra_create_vol_geom(det_count, det_count);

% create sinogram
[sinogram_id, sinogram] = astra_create_sino_cuda(I, proj_geom, vol_geom);
%astra_mex_data2d('delete', sinogram_id);


%show generated sinogram
figure, imagesc(sinogram), colormap('gray');

