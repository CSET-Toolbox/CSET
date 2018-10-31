%--------------------------------------------------------------------------
% This file is part of the ASTRA Toolbox
%
% Copyright: 2010-2014, iMinds-Vision Lab, University of Antwerp
%                 2014, CWI, Amsterdam
% License: Open Source under GPLv3
% Contact: astra@uantwerpen.be
% Website: http://sf.net/projects/astra-toolbox
%--------------------------------------------------------------------------

addpath('../');

%% Example 2: 3D parallel beam, cuda

% configuration
proj_count		= 20;
dart_iterations = 20;
outdir			= './';
prefix			= 'example2';
rho				= [0, 0.5, 1];
tau				= [0.25, 0.75];
gpu_core		= 0;

% load phantom
load('phantom3d');

% create projection and volume geometries
det_count = size(I, 1);
slice_count = size(I,3);
angles = linspace2(0, pi, proj_count);
proj_geom = astra_create_proj_geom('parallel3d', 1, 1, slice_count, det_count, angles);
vol_geom = astra_create_vol_geom(size(I));

% create sinogram
[sinogram_id, sinogram] = astra_create_sino3d_cuda(I, proj_geom, vol_geom);
astra_mex_data3d('delete', sinogram_id);

% DART
D						= DARTalgorithm(sinogram, proj_geom);
D.t0					= 100;
D.t						= 10;

D.tomography			= IterativeTomography3D();
D.tomography.method		= 'SIRT3D_CUDA';
D.tomography.gpu_core	= gpu_core;
D.tomography.use_minc	= 'yes';

D.segmentation.rho		= rho;
D.segmentation.tau		= tau;

D.smoothing.b			= 0.1;
D.smoothing.full3d		= 'yes';
D.smoothing.gpu_core	= gpu_core;
 
D.masking.random		= 0.1;
D.masking.conn			= 4;
D.masking.gpu_core		= gpu_core;

D.output.directory		= outdir;
D.output.pre 			= [prefix '_'];
D.output.save_images	= 'no';
D.output.save_results	= {'stats', 'settings', 'S', 'V'};
D.output.save_interval	= dart_iterations;
D.output.verbose		= 'yes';

D.statistics.proj_diff = 'no';

D.initialize();

D.iterate(dart_iterations);

% save the central slice of the reconstruction and the segmentation to file
imwritesc(D.S(:,:,64), [outdir '/' prefix '_S_slice_64.png']);
imwritesc(D.V(:,:,64), [outdir '/' prefix '_V_slice_64.png']);
