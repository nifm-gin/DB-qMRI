path(path,'../Planeto/');
PATH_PLANETO='/local_scratch/deleforg/Data/Planeto/';

%%%%%% .lbl Synthetic Training Files %%%%
info=cub_info([PATH_PLANETO,'Ldata1.lbl']);
cube=cub_read2(info);
cubesize=size(cube.spectels);
Xtrain_synth=reshape(cube.spectels,[cubesize(1),cubesize(3)])';

info=cub_info([PATH_PLANETO,'Ldata1_noise.lbl']);
cube=cub_read2(info);
cubesize=size(cube.spectels);
Xtrain_synth_noise=reshape(cube.spectels,[cubesize(1),cubesize(3)])';

info=cub_info([PATH_PLANETO,'Ldata1_param.lbl']);
cube=cub_read2(info);
cubesize=size(cube.spectels);
Ytrain_synth=reshape(cube.spectels,[cubesize(1),cubesize(3)])';

%%%%%% .mat Synthetic Test Files %%%%%%
load([PATH_PLANETO,'XBaseAleatoire']);
Xtest_synth=XBaseAleatoire';
clear XBaseAleatoire;
load([PATH_PLANETO,'XBaseAleatoireBruit']);
Xtest_synth_noise=XBaseAleatoireBruit';
clear XBaseAleatoireBruit;
load([PATH_PLANETO,'YBaseAleatoire']);
Ytest_synth=YBaseAleatoire';
clear YBaseAleatoire;

%%%%%% .lbl Big Synthetic Test Files %%%%
info=cub_info([PATH_PLANETO,'Ldata2SGMMDT.lbl']);
cube=cub_read2(info);
cubesize=size(cube.spectels);
Xtest_synth2=reshape(cube.spectels,[cubesize(1),cubesize(3)])';

info=cub_info([PATH_PLANETO,'Ldata2SGMMDT_noise.lbl']);
cube=cub_read2(info);
cubesize=size(cube.spectels);
Xtest_synth_noise2=reshape(cube.spectels,[cubesize(1),cubesize(3)])';

info=cub_info([PATH_PLANETO,'Ldata2SGMMDT_param.lbl']);
cube=cub_read2(info);
cubesize=size(cube.spectels);
Ytest_synth2=reshape(cube.spectels,[cubesize(1),cubesize(3)])';

%%%%%% .lbl Real Training Files %%%%
info=cub_info([PATH_PLANETO,'Ldata2SGMM41.lbl']);
cube=cub_read2(info);
cubesize=size(cube.spectels);
Xtrain_real=reshape(cube.spectels,[cubesize(1),cubesize(3)])';

info=cub_info([PATH_PLANETO,'Ldata2SGMM41_noise.lbl']);
cube=cub_read2(info);
cubesize=size(cube.spectels);
Xtrain_real_noise=reshape(cube.spectels,[cubesize(1),cubesize(3)])';

info=cub_info([PATH_PLANETO,'Ldata2SGMM41_param.lbl']);
cube=cub_read2(info);
cubesize=size(cube.spectels);
Ytrain_real=reshape(cube.spectels,[cubesize(1),cubesize(3)])';

%%%%%% .lbl Real Test Files %%%%
info=cub_info([PATH_PLANETO,'im_41_GMM.lbl']);
cube=cub_read2(info);
image41=cube.spectels;

info=cub_info([PATH_PLANETO,'mask_41_GMM.lbl']);
cube=cub_read2(info);
image_mask41=cube.spectels;

info=cub_info([PATH_PLANETO,'im_61_GMM.lbl']);
cube=cub_read2(info);
image61=cube.spectels;

info=cub_info([PATH_PLANETO,'mask_61_GMM.lbl']);
cube=cub_read2(info);
image_mask61=cube.spectels;

info=cub_info([PATH_PLANETO,'im_103_GMM.lbl']);
cube=cub_read2(info);
image103=cube.spectels;

info=cub_info([PATH_PLANETO,'mask_103_GMM.lbl']);
cube=cub_read2(info);
image_mask103=cube.spectels;

images={image41,image61,image103};
image_masks={image_mask41,image_mask61,image_mask103};


clear info;
clear cube;
clear cubesize;
clear PATH_PLANETO;