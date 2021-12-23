%%Demonstration Code for Data Visualization Seminar Offered by the Chair of Computer Graphics and Visualization at TUM
%%Seminor Topic: Glyph-based Tensor Field Visualization
%%Author: Junpeng Wang (junpeng.wang@tum.de)
%%Date: 2021-09-14
clear
clc

stressfileName = './data/cantilever2D_R500_iLoad5.carti';
samplingSpace = 20; %%control the density of the tensor glyphs (negatively-related) (evenly-spaced seeding per-elements)
scalingGlyph = 1.0; %%control the size of tensor glyphs (positively-related)

EllipticTensorGlyph2D(stressfileName,samplingSpace,scalingGlyph);
