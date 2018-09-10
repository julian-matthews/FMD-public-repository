function [ gabortex, gaborprops] = make_gabor( Cfg , texsize, contrast, phase )
%% MAKE_GABOR [gabor_texture]=make_gabor(Cfg,texture_size,cycles)
% Makes gabor-texture grating and transparency mask
% Input Cfg details, size of required texture

% Dimension of the region where will draw the Gabor in pixels
gaborDimPix = texsize;

% Sigma of Gaussian
sigma = Cfg.pixelsPerDegree / 2;

% Obvious Parameters
aspectRatio = 1.0;

% Spatial Frequency (Cycles Per Pixel)
% One Cycle = Grey-Black-Grey-White-Grey i.e. One Black and One White Lobe
numCycles = 2; % As per Sherman et al. (2015) - 2 cycles per degree
freq = numCycles / Cfg.pixelsPerDegree; % 

% Build a procedural gabor texture (Note: to get a "standard" Gabor patch
% we set a grey background offset, disable normalisation, and set a
% pre-contrast multiplier of 0.5.
% For full details see:
% https://groups.yahoo.com/neo/groups/psychtoolbox/conversations/topics/9174
backgroundOffset = [0.5 0.5 0.5 0.0];
disableNorm = 1;
preContrastMultiplier = 0.5;
gabortex = CreateProceduralGabor(Cfg.windowPtr, gaborDimPix, gaborDimPix, [],...
    backgroundOffset, disableNorm, preContrastMultiplier);

% Randomise the phase of the Gabors and make a properties matrix.
gaborprops = [phase, freq, sigma, contrast, aspectRatio, 0, 0, 0];

end

