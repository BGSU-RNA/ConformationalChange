% This file specifies two RNA 3D structure files, calculates conformational changes, and interactively displays local regions
clear Comparison
File1 = '2UUB';
Indices1 = 1:1512;      % all nucleotides in chain A
File2 = '2UUC';
Indices2 = 1:1512;      % all nucleotides in chain A

Comparison.Coloring = 'center';    % 'position' for position in the sequence, 'center' for distance from the center of the molecule, or 'discrepancy'
Comparison.RotationFigure = 1;       % figure number for rotation data points, or 0 for no figure
Comparison.TranslationFigure = 2;    % figure number for translation data points, or 0 for no figure
Comparison.RotationTranslationFigure = 3;  % figure number for angle of rotation versus norm of translation
Comparison.CoordinateFigure = 4;     % figure to show 3D coordinates in
Comparison.Interactive = 1;          % after plotting, allow the user to click on the RotationTranslationFigure to select a nucleotide and its neighborhood
Comparison.AminoAcids= 0;            % do not show amino acids when displaying local superpositions

Comparison.RotationData = 1;         % 1 means to return data in Comparison.Rcell, 2 means to also write to Excel
Comparison.TranslationData = 1;      % 1 means to return data in Comparison.Tcell, 2 means to also write to Excel
Comparison.RotationCutoff = 15;      % label nucleotides whose rotation vector has standardized distance above this number
Comparison.TranslationCutoff = 15;   % label nucleotides whose translation vector has standardized distance above this number
Comparison.UseR3DAlign = 0;
Comparison.UseR3DandNW = 0;

% if Matlab or Octave is in the FR3D directory, the following commands will add useful directories to the path

warning off
addpath('FR3DSource');
addpath('PrecomputedData');
addpath('SearchSaveFiles');
addpath('PDBFiles');

% evaluate the structures for conformational change

[Comparison] = CompareStructures(File1,Indices1,File2,Indices2,Comparison);

save(['Comparison_2UUB_2UUC.mat'],'Comparison');

% the following line shows how to load a previously calculated comparison between structures to avoid re-calculating
% load(['Comparison_2UUB_2UUC.mat']);

Comparison = Interact(Comparison);

% run just the line below to re-do the graphs without additional computation.

% [Comparison] = CompareStructures(Comparison.File1,Comparison.Indices1,Comparison.File2,Comparison.Indices2,Comparison);

