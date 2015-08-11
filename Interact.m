% Interact uses data from rEvaluateStructures to interactively display 3D coordinates of local aligned regions

function Comparison = Interact(Comparison)

if exist('OCTAVE_VERSION') ~= 0,
  Octave = 1;
  DotSize = 10;
  more off                           % turn off paging in Octave
  warning('off');
else
  Octave = 0;
  DotSize = 12;
end

% ---------------------------------------------- Set defaults

if ~isfield(Comparison,'AminoAcids'),
  Comparison.AminoAcids = 1;
end

AminoAcids = Comparison.AminoAcids;

% ---------------------------------------------- Set coloring of points

if strcmpi(Comparison.Coloring,'center'),
  fprintf('Coloring points by distance to the center of the molecule\n');
  zFlushOutput;
  colo = Comparison.CenterColoring;
elseif strcmpi(Comparison.Coloring,'position'),
  fprintf('Coloring points by position in the sequence\n');
  zFlushOutput;
  colo = Comparison.PositionColoring;
elseif strcmpi(Comparison.Coloring,'discrepancy'),
  fprintf('Coloring points by discrepancy between aligned local regions\n');
  zFlushOutput;
  colo = Comparison.DiscrepancyColoring;
else
  fprintf('Unrecognized coloring option, coloring by position in the sequence\n');
  zFlushOutput;
  colo = Comparison.PositionColoring;
end

% ---------------------------------------------- Plot rotation angle versus translation norm

if Comparison.RotationTranslationFigure > 0
  figure(Comparison.RotationTranslationFigure);
  clf

  scatter(Comparison.TranslationNorm,abs(Comparison.Angle),DotSize,colo,'filled');
  ax=axis;

  if Comparison.RotationTranslationColorbar > 0,
    if Octave == 1,                      % colorbar seems to confuse Octave's ability to read coordinates from the screen
    else
      colorbar('eastoutside')
    end
  end
  colormap('default');
  map = colormap;
  colormap(map(8:56,:))

  if Octave == 0,
    set(gcf,'renderer','zbuffer');
  end

  title(['Differences between ' Comparison.File1.Filename ' and ' Comparison.File2.Filename ', numbers from ' Comparison.File1.Filename]);
  hold on
  [y,k] = sort(Comparison.RotationStandardDistances,'descend');
  Numbers = '';
  I = 1;

  for j = 1:length(k),
    i = k(j);
    if Comparison.RotationStandardDistances(i) > Comparison.RotationCutoff || Comparison.TranslationStandardDistances(i) > Comparison.TranslationCutoff,
      if j > 1 && abs(diff(Comparison.RotationStandardDistances([i k(j-1)]))) < 0.001 && abs(diff(Comparison.TranslationStandardDistances([i k(j-1)]))) < 0.001,
        Numbers = [Numbers ', ' Comparison.File1.NT(Comparison.Indices1(i)).Number];
      else
        text(Comparison.TranslationNorm(I)+(ax(2)-ax(1))/100,abs(Comparison.Angle(I)),Numbers,'Fontsize',7);
        Numbers = Comparison.File1.NT(Comparison.Indices1(i)).Number;
        I = i;
      end
    end
  end
  xlabel('Translation norm');
  ylabel('Rotation angle');
  title(['Comparison of ' Comparison.File1.Filename ' and ' Comparison.File2.Filename]);
end

if isfield(Comparison,'NumNeighbors'),
  NumNeighbors = Comparison.NumNeighbors;
else
  NumNeighbors = 5;                    % number of nucleotides to show, including the current one
end

Centers1 = cat(1,Comparison.File1.NT(Comparison.Indices1).Center);      % nucleotide centers
Centers2 = cat(1,Comparison.File2.NT(Comparison.Indices2).Center);      % nucleotide centers
AACenters1 = cat(1,Comparison.File1.AA.Center);

AADist1 = zDistance(Centers1,cat(1,Comparison.File1.AA.Center));
AADist2 = zDistance(Centers2,cat(1,Comparison.File2.AA.Center));

if Comparison.RotationTranslationFigure > 0 && Comparison.Interactive > 0,
  Stop = 0;
  while Stop == 0,
    figure(Comparison.RotationTranslationFigure)

    fprintf('You can use Figure %d to select a nucleotide and its neighborhood to display.\n',Comparison.RotationTranslationFigure);
    if Octave == 1,
      fprintf('Press Enter in this window, then click a point on Figure %d.\n',Comparison.RotationTranslationFigure);
      fprintf('Right click and drag to zoom in on Figure %d.\n',Comparison.RotationTranslationFigure);
    else
      fprintf('Click a point on Figure %d and then press Enter in this window.\n',Comparison.RotationTranslationFigure);
    end
    fprintf('Or, enter a nucleotide number or range from %s; example: 12:17 or 12:17,25:31\n', Comparison.File1.Filename);
    fprintf('Press +, Enter to increase the viewing scope, -, Enter to decrease, q to quit.\n');
    fprintf('Press a to toggle the display of amino acids on or off\n');
    zFlushOutput

    k = input('','s');

    if Octave == 1,
%      fprintf('Input received, please be patient.  %s\n',k);
      zFlushOutput
    end

    if length(k) == 0,
      if Octave == 1,
        fprintf('Now click a point in Figure %d\n',Comparison.RotationTranslationFigure);
        zFlushOutput;
        [x,y,b] = ginput(1);             % apparently Octave waits now for a point to be clicked
        pt = [x y];
      else
        pt = get(gca,'CurrentPoint');    % in Matlab, the last-clicked point will be used
        x = pt(1,1);
        y = pt(1,2);
      end
      [vmin,imin] = min((Comparison.TranslationNorm - x).^2 + (abs(Comparison.Angle) - y).^2);
      i = Comparison.Indices1(imin);

      NT = Comparison.File1.NT(i);
      fprintf('Point clicked was (x,y) = (%8.4f,%8.4f), corresponding to nucleotide %s%s in chain %s\n',x,y,NT.Base,NT.Number,NT.Chain);
      if isempty(i),
        fprintf('No nucleotides selected, try again.\n');
      end

      zFlushOutput
    elseif strcmp(k,'+'),
      if NumNeighbors == 1,
        NumNeighbors = 5;
      else
        NumNeighbors = min(NumNeighbors + 5,length(Comparison.Indices1));
      end
    elseif strcmp(k,'-'),
      NumNeighbors = max(1,NumNeighbors-5);
    elseif strcmp(k,'q'),
      Stop = 1;
    elseif strcmp(k,'a'),
      AminoAcids = 1 - AminoAcids;
    elseif length(k) > 0 && ischar(k),
      i = zIndexLookup(Comparison.File1,k);
      if length(i) == 1,
        NumNeighbors = 5;
      else
        NumNeighbors = 1;
      end
      if isempty(i),
        fprintf('No nucleotides selected; perhaps the nucleotide number is not in the 3D structure file.\n')
      end
    else
      i = [];
    end

    zFlushOutput

    % find neighbors of the indicated nucleotide(s)

    if ~isempty(i) && Stop == 0,
      AllIndices1 = [];
      AllIndices2 = [];
      maxd = 0;
      for j = i,
        [a b] = sort(Comparison.MutualDistances(j,:));         % find the nearest nucleotides in File1
        AllIndices1 = [AllIndices1 Comparison.Indices1(b(1:NumNeighbors))];
        AllIndices2 = [AllIndices2 Comparison.Indices2(b(1:NumNeighbors))];
        maxd = max(maxd,a(NumNeighbors));                      % maximum distance
      end
      AllIndices1 = unique(AllIndices1);
      AllIndices2 = unique(AllIndices2);

      for w = AllIndices1,
        NT1 = Comparison.File1.NT(w);
        k = find(Comparison.Indices1 == w);
        NT2 = Comparison.File2.NT(Comparison.Indices2(k));
        fprintf('Nucleotide %s|%s|%s|%s aligns to %s|%s|%s|%s\n',Comparison.File1.Filename,NT1.Chain,NT1.Base,NT1.Number,Comparison.File2.Filename,NT2.Chain,NT2.Base,NT2.Number);
      end
      zFlushOutput

      AllAAIndices1 = [];
      AllAAIndices2 = [];
      for j = i,
        AllAAIndices1 = [AllAAIndices1 find(AADist1(j,:) < maxd)];
        AllAAIndices2 = [AllAAIndices2 find(AADist2(j,:) < maxd)];
      end
      AllAAIndices1 = unique(AllAAIndices1);
      AllAAIndices2 = unique(AllAAIndices2);

      VP2.Rotation = Comparison.OverallRotation';
      VP2.Shift = Comparison.OverallTranslation;
      VP1.Sugar = 1;
      VP2.Sugar = 1;
      VP1.Color = [1 0 0];            % red
      VP2.Color = [205 173 0]/255;    % gold
      VP1.LineThickness = 3;
      VP2.LineThickness = 3;

      figure(Comparison.CoordinateFigure)
      clf
      zDisplayNT(Comparison.File1,AllIndices1,VP1);
      zDisplayNT(Comparison.File2,AllIndices2,VP2);

      if AminoAcids > 0 && length(AllAAIndices1) > 0,
        VPA1.LabelBases = 10;
        VPA1.Color = 0.7*[1 0 0];
        VPA1.NumberColor = 0.7*[1 0 0];
        zDisplayAA(Comparison.File1,AllAAIndices1,VPA1);
      end
      if AminoAcids > 0 && length(AllAAIndices2) > 0,
        VPA2.LabelBases = 10;
        VPA2.Color = 0.7*[205 173 0]/255;    % gold
        VPA2.NumberColor = 0.7*[205 173 0]/255;    % gold
        VPA2.Rotation = Comparison.OverallRotation';
        VPA2.Shift = Comparison.OverallTranslation;
        zDisplayAA(Comparison.File2,AllAAIndices2,VPA2);
      end

      title([Comparison.File1.Filename ' in red, ' Comparison.File2.Filename ' in gold']);
      if Octave == 0,
        rotate3d on
      end

      zFlushOutput

    end
  end
end

zFlushOutput