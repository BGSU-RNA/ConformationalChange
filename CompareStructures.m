% rEvaluateStructures compares two RNA-containing 3D structures, indicating local translational and rotational differences
% File1 is a PDB ID or FR3D data structure for the first RNA 3D structure
% Indices1 is a text description of the indices from File1 to use, or a vector of FR3D indices
% File2 is a PDB ID or FR3D data structure for the second RNA 3D structure
% Indices2 is a text description of the indices from File1 to use, or a vector of FR3D indices
% Comparison is a structured variable used to set user choices
% Comparison.Coloring can be 'position' (position in chain), 'center' (distance from center of the molecule), 'discrepancy' (local discrepancy)
% Comparison.RotationFigure is an integer telling in what figure window to show local rotational changes, 0 for no figure
% Comparison.RotationColorbar is 1 to show a colorbar in RotationFigure or 0 for no bar
% Comparison.TranslationFigure is an integer telling in what figure window to show local translational changes, 0 for no figure
% Comparison.TranslationColorbar is 1 to show a colorbar in TranslationFigure or 0 for no bar
% Comparison.RotationTranslationFigure is an integer telling in what figure window to show rotation angle and translation vector norm, 0 for no figure
% Comparison.RotationTranslationColorbar is 1 to show a colorbar in TranslationFigure or 0 for no bar
% Comparison.RotationData is 1 to to calculate and return additional data about rotations, returned in Comparison.Rcell; 0 for no calculations; this feature requires the statistics toolbox
% Comparison.TranslationData is 1 to to calculate and return additional data about translations, returned in Comparison.Tcell; 0 for no calculations; this feature requires the statistics toolbox
% Comparison.RotationCutoff tells the standardized rotation score above which nucleotide numbers should be shown
% Comparison.TranslationCutoff tells the standardized translation score above which nucleotide numbers should be shown
% Comparison.Interactive is 1 to allow clicking on the RotationTranslationFigure to see 3D coordinates of aligned nucleotides
% Comparison.CoordinateFigure tells in what figure window to display 3D coordinates
% Comparison.AminoAcids tells whether or not to display nearby amino acids in CoordinateFigure

% Returned variables:
% Comparison.RotationCovariance is the variance-covariance matrix using rotation data
% Comparison.TranslationCovariance is the variance-covariance matrix using translation data

function [Comparison] = rEvaluateStructures(File1,Indices1,File2,Indices2,Comparison)

if exist('OCTAVE_VERSION') ~= 0,
  Octave = 1;
  DotSize = 8;
  more off                                       % turn off paging
  warning('off');
else
  Octave = 0;
  DotSize = 12;
end

% set default values

if ~isfield(Comparison,'Coloring'),
  Comparison.Coloring = 'position';
end
if ~isfield(Comparison,'RotationFigure'),
  Comparison.RotationFigure = 1;
end
if ~isfield(Comparison,'RotationColorbar'),
  Comparison.RotationColorbar = 1;
end
if ~isfield(Comparison,'TranslationFigure'),
  Comparison.TranslationFigure = 2;
end
if ~isfield(Comparison,'TranslationColorbar'),
  Comparison.TranslationColorbar = 1;
end
if ~isfield(Comparison,'RotationTranslationFigure'),
  Comparison.RotationTranslationFigure = 3;
end
if ~isfield(Comparison,'RotationTranslationColorbar'),
  Comparison.RotationTranslationColorbar = 1;
end
if ~isfield(Comparison,'RotationData'),
  Comparison.RotationData = 1;
end
if ~isfield(Comparison,'TranslationData'),
  Comparison.TranslationData = 1;
end
if ~isfield(Comparison,'RotationCutoff'),
  Comparison.RotationCutoff = 15;
end
if ~isfield(Comparison,'TranslationCutoff'),
  Comparison.TranslationCutoff = 15;
end
if ~isfield(Comparison,'Interactive'),
  Comparison.Interactive = 1;
end
if ~isfield(Comparison,'CoordinateFigure'),
  Comparison.CoordinateFigure = 1;
end

if Comparison.RotationTranslationFigure == 0 && Comparison.Interactive > 0,
  fprintf('rEvaluateStructures:  Set Comparison.RotationTranslationFigure to enable interactive mode\n');
end
if Comparison.CoordinateFigure == 0 && Comparison.Interactive > 0,
  fprintf('rEvaluateStructures:  Set Comparison.CoordinateFigure to enable interactive mode\n');
end

% ------------------------------ Load File1 and File2 and identify desired nucleotides within them

if ischar(File1),
  Filename = File1;
  File1 = zAddNTData(Filename,0);
end

if ischar(Indices1),
  NTList1 = Indices1;
  Indices1 = zIndexLookup(File1,Indices1);
  fprintf('Selected %d nucleotides from %s\n',length(Indices1),File1.Filename);
  zFlushOutput;
else
  NTList1 = [];
end

if ischar(File2),
  Filename = File2;
  File2 = zAddNTData(Filename,0);
end

if ischar(Indices2),
  NTList2 = Indices2;
  Indices2 = zIndexLookup(File2,Indices2);
  fprintf('Selected %d nucleotides from %s\n',length(Indices2),File2.Filename);
  zFlushOutput;
else
  NTList2 = [];
end

if ~ischar(NTList1) && ~ischar(NTList2) && length(Indices1) ~= length(Indices2),
  fprintf('Specified index lists do not have the same lengths\n');
  zFlushOutput;
end

if ischar(NTList1) || ischar(NTList2) || length(Indices1) ~= length(Indices2),
  fprintf('Using Needleman-Wunsch to align nucleotide sequences\n');
  zFlushOutput;
  Sequence1 = cat(2,File1.NT(Indices1).Base);
  Sequence2 = cat(2,File2.NT(Indices2).Base);
  [matches,align1,align2,s1,s2] = zNeedlemanWunschAffineGap(Sequence1,Sequence2);
  fprintf('Aligned %d nucleotides\n',length(align1));
  zFlushOutput;
  Indices1 = Indices1(align1);
  Indices2 = Indices2(align2);
end

% ------------------------------------ Superimpose local neighborhoods

Centers1 = cat(1,File1.NT(Indices1).Center);      % nucleotide centers
Centers2 = cat(1,File2.NT(Indices2).Center);      % nucleotide centers

DM = full(zMutualDistance(Centers1,Inf));         % mutual distances between nucleotide centers in File1

n = length(Indices1);
RR = zBestRotation(Centers2-ones(n,1)*mean(Centers2),Centers1-ones(n,1)*mean(Centers1));
tt = mean(Centers2) - mean(Centers1)*RR;          % shift needed to move points from File2 to File1 when plotting

TransformedCenters2 = (Centers2 - (ones(length(Indices1),1) * mean(Centers2) - ones(length(Indices1),1) * mean(Centers1) * RR)) * RR';

for i=1:length(Indices1)                          % loop through nucleotides in File1
  [a b]=sort(DM(i,:));                            % find the nearest nucleotides in File1
  [Discrepancy(i), Rot] = xDiscrepancy(File1,Indices1(b(1:5)),File2,Indices2(b(1:5))); % local discrepancy and rotation

  R{i} = (RR*Rot)';                               % calculate rotation matrix relative to overall rotation

  [AX,ANG] = zAxisAngle(R{i});                    % relative rotation axis and angle

  Xr(i)=AX(1)*ANG;                                % three components of the axis-angle vector, for plotting
  Yr(i)=AX(2)*ANG;
  Zr(i)=AX(3)*ANG;

  Comparison.Angle(i) = abs(ANG);                 % store to return
  Comparison.Axis{i} = AX;

  t{i} = mean(TransformedCenters2(b(1:5),:)) - mean(Centers1(b(1:5),:)); % relative translation

  Comparison.TranslationNorm(i) = norm(t{i});     % store to return

  Xt(i)=t{i}(1);                                  % three components of the translation, for plotting
  Yt(i)=t{i}(2);
  Zt(i)=t{i}(3);

  CenterDistance(i) = norm(mean(Centers1(Indices1(b(1:5))) - mean(Centers1(Indices1))));
end

% ---------------------------------------------- Set coloring of points

Comparison.CenterColoring = CenterDistance;
Comparison.PositionColoring = 1:n;
Comparison.DiscrepancyColoring = Discrepancy;

if strcmpi(Comparison.Coloring,'center'),
  ColorText =  'Coloring points by distance to the center of the molecule';
  colo = Comparison.CenterColoring;
elseif strcmpi(Comparison.Coloring,'position'),
  ColorText =  'Coloring points by position in the sequence';
  colo = Comparison.PositionColoring;
elseif strcmpi(Comparison.Coloring,'discrepancy'),
  ColorText =  'Coloring points by discrepancy between aligned local regions';
  colo = Comparison.DiscrepancyColoring;
else
  ColorText =  'Unrecognized coloring option, coloring by position in the sequence';
  colo = Comparison.PositionColoring;
end

% ---------------------------------------------- Plot rotation vectors

if Comparison.RotationFigure > 0
  fprintf('%s\n',ColorText);
  zFlushOutput;
  figure(ceil(Comparison.RotationFigure));
  clf

  scatter3(Xr,Yr,Zr,DotSize,colo,'filled');
  if Comparison.RotationColorbar > 0,
    colorbar('eastoutside')
  end
  colormap('default');
  map = colormap;
  colormap(map(8:56,:))

  if Octave == 0,
    set(gcf,'renderer','zbuffer');
  end

  title(['Rotation differences between ' File1.Filename ' and ' File2.Filename ', numbers from ' File1.Filename]);
  if Octave == 0,
    rotate3d on
  end
  hold on

  ax = axis;
end

n=length(Xr);           %number of observation
Rs=[Xr;Yr;Zr];          %3xn matrix of all observations
XrBar=mean(Xr);         %scalar
YrBar=mean(Yr);         %scalar
ZrBar=mean(Zr);         %scalar
XrSD=std(Xr);           %scalar
YrSD=std(Yr);           %scalar
ZrSD=std(Zr);           %scalar
RBar=[XrBar;YrBar;ZrBar];            %3x1 vector
Rm=RBar * ones(1,length(Rs));        %3xn matrix where each column is the value of RBar
Sr=(Rs-Rm)*(Rs-Rm)'/length(Rs);   %3x3 covariance matrix
SrINV=inv(Sr);                    %3x3 inverse covariance matrix
CorrR=inv(sqrt(diag(diag(Sr))))*Sr*inv(sqrt(diag(diag(Sr))));
RD=Rs-Rm;                        %3xn matrix-subtract RBar from each entry
DSQ=RD'*SrINV*RD;                %nxn where diagonal elements are the standardized squared distance...
D=diag(diag(DSQ));               %nxn diagonal matrix with di-sq's on diagonal
D2r=D*ones(n,1);                  %nx1 vector of di-sq's
[a b]=sort(D2r,'descend');

Rcell={};
Numbers=' ';
I=1;
Rcell{1,2}='D2';
for j=1:n
   Rcell{j+1,2}=a(j);
   Rcell{j+1,3}=[File1.NT(Indices1(b(j))).Base File1.NT(Indices1(b(j))).Number];
   if Comparison.RotationFigure > 0
      if D2r(b(j)) > Comparison.RotationCutoff,
         if j>1 && abs(D2r(b(j))-D2r(b(j-1)))<.0001
            Numbers=[Numbers ', ' File1.NT(Indices1(b(j))).Number];
         else
            text(Xr(b(I))+(ax(2)-ax(1))/100,Yr(b(I))+(ax(4)-ax(3))/100,Zr(b(I)),Numbers,'Fontsize',7);
            Numbers=File1.NT(Indices1(b(j))).Number;
            I=j;
         end
      end
   end
end

zFlushOutput;

if Comparison.RotationData > 0,

  if exist('mvncdf') && exist('normcdf'),
    y = mvncdf(Rs',RBar',Sr);
    z = y;
    z(z>.5)=1-z(z>.5);
    [a b]=sort(z,'ascend');
    for j=1:n
       Rcell{j+1,5}=y(b(j));
       Rcell{j+1,6}=File1.NT(Indices1(b(j))).Number;
    end

    y1=normcdf(Xr,XrBar,XrSD);
    z1=y1;
    z1(z1>.5)=1-z1(z1>.5);
    [a b]=sort(z1,'ascend');
    for j=1:n
       Rcell{j+1,8}=y1(b(j));
       Rcell{j+1,9}=File1.NT(Indices1(b(j))).Number;
    end

    y2=normcdf(Yr,YrBar,YrSD);
    z2=y2;
    z2(z2>.5)=1-z2(z2>.5);
    [a b]=sort(z2,'ascend');
    for j=1:n
       Rcell{j+1,11}=y2(b(j));
       Rcell{j+1,12}=File1.NT(Indices1(b(j))).Number;
    end

    y3=normcdf(Zr,ZrBar,ZrSD);
    z3=y3;
    z3(z3>.5)=1-z3(z3>.5);
    [a b]=sort(z3,'ascend');

    for j=1:n
       Rcell{j+1,14}=y3(b(j));
       Rcell{j+1,15}=File1.NT(Indices1(b(j))).Number;
    end
  end

  if exist('rFindAlignmentDiscrepancies'),
    c=rFindAlignmentDiscrepancies(File1,Indices1,File2,Indices2,'nearest4');
    [a b]=sort(c,'descend');
    Rcell{1,17}='Discrep';
    for j=1:n
        Rcell{j+1,17}=c(b(j));
        Rcell{j+1,18}=[File1.NT(Indices1(b(j))).Base File1.NT(Indices1(b(j))).Number];
    end
  end

  if Comparison.RotationData > 1,
    xlswrite([File1.Filename ' ' File2.Filename ' Rotations_D2'],Rcell);
  end
end

% ---------------------------------------------- Plot translation vectors

if Comparison.TranslationFigure > 0
  fprintf('%s\n',ColorText);
  zFlushOutput;
  figure(Comparison.TranslationFigure);
  clf

  scatter3(Xt,Yt,Zt,DotSize,colo,'filled')
  if Comparison.TranslationColorbar > 0,
    colorbar('eastoutside')
  end
  colormap('default');
  map = colormap;
  colormap(map(8:56,:))

  if Octave == 0,
    set(gcf,'renderer','zbuffer');
  end

  title(['Translation differences between ' File1.Filename ' and ' File2.Filename ', numbers from ' File1.Filename]);
  if Octave == 0,
    rotate3d on
  end
  hold on
  ax=axis;
end

n=length(Xt);           %number of observation
Ts=[Xt;Yt;Zt];          %3xn matrix of all observations
XtBar=mean(Xt);         %scalar
YtBar=mean(Yt);         %scalar
ZtBar=mean(Zt);         %scalar
XtSD=std(Xt);           %scalar
YtSD=std(Yt);           %scalar
ZtSD=std(Zt);           %scalar
TBar=[XtBar;YtBar;ZtBar];            %3x1 vector
Tm=TBar * ones(1,length(Ts));        %3xn matrix where each column is the value of RBar
St=(Ts-Tm)*(Ts-Tm)'/length(Ts);   %3x3 covariance matrix
StINV=inv(St);                    %3x3 inverse covariance matrix
CorrT=inv(sqrt(diag(diag(St))))*St*inv(sqrt(diag(diag(St))));
TD=Ts-Tm;                        %3xn matrix-subtract RBar from each entry
DSQ=TD'*StINV*TD;                %nxn where diagonal elements are the standardized squared distance...
D=diag(diag(DSQ));               %nxn diagonal matrix with di-sq's on diagonal
D2t=D*ones(n,1);                  %nx1 vector of di-sq's
[a b]=sort(D2t,'descend');
Dist=sqrt(Xt.*Xt+Yt.*Yt+Zt.*Zt);
Numbers=' ';
I=1;
Tcell{1,2}='D2';
for j=1:n
   Tcell{j+1,2}=a(j);
   Tcell{j+1,3}=[File1.NT(Indices1(b(j))).Base File1.NT(Indices1(b(j))).Number];
   if Comparison.TranslationFigure > 0
      if D2t(b(j)) > Comparison.TranslationCutoff,
         if j>1 && abs(D2t(b(j))-D2t(b(j-1)))<.0001
            Numbers=[Numbers ', ' File1.NT(Indices1(b(j))).Number];
         else
            text(Xt(b(I))+(ax(2)-ax(1))/100,Yt(b(I))+(ax(4)-ax(3))/100,Zt(b(I)),Numbers,'Fontsize',7);
            Numbers=File1.NT(Indices1(b(j))).Number;
            I=j;
         end
      end
   end
end

if Comparison.TranslationData > 0,
    Dist=sqrt(Xt.*Xt+Yt.*Yt+Zt.*Zt);
    [a b]=sort(Dist,'descend');
    Tcell{1,5}='Dist';
    for j=1:n
      Tcell{j+1,5}=a(j);
      Tcell{j+1,6}=[File1.NT(Indices1(b(j))).Base File1.NT(Indices1(b(j))).Number];
    end

    if exist('mvncdf') && exist('normcdf'),
      y = mvncdf(Ts',TBar',St);
      z = y;
      z(z>.5)=1-z(z>.5);
      [a b]=sort(z,'ascend');
      for j=1:n
        Tcell{j+1,5}=y(b(j));
        Tcell{j+1,6}=File1.NT(Indices1(b(j))).Number;
      end

      y1=normcdf(Xt,XtBar,XtSD);
      z1=y1;
      z1(z1>.5)=1-z1(z1>.5);
      [a b]=sort(z1,'ascend');
      for j=1:n
        Tcell{j+1,8}=y1(b(j));
        Tcell{j+1,9}=File1.NT(Indices1(b(j))).Number;
      end

      y2=normcdf(Yt,YtBar,YtSD);
      z2=y2;
      z2(z2>.5)=1-z2(z2>.5);
      [a b]=sort(z2,'ascend');
      for j=1:n
        Tcell{j+1,11}=y2(b(j));
        Tcell{j+1,12}=File1.NT(Indices1(b(j))).Number;
      end

      y3=normcdf(Zt,ZtBar,ZtSD);
      z3=y3;
      z3(z3>.5)=1-z3(z3>.5);
      [a b]=sort(z3,'ascend');
      for j=1:n
        Tcell{j+1,14}=y3(b(j));
        Tcell{j+1,15}=File1.NT(Indices1(b(j))).Number;
      end
    end

    if exist('rFindAlignmentDiscrepancies'),
      c=rFindAlignmentDiscrepancies(File1,Indices1,File2,Indices2,'nearest4');
      [a b]=sort(c,'descend');
      Tcell{1,17}='Discrep';
      for j=1:n
        Tcell{j+1,17}=c(b(j));
        Tcell{j+1,18}=[File1.NT(Indices1(b(j))).Base File1.NT(Indices1(b(j))).Number];
      end
    end

    if Comparison.TranslationData > 1,
      xlswrite([File1.Filename ' ' File2.Filename ' Translations_D2_2'],Tcell);
    end
end

% ---------------------------------------------- Save data to pass back

Comparison.File1 = File1;
Comparison.Indices1 = Indices1;
Comparison.File2 = File2;
Comparison.Indices2 = Indices2;
Comparison.MutualDistances = DM;

Comparison.OverallRotation = RR;
Comparison.OverallTranslation = tt;
Comparison.Rotations = R;
Comparison.Translations = t;
Comparison.Rcell = Rcell;
Comparison.Tcell = Tcell;
Comparison.RotationCovariance = Sr;
Comparison.TranslationCovariance = St;
Comparison.RotationStandardDistances = D2r;
Comparison.TranslationStandardDistances = D2t;
