% The code extracts results from CNMF algorithm.
% Author: Weijian Yang and Darcy S. Peterka

% toolbox:
% plotROIContour: plot the ROI contour ()
% signalExtraction: extract the raw, filtered, and inferred signal (dfof, F, df)
% plotActivityTrace: plot the signal

addpath(genpath('utilities'));

%% plot ROI contour
% [ CC, info ] = plotROIContour( A, d1, d2, plotControl )
% input parameters
% A: spatial component from CNMF 
% d1, d2: dimension of the data
% plotControl.Cn: background image
% plotControl.displayLabel: if it is 1, the ROI number will be labelled on the plot
% plotControl.ROIID: id of ROIs to be plotted
% plotControl.option: 1: The contour is drawn around the value above which a specified fraction of energy is explained (default 99%)
%                     2: The contour is drawn around the exterior boundary 
% plotControl.thr: energy fraction included in the drawn contour in option 1
% plotControl.cmap: colormap of the background image
% plotControl.lineColor: line color of the drawn contour

% output parameters
% CC: outline of the contour
% info: contour information, contains ROI id, weighted contour map, boundary box, and centroid

% Author: Eftychios A. Pnevmatikakis, Weijian Yang and Darcy S. Peterka 

plotControl.Cn=Cn;
plotControl.thr=0.95;
plotControl.option=1;
d1=options.d1;
d2=options.d2;
[coor, contourInfo] = plotROIContour( A2, d1, d2, plotControl );

%% extract the signal
% [ inferred, filtered, raw ] = signalExtraction(Y,A,C,b,f,d1,d2,extractControl)
% inputs: Y raw data (d X T matrix, d # number of pixels, T # of timesteps)
%         A matrix of spatial components (d x K matrix, K # of components)
%         C matrix of temporal components (K x T matrix)
%         S matrix of deconvolved activity ((K-1) x T matrix) (optional)
%         b matrix of background spatial components (d x options.nb, options.nb # of components)
%         f matrix of background temporal components (options.nb x T matrix, options.nb # of components) 
%         extractControl.baselineRatio: residue baseline of signals from CNMF (default:0.25); set to 0 to ignore this residue baseline 
%         extractControl.thr: energy fraction included in the raw signal contour (default:0.99)

% inferred: spatial weighting on the ROI pixels, unmixing, background substraction, and denoising
% filtered: spatial weighting on the ROI pixels, unmixing, background substraction
% raw: uniform spatial weighting on the ROI pixels (with threshold to remove the very low energy pixels), background substraction
% inferred, filtered, raw contains structure fields: df, F, dfof, with definition below
% df: df signal
% F: baseline
% dfof: dfof
% in addition, raw contains a structure field named "CR", with definition below
% CR: information of the cell contour used in raw signal calculation

% Author: Eftychios A. Pnevmatikakis, and Weijian Yang

extractControl.baselineRatio = 0.25;
extractControl.thr = 0.95;
[inferred, filtered, raw] = signalExtraction(Yr,A2,C2,b2,f2,d1,d2,extractControl);

%% plot the signals
% plotActivityTrace( inferred, filtered, raw, plotControl)
% inputs: inferred, filtered, raw (from the output of "signalExtraction.m")
%         plotControl.frameRate: frame rate in fps
%         plotControl.normalization: normalize the amplitude of the temporal trace for plotting purpose
%         plotControl.sep: vertical (amplitude) offset between adjacent temporal traces 
%         plotControl.displayLabel: whether to display ROI ID for the temporal traces 
%         plotControl.plotInferred: whether the inferred signal is plotted
%         plotControl.plotFiltered: whether the filtered signal is plotted
%         plotControl.plotRaw: whether the raw signal is plotted
%         plotControl.rollingView: whether the display is set such that only 20 temporal traces are displayed at a time, followed by any key to move onto the next 20 traces, and etc.

% Author: Weijian Yang

plotControl2.frameRate=6.7468;
plotControl2.normalization=1;   plotControl2.sep=1.2;   plotControl2.displayLabel=1;  plotControl2.rollingView=1;
plotControl2.plotInferred=1;    plotControl2.plotFiltered=1;   plotControl2.plotRaw=0;
plotActivityTrace( inferred.dfof, filtered.dfof, raw.dfof, plotControl2 );
