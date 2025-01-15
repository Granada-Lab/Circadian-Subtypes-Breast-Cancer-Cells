function cwt_analysis_pipeline

%Carolin Ector, 02.11.2023

%Function runs analysis on continuous wavelet analysis readout files extracted with "pyBOAT_MRA_pipeline.py" (run beforehand)

%Runs the following functions sequentially
%cwt_extract_ridge.m to extract a continuous ridge from the wavelet spectrum and determine its length
%cwt_plot_ridge_readout.m to plot the readout along the ridge as a function of time
%cwt_phase_difference_weighted.m to calculate and plot the phase difference between Bmal1 and Per2 reporters using circular statistics
%cwt_boxplot_rankings.m to rank cell lines by their circadian
%parameters obtained from continuous wavelet transform

%Clock-TNBC Manuscript Fig. 2; the extracted parameters are used throughout the manuscript

%input: stored in "cwt_analysis_pipeline.mat"
% pathtocwt_readout: path to folder that contains the continuous wavelet transform ridge readout files
% celllinenames_file: names of the cell lines being analysed, as written in the file names
% reporters: Luciferase reporters being analyzed 
% subtype_colors: color coding for each subtype, cellline-by-cellline
% reporter_colors: colors for each Luc-reporter

load('cwt_analysis_pipeline.mat')

%define different names for Bmal1 and Per2 in the files
bmal = {'BLH','BMAL','Bmal','bmal'};
per = {'PER','Per','per'};

%% load data

%define folder names that contain the ridge readout files.

% global = 1/4 of median half-maximal wavelet power across all samples(fixed threshold = 118.4/4) --> amplitude, ridge length
folderPath_global = append(pathtocwt_readout,'/','cwt_ridge_readout_unnormalized_detrended'); 
% adaptive = 1/4 of maximal wavelet power per sample (adaptive threshold)--> period, phase difference
folderPath_adaptive = append(pathtocwt_readout,'/','cwt_ridge_readout_normalized_detrended'); 

%list file names.
fileList_global = dir(folderPath_global);
fileList_adaptive = dir(folderPath_adaptive);

%convert all file names to lowercase for case-insensitive comparison
AllFiles_global = ({fileList_global.name});
AllFiles_adaptive = ({fileList_adaptive.name});
is_bmal = false(1, length(AllFiles_adaptive));
is_per = false(1, length(AllFiles_adaptive));

% check for reporter names
for b = 1:numel(bmal)
    is_bmal = is_bmal | contains(AllFiles_global, bmal{b});
end
for p = 1:numel(per)
    is_per = is_per | contains(AllFiles_global, per{p});
end

%% process ridge lengths (identify discontiunity), extract ridge length parameter.
cwt1_extract_ridge(folderPath_global,folderPath_adaptive,AllFiles_global,AllFiles_adaptive,is_bmal,is_per,celllinenames_file,reporters)

%% extract amplitude and period parameters; sort continuous amplitudes and periods per cell model (all replicates).
cwt2_plot_ridge_readout(folderPath_global,folderPath_adaptive,AllFiles_global,AllFiles_adaptive,is_bmal,is_per,celllinenames_file,reporters)

%% calculate continuous Bmal1-Per2 phase differences; extract phase difference parameters.
cwt3_extract_phase_difference(folderPath_adaptive,bmal,per,celllinenames_file);

%% calculate weights for continuous wavelet transform parameters.
cwt4_weighting_of_cwt_parameters(celllinenames_file,reporter_colors)

%% plot histogram to see distribution of ridge lengths for respective threshold
cwt5_plot_ridge_readout_distribution(reporters,reporter_colors,celllinenames_file,reporter_colors_dark)

%% calculate weighted phase difference between Bmal1- and Per2-reporters (by cell line and by subtype)
cwt6_plot_weighted_phase_difference(celllinenames_file)

end %function