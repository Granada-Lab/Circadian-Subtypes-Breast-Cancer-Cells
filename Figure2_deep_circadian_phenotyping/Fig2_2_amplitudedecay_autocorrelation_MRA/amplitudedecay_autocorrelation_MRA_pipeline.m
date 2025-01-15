function amplitudedecay_autocorrelation_MRA_pipeline

%Carolin Ector, 13.11.2023

% Function loads circadian luciferase signals from different cell line models and circadain reporters (Bmal1 and Per2) and runs the following functions sequentially:
% function "f1_sort_lumicycle_and_mra_data.m": sorts datafiles to read signals for all replicates per cell line.
% function "f2_overlay_circadian_signals_mean.m": overlays unnormalized or normalized detrended signals from the Bmal1 and Per2 luciferase reporters per cell line, averaged for all replicates.
% function "f3_fit_exponential_decay_amplitude.m": fits an exponential function to a circadian signal's envelope and calculates the exponential decay rate.
% function "f4_autocorrelation_calculation.m": calculates the autocorrelation of a signal.

%Clock-TNBC Manuscript Fig. 2

%input: stored in "amplitudedecay_autocorrelation_MRA_pipeline.mat"
% lumicycle_datafile: excel file where processed circadian time-series data is stored and multi-resolution values (output from python function)
% celllinenames_file: names of the cell lines being analysed, as written in the file names
% reporters: circadian gene names for the luciferase reporters
% reporter_colors: colors used in the graphs for the two circadian clock luciferase reporters
% reporter_colors_dark: darker version of colorbp to use for text and thin lines
% axOpt: setting for appearance of axes of a plot
% recordingtime_all: total recording time as a time series (1-137.7 hours)

load(['amplitudedecay_autocorrelation_MRA_pipeline.mat'])

%% sort data --> extract all replicates of the same reporter cell line and save to cell array 
[detrended_all, envelope_all, norm_all] = f1_sort_lumicycle_and_mra_data(lumicycle_datafile,celllinenames_file,reporters);

%% overlay detrended signals from the Bmal1- and Per2-Luc reporters (mean all replicates)
f2_overlay_circadian_signals_mean(detrended_all,norm_all,celllinenames_file,reporters,reporter_colors,reporter_colors_dark,axOpt,recordingtime_all)

%% fit an exponential decaying function to the amplitude envelope to extract decay rates (sample-by-sample)
f3_fit_exponential_decay_amplitude(detrended_all,envelope_all,celllinenames_file,reporters,reporter_colors,axOpt,recordingtime_all)

%% calculate autocorrelation function
f4_autocorrelation_calculation(detrended_all,reporters,celllinenames_file,reporter_colors,axOpt)

%% create scatterplots for autocorrelation and MRA values
f5_scatterplot_autocorrelation_mra(reporters,celllinenames_file,subtype_colors,markers)

end %function
