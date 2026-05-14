function y = sp26_post_ica(datasetnumb, lang, phase, event)
dataset     = strcat(datasetnumb,'_',lang,'.set');
datasetname = strcat(datasetnumb,'_',lang);
savename    = strcat(datasetname,'_',phase,'_',event,'.set');
loadpath    = '/Users/kayleefernandez/Desktop/SP26_data/1_eeg_data/preprocessed/';
savepath    = '/Users/kayleefernandez/Desktop/SP26_data/1_eeg_data/preprocessed/post_ica/';
plotpath    = '/Users/kayleefernandez/Desktop/SP26_data/1_eeg_data/preprocessed/rejection_plots/';
EEG = pop_loadset('filename', dataset, 'filepath', loadpath);
event_offsets = struct('sentence',0,'noun_onset',1,'noun_F0',2,'verb_onset',3,'verb_F0',4,'suffix',5);
event_windows = struct('sentence',[-0.2 2.5],'noun_onset',[-0.2 1.0],'noun_F0',[-0.2 1.0],'verb_onset',[-0.2 1.0],'verb_F0',[-0.2 1.0],'suffix',[-0.2 1.0]);
offset = event_offsets.(event);
window = event_windows.(event);
if strcmp(phase, 'testing'), cond_bases = [10 20 30 40 50 60]; else, cond_bases = [70 80]; end
trigger_codes = cond_bases + offset;
trigger_strings = arrayfun(@(x) sprintf('S %d', x), trigger_codes, 'UniformOutput', false);
EEG = pop_epoch(EEG, trigger_strings, window, 'newname', [datasetname '_' phase '_' event], 'epochinfo', 'yes');
EEG_before_rej = EEG;
EEG = pop_eegthresh(EEG, 1, 1:30, -100, 100, window(1), window(2)-0.002, 1, 0);
sp26_save_rejection_plot(EEG_before_rej, EEG, datasetname, phase, event, plotpath);
eeglab redraw;
EEG.setname = strcat(datasetname,'_',phase,'_',event);
EEG = pop_saveset(EEG, 'filename', savename, 'filepath', savepath);
y = EEG;
