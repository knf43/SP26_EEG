function sp26_save_rejection_plot(EEG_before, EEG_after, datasetname, phase, event, plotpath)
% Generate a PNG showing which epochs got rejected by pop_eegthresh.
% Top: summary counts. Bottom: waveforms of rejected epochs (if any).
% Saves as <datasetname>_<phase>_<event>_rejections.png at 150 DPI.

% Identify rejected trial indices by diffing epoch counts
n_before = size(EEG_before.data, 3);
n_after  = size(EEG_after.data, 3);
n_rejected = n_before - n_after;

% Try to recover the rejection mask from EEG_before.reject (set by pop_eegthresh
% even with auto-reject on, the indices are preserved before deletion).
% Fall back to "first n_rejected epochs" if the field isn't there.
if isfield(EEG_before, 'reject') && isfield(EEG_before.reject, 'rejthresh') ...
        && ~isempty(EEG_before.reject.rejthresh)
    rejected_idx = find(EEG_before.reject.rejthresh);
else
    rejected_idx = []; % will be backfilled by re-running threshold check below
end

% If we still don't have the mask, just re-run thresh detection without rejecting
if isempty(rejected_idx) && n_rejected > 0
    [~, ~, ~, rej] = pop_eegthresh(EEG_before, 1, 1:size(EEG_before.data,1), ...
        -100, 100, EEG_before.xmin, EEG_before.xmax-0.002, 0, 0);
    rejected_idx = find(rej);
end

fig = figure('Visible', 'off', 'Position', [100 100 1200 800]);

% --- Top panel: summary text ---
subplot('Position', [0.05 0.85 0.9 0.12]);
axis off;
summary_lines = {
    sprintf('Participant: %s   |   Phase: %s   |   Event: %s', datasetname, phase, event), ...
    sprintf('Epochs before rejection: %d', n_before), ...
    sprintf('Epochs after rejection:  %d', n_after), ...
    sprintf('Epochs rejected:         %d  (%.1f%%)', n_rejected, 100 * n_rejected / max(n_before, 1))};
text(0, 0.5, summary_lines, 'FontSize', 12, 'FontName', 'FixedWidth', ...
     'VerticalAlignment', 'middle');

% --- Bottom panel: rejected epoch waveforms (or placeholder) ---
if n_rejected == 0
    % Empty placeholder
    subplot('Position', [0.05 0.05 0.9 0.75]);
    axis off;
    text(0.5, 0.5, 'No epochs rejected', 'FontSize', 24, ...
         'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
else
    % Plot up to first 40 rejected epochs in a grid
    n_to_plot = min(n_rejected, 40);
    n_cols = 8;
    n_rows = ceil(n_to_plot / n_cols);
    times = EEG_before.times;

    for i = 1:n_to_plot
        ep = rejected_idx(i);
        ax = subplot(n_rows + 1, n_cols, n_cols + i);
        plot(times, squeeze(EEG_before.data(:, :, ep))', 'LineWidth', 0.5);
        ylim([-150 150]);
        xlim([times(1) times(end)]);
        title(sprintf('Ep %d', ep), 'FontSize', 8);
        if mod(i-1, n_cols) ~= 0, set(ax, 'YTickLabel', []); end
        if i <= n_to_plot - n_cols, set(ax, 'XTickLabel', []); end
        set(ax, 'FontSize', 7);
        hline = refline(0, 100); hline.Color = [0.8 0.2 0.2]; hline.LineStyle = ':';
        hline = refline(0, -100); hline.Color = [0.8 0.2 0.2]; hline.LineStyle = ':';
    end

    if n_rejected > n_to_plot
        annotation('textbox', [0.05 0.01 0.9 0.03], ...
            'String', sprintf('(Showing first %d of %d rejected epochs)', n_to_plot, n_rejected), ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
    end
end

% --- Save ---
filename = strcat(datasetname,'_',phase,'_',event,'_rejections.png');
fullpath = fullfile(plotpath, filename);
print(fig, fullpath, '-dpng', '-r150');
close(fig);
fprintf('Saved rejection plot: %s\n', filename);
end
