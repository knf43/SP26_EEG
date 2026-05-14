% sp26_post_ica_batch.m
% Runs sp26_post_ica for every participant x phase x event.
% Stops on first error. Edit the CONFIG block below to change behavior.

%% ===== CONFIG =====
participants = {
    {'0101','L1'}, {'0103','L1'}, {'0104','L1'}, {'0105','L1'}, ...
    {'0106','L1'}, {'0107','L1'}, {'0108','L1'}, {'0109','L1'}, ...
    {'0110','L1'}, {'0111','L1'}, ...
    {'0201','L2'}, {'0202','L2'}, {'0203','L2'}, {'0204','L2'}, ...
    {'0205','L2'}, {'0206','L2'}, {'0207','L2'}, {'0208','L2'}, ...
    {'0209','L2'} ...
};

phases = {'testing','training'};
events = {'sentence','noun_onset','noun_F0','verb_onset','verb_F0','suffix'};

FORCE_REPROCESS = false;   % true = redo files that already exist
%% ===================

savepath = '/Users/kayleefernandez/Desktop/SP26_data/1_eeg_data/preprocessed/post_ica/';

total = length(participants) * length(phases) * length(events);
count = 0;
skipped = 0;

fprintf('\n=== Starting batch: %d total files to process ===\n\n', total);
tic;

for p = 1:length(participants)
    pid  = participants{p}{1};
    lang = participants{p}{2};

    for ph = 1:length(phases)
        phase = phases{ph};

        for e = 1:length(events)
            event = events{e};
            count = count + 1;

            savename = strcat(pid,'_',lang,'_',phase,'_',event,'.set');
            fullpath = fullfile(savepath, savename);

            fprintf('[%d/%d] %s ... ', count, total, savename);

            if ~FORCE_REPROCESS && exist(fullpath, 'file')
                fprintf('SKIP (exists)\n');
                skipped = skipped + 1;
                continue;
            end

            sp26_post_ica(pid, lang, phase, event);
            fprintf('done\n');
        end
    end
end

elapsed = toc;
fprintf('\n=== Batch complete: %d processed, %d skipped, %.1f min ===\n', ...
        count - skipped, skipped, elapsed/60);
