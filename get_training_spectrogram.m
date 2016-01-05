function [spectrogram, window_size, hop_size, sr, T] = get_training_spectrogram(speaker_id, training_folder, training_samples_proportion)

[audio, sr] = get_training_audio(speaker_id, training_folder, training_samples_proportion);
T = length(audio);

% soundsc(audio, sr); % to play what is sounds like (to stop it: clear all)
% NB: while the first one will keep playing, the second will start as soon
% as it's computed.

% Fast STFT
window_size = 1024;
hop_size = window_size/4;
g = gabwin({'tight', 'hann'}, hop_size, window_size, window_size);
[G, ~] = dgtreal(audio, g, hop_size, window_size);

spectrogram = abs(G);
%spectrogram = 20*log(abs(G));

end

function [y, sr] = get_training_audio(speaker_id, training_folder, proportion)
%GET_TRAINING_WAV

dirname = [training_folder, 'id', int2str(speaker_id), '/'];

files = dir([dirname, '/*.wav']);
if proportion > 0
    nb_files = length(files);
    selection = randsample(nb_files, ceil(proportion*nb_files), false);
    selected_files = files(selection);
else
    selected_files = files;
end

y = [];
sr_old = 0;

for file = selected_files'
    filename = [dirname, file.name];
    [y_file, sr] = audioread(filename);
    
    % Check sample rate
    if sr_old == 0 % first iteration
        sr_old = sr;
    elseif sr ~= sr_old
        error('All audio files should have the same sampling frequency!');
    end
    
    % Conversion to mono when required
    if size(y_file, 2) > 1
        y_file = mean(y_file, 2);
    end
    
    % Add bit to the training audio file
    y = [y; y_file];
end

end