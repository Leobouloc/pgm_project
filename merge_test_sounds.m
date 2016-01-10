function [mixture_spectrogram, audio_mixture, ...
          spectrogram1, audio1, spectrogram2, audio2, sr, T] = ...
          merge_test_sounds(speaker1_id, speaker2_id, test_folder, randomize_bool, test_samples_proportion)
%MERGE_TEST_SOUNDS

[audio1, sr1] = get_test_audio(speaker1_id, test_folder, randomize_bool, test_samples_proportion);
[audio2, sr2] = get_test_audio(speaker2_id, test_folder, randomize_bool, test_samples_proportion);

if sr1 == sr2
    sr = sr1;
else
    error('Test audio files were not formatted the same way for both speakers!');
end

T1 = length(audio1);
T2 = length(audio2);

if T1 > T2
    audio2 = [audio2; zeros(T1 - T2, 1)];
    T = T1;
else
    audio1 = [audio1; zeros(T2 - T1, 1)];
    T = T2;
end
audio_mixture = (audio1 + audio2)/2;

% Fast STFTs
window_size = 1024;
hop_size = window_size/4;
g = gabwin({'tight', 'hann'}, hop_size, window_size, window_size);
[G1, ~] = dgtreal(audio1, g, hop_size, window_size);
[G2, ~] = dgtreal(audio2, g, hop_size, window_size);
G = dgtreal(audio_mixture, g, hop_size, window_size);

spectrogram1 = abs(G1);
%spectrogram1 = 20*log(abs(G1));
spectrogram2 = abs(G2);
%spectrogram2 = 20*log(abs(G2));
mixture_spectrogram = abs(G);
%mixture_spectrogram = 20*log(abs(G));

end

function [y, sr] = get_test_audio(speaker_id, test_folder, randomize_bool, proportion)
%GET_TRAINING_WAV

prefix = [test_folder, 's', int2str(speaker_id)];

files = dir([prefix, '_*.wav']);
nb_files = length(files);
if randomize_bool
    permutation = randsample(nb_files, ceil(proportion*nb_files), false);
    perm_files = files(permutation);
else
    perm_files = files(1:ceil(proportion*nb_files));
end

y = [];
sr_old = 0;

for file = perm_files'
    filename = [test_folder, file.name];
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


