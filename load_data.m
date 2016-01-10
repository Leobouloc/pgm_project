function [int_train_spectrogram1, train_audio1, int_train_spectrogram2, train_audio2, ...
          int_test_spectrogram1, test_audio1, int_test_spectrogram2, test_audio2, ...
          int_mixture_spectrogram, audio_mixture, ...
          int_size, complex_mixture_spectrogram, sr, g, hop_size, window_size] = load_data(display_bool)

close all

addpath('ltfat/');
ltfatstart; % used for fast STFT

%%%%%%%%%%%
% Options %
%%%%%%%%%%%

%int_size = 8;
int_size = 16;

training_samples_proportion = 0.01; % proportion of training samples to use (0 means no shuffle)
test_samples_proportion = 0.5;
randomize_mixture_order = true;

training_folder = '../PCCdata16kHz_train/train/reverberated/';
test_folder = '../PCCdata16kHz_test/devel/isolated/clean/';

if nargin < 1
    display = 1;
else
    display = display_bool;
end

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%

speaker1 = 8; % male
speaker2 = 25; % Emma <3

%%%%%%%%%%%%%%%%%
% Training data %
%%%%%%%%%%%%%%%%%

[train_spectrogram1, train_audio1, ~, ~, sr1, T1] = get_training_spectrogram(speaker1, training_folder, training_samples_proportion);
[train_spectrogram2, train_audio2, ~, ~, sr2, T2] = get_training_spectrogram(speaker2, training_folder, training_samples_proportion);

if sr1 == sr2
    sr = sr1;
else
    error('Training audio files were not formatted the same way for both speakers!');
end

if display
    title1 = sprintf('Training spectrogram for speaker %d', speaker1);
    title2 = sprintf('Training spectrogram for speaker %d', speaker2);
    display_spectrogram(train_spectrogram1, T1, sr, title1);
    display_spectrogram(train_spectrogram2, T2, sr, title2);
end

switch(int_size)
    case 8
        int_train_spectrogram1 = im2uint8(train_spectrogram1);
        int_train_spectrogram2 = im2uint8(train_spectrogram2);
    case 16
        int_train_spectrogram1 = im2uint16(train_spectrogram1);
        int_train_spectrogram2 = im2uint16(train_spectrogram2);
    otherwise
        error('Only depths of 8 bits and 16 bits were implemented.');
end

%%%%%%%%%%%%%
% Test data %
%%%%%%%%%%%%%

[mixture_spectrogram, complex_mixture_spectrogram, audio_mixture, test_spectrogram1, test_audio1, test_spectrogram2, test_audio2, sr, T, g, hop_size, window_size] = merge_test_sounds(speaker1, speaker2, test_folder, randomize_mixture_order, test_samples_proportion);

if display
    % Listen to the sounds (plays all at once, 'clear all' to stop)
    soundsc(audio_mixture, sr);
    %soundsc(test_audio1, sr);
    %soundsc(test_audio2, sr);
    display_spectrogram(mixture_spectrogram, T, sr, 'Mixture spectrogram for the test data');
end

switch(int_size)
    case 8
        int_test_spectrogram1 = im2uint8(test_spectrogram1);
        int_test_spectrogram2 = im2uint8(test_spectrogram2);
        %int_mixture_spectrogram = im2uint8(mixture_spectrogram);
    case 16
        int_test_spectrogram1 = im2uint16(test_spectrogram1);
        int_test_spectrogram2 = im2uint16(test_spectrogram2);
        %int_mixture_spectrogram = im2uint16(mixture_spectrogram);
    otherwise
        error('Only depths of 8 bits and 16 bits were implemented.');
end

int_mixture_spectrogram = int_test_spectrogram1 + int_test_spectrogram2;

end
