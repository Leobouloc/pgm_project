function [audio1, audio2, sr] = mixture_weights_to_audio(P_final_s1, complex_mixture, sr, g, hop_size, window_size)
%MIXTURE_WEIGHTS_TO_AUDIO

mixture_weights1 = P_final_s1';

complex_spectrogram1 = mixture_weights1() .* complex_mixture;
complex_spectrogram2 = (1 - mixture_weights1) .* complex_mixture;

audio1 = idgtreal(complex_spectrogram1, g, hop_size, window_size);
audio2 = idgtreal(complex_spectrogram2, g, hop_size, window_size);

soundsc(audio1, sr);
soundsc(audio2, sr);
end

