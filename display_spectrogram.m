function display_spectrogram(spectrogram, T, sr, title)
%PLOT_SPECTROGRAM 

if nargin < 3
    title = 'Spectrogram';
end

figure('name', title);
t_range = linspace(0, T/sr, size(spectrogram, 2));
f_range = linspace(0, sr/2, size(spectrogram, 1));
imagesc(t_range, f_range, 20*log(abs(spectrogram)))
set(gca,'YDir','normal');
drawnow;

end