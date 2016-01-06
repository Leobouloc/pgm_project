function display_dictionaries(P4_mat, sr)
%DISPLAY_DICTIONARIES displays dictionaries from probabilities stores in
%matrix P4_mat = p(f|z, q) with size num_freq x K x num_dicts

num_freq = size(P4_mat, 1);
K = size(P4_mat, 2);
num_dicts = size(P4_mat, 3);

f_range = linspace(0, sr/2, num_freq);

figure('name', 'Dictionaries')
subplot1(1, num_dicts);
for q = 1:num_dicts
    subplot1(q);
    imagesc(1:K, f_range, 20*log(abs(P4_mat(:, :, q)))); % displayed in magnitude
    axis([1 K 0 sr/2]);
    set(gca,'YDir','normal');
    drawnow;
end

end

