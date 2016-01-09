Vft = int_mixture_spectrogram;

% Set Parameters
T = size(Vft, 1); % Number of time frames (remove +1 if no white noise)
K = 2; % Number of vectors per dictionary ; num latent components ; num phonems
num_dicts1 = 20; % Number of dictionaries for source 1
num_dicts2 = 20; % Number of dictionaries for source 2
num_freq = size(Vft, 2); % Number of frequencies in spectrogram

% Adding white noise at end to avoid zero probabilities
Vft(T + 1, :) = 1; % WARNING : change T if you remove this
T = T+1;

% Compute Energy
v = sum(Vft, 2); % Observed energy

% 1 - Initialise transition probabilites
Pi1 = 1/num_dicts1 * ones(num_dicts1, 1); % Initial state for source 1
Pi2 = 1/num_dicts2 * ones(num_dicts2, 1); % Initial state for source 2
A1 = 1/num_dicts1 * ones(num_dicts1, num_dicts1); % State transition for source 1
A2 = 1/num_dicts2 * ones(num_dicts2, num_dicts2); % State transition for source 2
mu1_vq = mean(v)/2 * ones(num_dicts1, 1); % Average power of source 1
sigma1_vq = std(v) * ones(num_dicts1, 1);
mu2_vq = mean(v)/2 * ones(num_dicts2, 1); % Average power of source 2
sigma2_vq = std(v) * ones(num_dicts2, 1);

% 2.2 - Compute alpha, beta
log_alpha = zeros(T, num_dicts1, num_dicts2);
log_beta = zeros(T, num_dicts1, num_dicts2);

log_alpha(1, :) = log(pi) .* P3_mat(1,:)';
log_beta(T, :) = log(1);

for t = 2:T
    for qt = 1:num_dicts
        temp = log(A(:, qt)') + log_alpha(t-1, :);
        [max_val, max_ind] = max(temp);
        temp_but_max = [temp(1:max_ind-1), temp(max_ind+1:end)];      
        log_alpha(t, qt) = P3_mat(t, qt) ...
                    + max_val + log(1 + sum(exp(temp_but_max - max_val)));
    end
end

for t = 1:T-1
    for qt=1:num_dicts
        temp = (log(A(qt, :)') + P3_mat(T-t+1, :)' + log_beta(T-t+1, :)')';
        [max_val, max_ind] = max(temp);
        temp_but_max = [temp(1:max_ind-1), temp(max_ind+1:end)];
        log_beta(T-t, qt) = max_val + log(1 + sum(exp(temp_but_max - max_val)));        
    end    
end