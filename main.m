%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non Negative Hiddem Markov Model of Audio with Source Separation

% NB : We tried to keep the same notations as in the original paper
%       Probabilites are numbered according to the number in the paper

%¨TO DO
% 1 - Initialise transition probabilites (P4_mat, P5_mat, mu_vq, sigma_vq, A, pi)
% 2 - Loop :
% 2.1 - Compute likelihood (P3_mat)
% 2.2 - Compute alpha, beta
% 2.3 - Compute probabilites (P1_mat, P2_mat)   
% 2.4 - Update transition probabilities (P4_mat, P5_mat, mu_vq, sigma_vq, A, pi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% Set Parameters
T = 50; % Number of time frames
K = 10; % Number of vectors per dictionary ; num latent components ; num phonems
num_dicts = 40; % Number of dictionaries
num_freq = 100; % Number of frequencies in spectrogram

% Data
Vft = rand(T, num_freq); % Spectrogram (over all time frames)
v = sum(Vft, 2); % Observed energy

% 1 - Initialise transition probabilites
pi = 1/num_dicts * ones(num_dicts, 1); % Initial state
A = ones(num_dicts, num_dicts); % State transition
mu_vq = mean(v) * ones(num_dicts, 1); % Average power
sigma_vq = std(v) * ones(num_dicts, 1);
P4_mat = 1/num_freq * ones(num_freq, K, num_dicts) ; % P(f|z,q)
P5_mat = 1/K * ones(T, K, num_dicts); % P(zt|qt)

tic()

% 2.1 - Compute likelihood (actually : log likelihood)
P3_mat = zeros(T, num_dicts); % P(ft_bold, vt | qt)
for t = 1:T
   for qt = 1:num_dicts
       P3_mat(t, qt) = P_vq(v(t), qt, mu_vq, sigma_vq); % first term
       for ft = 1:num_freq
           P3_mat(t, qt) = P3_mat(t, qt) * (P4_mat(ft, :, qt) * P5_mat(t, :, qt)') ^ Vft(t, ft);
       end
   end
end

% 2.2 - Compute alpha, beta
log_alpha = zeros(T, num_dicts);
log_beta = zeros(T, num_dicts);

log_alpha(1, :) = log(pi .* P3_mat(1,:)');
log_beta(T, :) = log(1);

for t = 2:T
    for qt = 1:num_dicts
        temp = log(A(:, qt)') + log_alpha(t-1, :);
        [max_val, max_ind] = max(temp);
        temp_but_max = [temp(1:max_ind-1), temp(max_ind+1:end)];      
        log_alpha(t, qt) = log(P3_mat(t, qt)) ...
                    + max_val + log(1 + sum(exp(temp_but_max - max_val)));
    end
end

for t = 1:T-1
    for qt=1:num_dicts
        temp = (log(A(qt, :)') + log(P3_mat(T-t, :)') + log_beta(T-t+1, :)')';
        [max_val, max_ind] = max(temp);
        temp_but_max = [temp(1:max_ind-1), temp(max_ind+1:end)];
        log_beta(T-t, qt) = max_val + log(1 + sum(exp(temp_but_max - max_val)));        
    end    
end

% 2.3 - Compute probabilites (P1_mat, P2_mat)

% Probability of having hidden state qt at time t
p_qt = zeros(T,num_dicts);
for t=1:T
    temp = log_alpha(t,:) + log_beta(t,:);
    for qt = 1:num_dicts
        p_qt(t, qt) = 1/sum(exp(temp - temp(qt)));
    end
end

% Probability of state transition
p_qt_qt_plus_one = zeros(T-1,num_dicts,num_dicts);
for t=1:T-1
    for qt = 1:num_dicts
        for qt_plus_one = 1:num_dicts
            p_qt_qt_plus_one(t, qt, qt_plus_one) = 1/sum(exp(log_alpha(t,:) + log_beta(t,:) - log_alpha(t,qt) - log_beta(t+1, qt_plus_one)))...
                     * A(qt, qt_plus_one) * P3_mat(t+1, qt_plus_one);
        end
    end
end


% P2_mat : P(zt, ft | qt)
P2_mat = zeros(T, K, num_freq, num_dicts);
for t = 1:T
    for ft = 1:num_freq
       for qt = 1:num_dicts
           P2_mat(t, :, ft, qt) = P5_mat(t, :, qt) .* P4_mat(ft, :, qt) / sum(P5_mat(t, :, qt) .* P4_mat(ft, :, qt));
       end
    end
end

% P1_mat : P(zt, qt | ft, f_bold, v_bold)
P1_mat = zeros(T, K, num_dicts, num_freq);
for t = 1:T
    for ft = 1:num_freq
       for qt = 1:num_dicts
           P1_mat(t, :, ft, qt) = p_qt(t, qt) * P2_mat(t, :, ft, qt);
       end
    end
end


% Update P4_mat
P4_mat = zeros(num_freq, K, num_dicts) ; % P(f|z,q)
for zt = 1:K
    for qt = 1:num_dicts
        for f = 1:num_freq
            P4(f, zt, qt) = Vft(:, f)' * P1_mat(:, zt, qt, f);
        end
        P4(:, zt, qt) = P4(:, zt, qt) / sum(P4(:, zt, qt));
    end
end


% Update P5_mat
P5_mat = zeros(T, K, num_dicts); % p(z|q)
for t=1:T
    for qt = 1:num_dicts
        for zt = 1:K
            P5_mat(t, zt, qt) = Vft(t,:) * squeeze(P1_mat(t, zt, qt, :));
        end
        P5_mat(t, :, qt) = P5_mat(t, :, qt) / sum(P5_mat(t, :, qt));
    end
end

% Update pi (start probability)
for q=1:K
    pi(q) = p_qt(1,q) / sum(p_qt(1,:));
end

% Update A (transition matrix)
for q1=1:K
    for q2=1:K
        A(q1,q2) = sum(p_qt_qt_plus_one(:, q1,q2)) / sum(sum(p_qt_qt_plus_one(:,q1,:)));
    end
end


toc()