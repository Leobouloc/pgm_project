%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non Negative Hidden Markov Model of Audio with Source Separation

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

% Data
[Vft, Vft2, int_test_spectrogram1, int_test_spectrogram2, int_mixture_spectrogram] = load_data(0);
Vft = double(Vft');
%Vft2 = double(Vft2');

%load('Vft.mat')

% Set Parameters
T = size(Vft, 1); % Number of time frames (remove +1 if no white noise)
K = 4; % Number of vectors per dictionary ; num latent components ; num phonems
num_dicts = 16; % Number of dictionaries
num_freq = size(Vft, 2); % Number of frequencies in spectrogram

% Adding white noise at end to avoid zero probabilities
Vft(T + 1, :) = 1; % WARNING : change T if you remove this
T = T+1;

% Filling blancs by white noise for 0/ NaN issues with P5_mat
Vft(all(Vft == 0, 2), :) = 1;

% Compute Energy
v = sum(Vft, 2); % Observed energy

% For smoothing
my_epsilon = 10^-10;

% 1 - Initialise transition probabilites
pi = log(1/num_dicts * ones(num_dicts, 1)); % Initial state
A = eye(num_dicts, num_dicts) + 1/10 * rand(num_dicts, num_dicts); % State transition
mu_vq = mean(v) * ones(num_dicts, 1); % Average power
sigma_vq = std(v) * ones(num_dicts, 1);
P4_mat = 1/num_freq * rand(num_freq, K, num_dicts); % P(f|z,q)
P5_mat = 1/K * rand(T, K, num_dicts); % P(zt|qt)

A = log(A ./ repmat(sum(A, 2), 1, num_dicts));
P4_mat = P4_mat ./ repmat(sum(P4_mat), num_freq, 1);
P5_mat = P5_mat ./ repmat(sum(P5_mat, 2), 1, K, 1);

log_likelihoods = [];

num_iterations = 20;
for xxx = 1:num_iterations
    tic()
    fprintf('***********************************\n - At iteration %d\n', xxx)
    % 2.1 - Compute likelihood (log-likelihood actually) (Vect)
    tic()
    P3_mat = zeros(T, num_dicts); % P(ft_bold, vt | qt)
    % first term
    for qt = 1:num_dicts
           P3_mat(:, qt) = log(P_vq(v, qt, mu_vq, sigma_vq)); 
       for ft = 1:num_freq
           P3_mat(:, qt) = P3_mat(:, qt)' +  Vft(:, ft)' .* log((P4_mat(ft, :, qt) * P5_mat(:, :, qt)'));
       end
    end
    if any(any(any(isnan(P3_mat))))    
        'h1'
        break
    end
    %tic()

    % 2.2 - Compute alpha, beta
    log_alpha = zeros(T, num_dicts);
    log_beta = zeros(T, num_dicts);

    log_alpha(1, :) = pi + P3_mat(1,:)';
    log_beta(T, :) = log(1);

    for t = 2:T
        for qt = 1:num_dicts
            temp = A(:, qt)' + log_alpha(t-1, :);
            [max_val, max_ind] = max(temp);
            temp_but_max = [temp(1:max_ind-1), temp(max_ind+1:end)];      
            log_alpha(t, qt) = P3_mat(t, qt) ...
                        + max_val + log(1 + sum(exp(temp_but_max - max_val)));
        end
    end

    for t = 1:T-1
        for qt=1:num_dicts
            temp = (A(qt, :)' + P3_mat(T-t+1, :)' + log_beta(T-t+1, :)')'; % +1 ???
            [max_val, max_ind] = max(temp);
            temp_but_max = [temp(1:max_ind-1), temp(max_ind+1:end)];
            log_beta(T-t, qt) = max_val + log(1 + sum(exp(temp_but_max - max_val)));        
        end    
    end
    
    if any(any(any(isnan(log_alpha))))
        '1'
        break
    end
    
    if any(any(any(isnan(log_beta))))
        '2'
        break
    end

    % 2.3 - Compute probabilites (P1_mat, P2_mat)
    %tic()
    % Probability of having hidden state qt at time t
    p_qt = zeros(T,num_dicts);
    %for t=1:T
    %    temp = log_alpha(t,:) + log_beta(t,:);
    %    for qt = 1:num_dicts
    %        p_qt(t, qt) = 1/sum(exp(temp - temp(qt)));
    %    end
    %end
    
    for t=1:T
        temp = log_alpha(t,:) + log_beta(t,:);
        [max_val, max_ind] = max(temp);
        temp_but_max = [temp(1:max_ind-1), temp(max_ind+1:end)];      
        for qt = 1:num_dicts
            p_qt(t,qt) = log_alpha(t, qt) + log_beta(t, qt) - max_val - log(1 + sum(exp(temp_but_max - max_val)));
        end
    end
    
    if any(any(any(isnan(p_qt))))
        '3'
        break
    end

    % Probability of state transition
    p_qt_qt_plus_one = zeros(T-1,num_dicts,num_dicts);
    for t=1:T-1
        temp = log_alpha(t,:) + log_beta(t,:);
        [max_val, max_ind] = max(temp);
        temp_but_max = [temp(1:max_ind-1), temp(max_ind+1:end)];     
        for qt = 1:num_dicts
            for qt_plus_one = 1:num_dicts
                p_qt_qt_plus_one(t, qt, qt_plus_one) = log_alpha(t, qt) ...
                    + log_beta(t+1, qt_plus_one) ...
                    + A(qt, qt_plus_one) ...
                    + P3_mat(t+1, qt_plus_one) ...
                    - max_val - log(1 + sum(exp(temp_but_max - max_val)));
            end
            %for qt_plus_one = 1:num_dicts
            %    p_qt_qt_plus_one(t, qt, qt_plus_one) = 1/sum(exp(log_alpha(t,:) + log_beta(t + 1,:) - log_alpha(t,qt) - log_beta(t+1, qt_plus_one) - P3_mat(t+1, qt_plus_one)))...
            %             * A(qt, qt_plus_one);
            %end
        end
        %p_qt_qt_plus_one(t, :, :) = p_qt_qt_plus_one(t, :, :) / sum(sum(p_qt_qt_plus_one(t, :, :)));
    end
    
    
    if any(any(any(isnan(p_qt_qt_plus_one))))
        '4'
        break
    end
    

    %tic()
    % P2_mat : P(zt, ft | qt (Vect)
    P2_mat = zeros(T, K, num_freq, num_dicts, 'double');
    %for ft = 1:num_freq
    %   for qt = 1:num_dicts
    %       temp = P5_mat(:, :, qt) .* repmat(P4_mat(ft, :, qt), T, 1);
    %       P2_mat(:, :, ft, qt) = temp ./ repmat(sum(temp, 2), 1, K);
    %   end
    %end
    
    for ft = 1:num_freq
       for qt = 1:num_dicts
           temp = log(P5_mat(:, :, qt)) + repmat(log(P4_mat(ft, :, qt)), T, 1);
           [max_val, max_ind] = max(temp');
           P2_mat(:, :, ft, qt) = temp - repmat(max_val' + log(sum(exp(temp - repmat(max_val', 1, K)), 2)), 1, K);
       end
    end 
    
    if any(any(any(any(isnan(P2_mat)))))
        '5'
        break
    end
    
    %tic()
    % P1_mat : P(zt, qt | ft, f_bold, v_bold)
    P1_mat = zeros(T, K, num_dicts, num_freq, 'double');
    for ft = 1:num_freq
       for qt = 1:num_dicts
           P1_mat(:, :, qt, ft) = repmat(p_qt(:, qt), 1, K) + P2_mat(:, :, ft, qt);
       end
    end
    
    if any(any(any(any(isnan(P1_mat)))))
        '6'
        break
    end

    %tic()
    % Update P4_mat
    %P4_mat = zeros(num_freq, K, num_dicts) ; % P(f|z,q)
    %for zt = 1:K
    %    for qt = 1:num_dicts
    %        for f = 1:num_freq
    %            P4_mat(f, zt, qt) = Vft(:, f)' * exp(P1_mat(:, zt, qt, f));
    %        end
    %        temp = log(P4_mat(:, zt, qt));
    %        [max_val, max_ind] = max(temp);
    %        temp_but_max = [temp(1:max_ind-1)', temp(max_ind+1:end)']; 
    %        P4_mat(:, zt, qt) = temp - max_val - log(1 + sum(exp(temp_but_max)));
    %    end
    %end
    
    % Update P4_mat
    P4_mat = zeros(num_freq, K, num_dicts) ; % P(f|z,q)
    for zt = 1:K
        for qt = 1:num_dicts
            for f = 1:num_freq
                P4_mat(f, zt, qt) = Vft(:, f)' * exp(P1_mat(:, zt, qt, f));
            end
            P4_mat(:, :, qt) = P4_mat(:, :, qt) + my_epsilon; % smoothing
            P4_mat(:, zt, qt) = P4_mat(:, zt, qt) / sum(P4_mat(:, zt, qt));
        end
    end
    
    
    if any(any(any(isnan(P4_mat))))
        '7'
        break
    end

    %tic()
    % Update P5_mat
    %P5_mat = zeros(T, K, num_dicts); % p(z|q)
    %for qt = 1:num_dicts
    %    for zt = 1:K
    %        P5_mat(:, zt, qt) = sum(Vft .* exp(squeeze(P1_mat(:, zt, qt, :))), 2);
    %    end
    %    temp = log(P5_mat(:, :, qt));
    %    [max_val, max_ind] = max(temp');
    %    P5_mat(:, :, qt) = temp - repmat(max_val' + log(sum(exp(temp - repmat(max_val', 1,K)), 2)), 1, K);
        
        %P5_mat(:, :, qt) = P5_mat(:, :, qt) + my_epsilon; % smoothing
        %P5_mat(:, :, qt) = P5_mat(:, :, qt) ./ repmat(sum(P5_mat(:, :, qt), 2), 1, K);
    %end
    
    P5_mat = zeros(T, K, num_dicts); % p(z|q)
    for qt = 1:num_dicts
        for zt = 1:K
            P5_mat(:, zt, qt) = sum(Vft .* squeeze(exp(P1_mat(:, zt, qt, :))), 2);
        end
        P5_mat(:, :, qt) = P5_mat(:, :, qt) + my_epsilon; % smoothing
        P5_mat(:, :, qt) = P5_mat(:, :, qt) ./ repmat(sum(P5_mat(:, :, qt), 2), 1, K);
    end
    
    
    if any(any(any(isnan(P5_mat))))
        '8'
        break
    end

    % Update mu_vq
    for q=1:num_dicts
        mu_vq(:) = (sum(repmat(v, 1, num_dicts) .* exp(p_qt)) ./ sum(exp(p_qt)))';
    end
    
    % Update sigma_vq
    for q=1:num_dicts
        sigma_vq(q) = sqrt(sum(exp(p_qt(:, q)).* (v - mu_vq(q)).^2) ./ sum(exp(p_qt(:, q))));
    end   
    
    % Update pi (start probability)
    for q=1:num_dicts
        temp = p_qt(1,:);
        [max_val, max_ind] = max(temp);
        temp_but_max = [temp(1:max_ind-1), temp(max_ind+1:end)];  
        pi(q) = temp(q) - max_val - log(1 + sum(exp(temp_but_max - max_val)));
    end
    
    % Update A (transition matrix)
    for q1=1:num_dicts
        temp1 = p_qt_qt_plus_one(:, q1,:);
        max_val1 = max(max(temp1));
        for q2=1:num_dicts
            temp2 = p_qt_qt_plus_one(:, q1,q2);
            max_val2 = max(temp2);
            A(q1,q2) = max_val2 - max_val1 + log(sum(exp(temp2 - max_val2))) - log(sum(sum(exp(temp1 - max_val1))));
        end
    end
    
    if any(any(isnan(A)))
        '10'
        break
    end
    
    [min(min(min(min(P1_mat)))), min(min(min(min(P2_mat)))), min(min(min(min(P3_mat)))), min(min(min(min(P4_mat)))), min(min(min(min(P5_mat))))]
    [sum(sum(sum(sum(P1_mat == 0)))), sum(sum(sum(sum(P2_mat== 0)))), sum(sum(sum(sum(P3_mat== 0)))), sum(sum(sum(sum(P4_mat == 0)))), sum(sum(sum(sum(P5_mat== 0))))]
    toc()   
    
    % Compute log likelihood
    yyy = 60;
    temp = log_alpha + log_beta;
    [max_val, max_ind] = max(temp(yyy,:));
    temp_but_max = [temp(1:max_ind-1), temp(max_ind+1:end)];   
    log_likelihood = max_val + log(1 + sum(exp(temp_but_max - max_val)));
    
    log_likelihoods = [log_likelihoods, log_likelihood];
    
    fprintf('Log likelihood: %d\n', log_likelihood)
    toc()

end

surf(exp(A));
view(0,90)
colormap gray
title('Transition Matrix A (between phonemes)', 'FontWeight', 'bold')

plot(log_likelihoods)
title('Evolution of Log-Likelihood', 'FontWeight', 'bold')
xlabel('Number of EM iterations')
ylabel('Log-Likelihood')