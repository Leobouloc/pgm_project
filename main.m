%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non Negative Hidden Markov Model of Audio with Source Separation

% NB : We tried to keep the same notations as in the original paper
%       Probabilites are numbered according to the number in the paper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

file_suffix = '1';

% Data (generate and save: 1, or load from save: 0)
generate = 1;
if generate
    
[Vft1, train_audio1, Vft2, train_audio2, ...
  Vft1_test, test_audio1, Vft2_test, test_audio2, ...
  Vft, audio_mixture, ...
  int_size, Vft_complex, sr, g, hop_size, window_size] = load_data(0);

    Vft1 = double(Vft1');
    Vft2 = double(Vft2');
    Vft1_test = double(Vft1_test');
    Vft2_test = double(Vft2_test');
    Vft = double(Vft');
    
    save('data_for_separation.mat', 'Vft1', 'train_audio1', 'Vft2', 'train_audio2', ...
                                  'Vft1_test', 'test_audio1', 'Vft2_test', 'test_audio2', ...
                                  'Vft', 'audio_mixture', ...
                                  'int_size', 'Vft_complex', 'sr', 'g', ...
                                  'hop_size', 'window_size')
else
    load('data_for_separation.mat')
end

% Pre-process Vft
Vft1 = preprocess_Vft(Vft1);
Vft2 = preprocess_Vft(Vft2);
Vft1_test = preprocess_Vft(Vft1_test);
Vft2_test = preprocess_Vft(Vft2_test);
Vft = preprocess_Vft(Vft);


% Set Parameters
K = 3; % Number of vectors per dictionary ; num latent components ; num phonems
num_dicts = 10; % Number of dictionaries

% For smoothing in single source EM
my_epsilon = 10^-10;

% Number of iterations for single source EM
num_iterations_single_source = 10;

% Display progression (Yes or No)
verbose = 1;


%%%%
compute_single_EM = 1   ;
if compute_single_EM
    % Actual EM single source here
    [P4_mat1, ~, mu_vq1, sigma_vq1, A1, Pi1, log_likelihoods1, T1, num_freq1] = EM_single_source(Vft1, K, num_dicts, num_iterations_single_source, my_epsilon, verbose);
    [P4_mat2, ~, mu_vq2, sigma_vq2, A2, Pi2, log_likelihoods2, T2, num_freq2] = EM_single_source(Vft2, K, num_dicts, num_iterations_single_source, my_epsilon, verbose);
    save('profiles.mat', 'P4_mat1', 'mu_vq1', 'sigma_vq1', 'A1', 'T1', 'num_freq1', 'P4_mat2', 'mu_vq2', 'sigma_vq2', 'A2',  'T2', 'num_freq2')
    clear Pi1 Pi2
else
    load('profiles.mat');
end
%%%%

fprintf('********\nNow Separating Sources\n************\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T1, num_freq1] = size(Vft1_test);
[T2, num_freq2] = size(Vft2_test);

assert(num_freq1 == num_freq2);
num_freq = num_freq1;

T = min(T1, T2);

Vft = Vft1_test(1:T, :) + Vft2_test(1:T, :);
v = sum(Vft, 2); % Observed energy

%%%%%

num_sources = 2;

% 1 - Initialise transition probabilites
Pi1 = log(1/num_dicts * ones(num_dicts, 1)); % Initial state
Pi2 = log(1/num_dicts * ones(num_dicts, 1)); % Initial state

% Dictionaries learned 
% Pu_mat = zeros(num_freq, K, num_sources, num_dicts); % P(ft|zt, st, qt)


% Initialise Transition probabilities
P11_mat = ones(T, K, num_sources, num_dicts, num_dicts) + rand(T, K, num_sources, num_dicts, num_dicts); % P(zt,st | q1, q2)
P11_mat = P11_mat ./ repmat(sum(sum(P11_mat, 2), 3), 1, K, num_sources);


% Compute log_likelihood
log_likelihoods = [];

num_iterations = 10;
for xxx = 1:num_iterations
    fprintf('*****************\n')
    tic()
    P10_mat = zeros(T, num_dicts, num_dicts); % P(ft_bold, vt | qt)
    % first term
    for q1 = 1:num_dicts
        for q2 = 1:num_dicts
           P10_mat(:, q1, q2) = log(p_vq1q2(v, q1, q2, mu_vq1, mu_vq2,sigma_vq1, sigma_vq2)); 
           for ft = 1:num_freq
               P10_mat(:, q1, q2) = P10_mat(:, q1, q2)' +  Vft(:, ft)' .* log((P4_mat1(ft, :, q1) * P11_mat(:,:, 1, q1, q2)' + P4_mat2(ft, :, q2) * P11_mat(:,:, 2, q1, q2)' ));
           end
        end
    end
    if any(any(any(any(any(any(P10_mat == -Inf))))))
        'log_beta has -Inf'
        assert(false)
    end    
    if any(any(any(any(any(any(isnan(P10_mat)))))))
        'log_beta has NaN'
        assert(false)
    end


    % 2.2 - Compute alpha, beta
    log_alpha = zeros(T, num_dicts, num_dicts);
    log_beta = zeros(T, num_dicts, num_dicts);

    log_alpha(1, :, :) = repmat(Pi1, 1, num_dicts) + repmat(Pi2', num_dicts, 1) + squeeze(P10_mat(1, :, :))';
    log_beta(T, :, :) = log(1);

    var_break = 0;
    for t = 2:T
        for q1 = 1:num_dicts
            for q2 = 1:num_dicts
                temp = repmat(A1(:, q1), 1, num_dicts) + repmat(A2(:, q2)', num_dicts, 1) + squeeze(log_alpha(t-1, :, :)); % WARNING MIGHT HAVE TO INVERSE TRANSPOSE
                temp = reshape(temp, num_dicts * num_dicts, 1)';
                [max_val, max_ind] = max(temp);
                temp_but_max = [temp(1:max_ind-1), temp(max_ind+1:end)];  

                log_alpha(t, q1, q2) = P10_mat(t, q1, q2) ...
                            + max_val + log(1 + sum(exp(temp_but_max - max_val)));
            end
        end
    end

    for t = 1:T-1
        for q1 = 1:num_dicts
            for q2 = 1:num_dicts
                temp = repmat(A1(q1, :)', 1, num_dicts) + repmat(A2(q2, :), num_dicts, 1) + squeeze(P10_mat(T-t+1, :, :)) + squeeze(log_beta(T-t+1, :, :)); 
                temp = reshape(temp, num_dicts * num_dicts, 1)';
                [max_val, max_ind] = max(temp);
                temp_but_max = [temp(1:max_ind-1), temp(max_ind+1:end)];
                log_beta(T-t, :, :) = max_val + log(1 + sum(exp(temp_but_max - max_val)));
            end
        end    
    end
    if any(any(any(any(any(any(isnan(log_alpha)))))))
        'log_alpha has NaN'
        assert(false)
    end    
    if any(any(any(any(any(any(isnan(log_beta)))))))
        'log_beta has NaN'
        assert(false)
    end


    % Reminder P4_mat : P(f|z,q) 

    % Compute transistion probabilities
    P9_mat = zeros(T, K, 2, num_freq, num_dicts, num_dicts); % P(zt, st| ft, q1, q2)
    %for q1=1:num_dicts
    %    for q2=1:num_dicts
    %        for f = 1:num_freq
    %            P9_mat(:, :, 1, f, q1, q2) = repmat(P4_mat1(f, :, q1), T, 1) .* P11_mat(:, :, 1, q1, q2);
    %            P9_mat(:, :, 2, f, q1, q2) = repmat(P4_mat2(f, :, q2), T, 1) .* P11_mat(:, :, 2, q1, q2);
    %            if any(any(P9_mat(:, :, 1, f, q1, q2) == 0))
    %               assert(false) 
    %            end            
    %        end
    %    end
    %end
    
    P9_part_1 = permute(repmat(P11_mat(:, :, 1, :, :), 1,1,1,1,1,num_freq), [1,2,6,4,5,3]);
    P9_part_2 = permute(repmat(P4_mat1, 1,1,1, T, num_dicts), [4,2,1,3,5]);
    P9_mat(:, :, 1, :, :, :) = P9_part_1 .* P9_part_2;
    
    P9_part_1 = permute(repmat(P11_mat(:, :, 2, :, :), 1,1,1,   1,1,num_freq), [1,2,6,4,5,3]);
    P9_part_2 = permute(repmat(P4_mat2, 1,1,1, T, num_dicts), [4,2,1,5,3]); % Inversion 5-3 for q1, q2
    P9_mat(:, :, 2, :, :, :) = P9_part_1 .* P9_part_2; 
    clear P9_part_1 P9_part_2

    P9_mat = P9_mat ./ repmat(sum(sum(P9_mat, 2), 3), 1, K, num_sources);

    if any(any(any(any(any(any(isnan(P9_mat)))))))
        'P9 has NaN'
        assert(false)
    end

    % Compute P8
    P8_mat = zeros(T, K, num_sources, num_dicts, num_dicts, num_freq); % P(zt, st, q1, q2 | ft)
    temp = log_alpha + log_beta;
    temp = reshape(temp, T, num_dicts * num_dicts);
    max_val = max(temp');

    P8_mat_part1 = log_alpha + log_beta - repmat((max_val' + log(sum(exp(temp - repmat(max_val', 1, num_dicts * num_dicts)), 2))), 1, num_dicts, num_dicts);
    P8_mat = permute(repmat(P8_mat_part1, 1, 1, 1, K, num_sources, num_freq), [1,4,5,6,2,3]) + log(P9_mat);
    P8_mat = permute(P8_mat, [1,2,3,5,6,4]);
    if any(any(any(any(any(any(isnan(P8_mat)))))))
        'P8 has NaN'
        assert(false)
    end

    % Update P11_mat
    P11_mat = sum(permute(repmat(Vft, 1,1, K, num_sources, num_dicts, num_dicts), [1,3,4,5,6,2]) .* exp(P8_mat), 6);
    P11_mat = P11_mat + my_epsilon;
    P11_mat = P11_mat ./ repmat(sum(sum(P11_mat, 2), 3), 1, K, num_sources);
    if any(any(any(any(any(any(isnan(P11_mat)))))))
        'P11 has NaN'
        assert(false)
    end

    % State probabities % p(q1, q2 | Observations)
    temp = log_alpha + log_beta;
    temp = reshape(temp, T, num_dicts * num_dicts);
    max_val = max(temp');
    p_q1_q2 = log_alpha + log_beta - repmat((max_val' + log(sum(exp(temp - repmat(max_val', 1, num_dicts * num_dicts)), 2))), 1, num_dicts, num_dicts);

    if any(any(any(any(any(any(isnan(p_q1_q2)))))))
        'p_q1_q2 has NaN'
        assert(false)
    end    
    
    % Update Pi1, Pi2
    temp = squeeze(p_q1_q2(1,:,:));
    max_val = max(temp, [], 2);
    Pi1 = max_val +  log(sum(exp(temp - repmat(max_val, 1, num_dicts)), 2));

    temp = squeeze(p_q1_q2(1,:,:))';
    max_val = max(temp, [], 2);
    Pi2 = max_val +  log(sum(exp(temp - repmat(max_val, 1, num_dicts)), 2));

    %temp = log_alpha + log_beta;
    %temp = reshape(temp, T, num_dicts * num_dicts);
    %max_val = max(temp');
    %log_likelihood = max_val' + log(sum(exp(temp - repmat(max_val', 1, num_dicts * num_dicts)), 2));
 
    yyy = 1;
    temp = log_alpha(yyy,:,:) + log_beta(yyy,:,:);
    temp = reshape(temp, 1, num_dicts * num_dicts);
    [max_val, max_ind] = max(temp);
    temp_but_max = [temp(1:max_ind-1), temp(max_ind+1:end)];
    log_likelihood = max_val + log(1 + sum(exp(temp_but_max - max_val)))    
    
    
    log_likelihoods = [log_likelihoods, log_likelihood(1)];
    fprintf('Log-likelihood at iteration %d is : %d\n', xxx, log_likelihood);
    toc()
end % end loop


%%%%%%%%%
% After the loop
%%%%%%%%%

p_final = zeros(T, num_sources, num_freq); % p(st|ft)

% s = 1
P_final_part_1 = permute(repmat(P11_mat(:, :, 1, :, :), 1,1,num_freq), [1,3,2,4,5]);
P_final_part_2 = permute(repmat(P4_mat1, 1, 1, 1, num_dicts, T), [5,1,2,3,4]);
P_final(:, 1, :) = sum(sum(permute(repmat(exp(p_q1_q2), 1, 1, 1, num_freq), [1, 4, 2, 3]) .* squeeze(sum(P_final_part_1 .* P_final_part_2, 3)), 4),3); 

% s = 2
P_final_part_1 = permute(repmat(P11_mat(:, :, 2, :, :), 1,1,num_freq), [1,3,2,4,5]);
P_final_part_2 = permute(repmat(P4_mat2, 1, 1, 1, num_dicts, T), [5,1,2,4,3]);
P_final(:, 2, :) = sum(sum(permute(repmat(exp(p_q1_q2), 1, 1, 1, num_freq), [1, 4, 2, 3]) .* squeeze(sum(P_final_part_1 .* P_final_part_2, 3)), 4),3); 
clear p_final_part_1 p_final_part_2

P_final = P_final ./ repmat(sum(P_final, 2), 1, 2);


save('my_separation.mat', 'Vft1_test', 'Vft2_test', 'Vft', 'P_final')
mixture_weights_to_audio(squeeze(P_final(1:T-1, 1, :)), Vft_complex, sr, g, hop_size, window_size)


if false
% Display dictionaries
display_dictionaries(P4_mat, 16000)

surf(exp(A));
view(0,90)
colormap gray
title('Transition Matrix A (between phonemes)', 'FontWeight', 'bold')

plot(log_likelihoods)
title('Evolution of Log-Likelihood', 'FontWeight', 'bold')
xlabel('Number of EM iterations')
ylabel('Log-Likelihood')
end