function Vft = preprocess_Vft(Vft)
    % Adding white noise at end to avoid zero probabilities
    Vft(end + 1, :) = 1;

    % Filling blancs by white noise for 0/ NaN issues with P5_mat
    Vft(all(Vft == 0, 2), :) = 1;
end