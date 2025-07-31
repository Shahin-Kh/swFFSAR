function peakiness = calculate_peakiness_new(signal)
    % Ensure the signal is a column vector
    signal = signal(:);

    % Ensure all values are non-negative, if W_i represents magnitude or power
    % This line is optional depending on your interpretation of W_i
    signal = abs(signal);

    % Calculate numerator (maximum value)
    max_val = max(signal);

    % Calculate denominator (sum of all values)
    sum_val = sum(signal);

    % Prevent division by zero
    if sum_val == 0
        peakiness = NaN;
        warning('Sum of signal is zero. Peakiness is undefined.');
        return;
    end

    % Peakiness calculation
    peakiness = (max_val / sum_val)*100;
end
