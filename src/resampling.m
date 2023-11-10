% Original signal
original_signal = linspace(0, 10, 10);

% Time vector for the original signal
original_time = linspace(0, 1, 10);

% New time vector for the desired signal
desired_time = linspace(-1, 1, 10);

% Interpolate the original signal to match the new time vector
desired_signal = interp1(original_time, original_signal, desired_time);

% Plotting
figure;
subplot(2,1,1);
stem(original_time, original_signal, 'LineWidth', 2);
title('Original Signal');
xlabel('Time');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
stem(desired_time, desired_signal, 'LineWidth', 2);
title('Desired Signal');
xlabel('Time');
ylabel('Amplitude');
grid on;
