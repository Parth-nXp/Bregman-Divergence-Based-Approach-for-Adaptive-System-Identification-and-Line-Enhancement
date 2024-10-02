close all;
clear all;
clc;

% Generate unknown system and noise input
audio_input=audioread('bass.wav'); %input .wav signal
channel_taps=16; %optimal filter order
audio_input = audio_input + 0.05*(abs(min(audio_input)))*ones(size(audio_input)); % positive shift
noise=0.007*randn(size(audio_input));    %random noise
z = audio_input + noise;     % noisy signal

subplot(5,2,1)
plot(audio_input(1:1000)*10) %plot of original signal for 1000 samples
title('Original input signal');
xlabel('Frequency');
ylabel('Amplitude');
subplot(5,2,2)
plot(z(1:1000)*10) %plot of noisy signal for 1000 samples
title('Noisy time domain signal');
xlabel('Frequency');
ylabel('Amplitude');

signal_length=length(z); %length of noisy signal
delay=channel_taps; % introducing the delay to desired signal equal to channel length
delta=0.05; %step size(convergence gain) of adaptive filter
delayed_noisy_signal=[zeros(1,delay-1) z(delay:signal_length)']; %delaying noisy signal by the filter order
input = delayed_noisy_signal(1:441000);
desired_input = z(1:441000)';
input_length = length(input);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',65));  %data1
filter_weights = abs(randn(1,channel_taps)); %initialize filter coefficients
w_LMS = filter_weights;
w_KLLMS = filter_weights;
w_ISLMS = filter_weights;
w_AALMS = filter_weights;
w_BLMS = filter_weights;
filter_output_LMS = zeros(1,input_length); % Initialize outputs
filter_output_KLLMS = zeros(1,input_length); % Initialize outputs
filter_output_ISLMS = zeros(1,input_length); % Initialize outputs
filter_output_AALMS = zeros(1,input_length); % Initialize outputs
filter_output_BLMS = zeros(1,input_length); % Initialize outputs
alpha = 2;
beta = 2;
for i=channel_taps:input_length
    x1 = input(i:-1:i-channel_taps+1); % Select input for convolution
    filter_output_LMS(i) = w_LMS * x1'; % Convolve (multiply)weights with input
    filter_output_KLLMS(i) = w_KLLMS * x1'; % Convolve (multiply)weights with input
    filter_output_ISLMS(i) = w_ISLMS * x1'; % Convolve (multiply)weights with input
    filter_output_AALMS(i) = w_AALMS * x1'; % Convolve (multiply)weights with input
    filter_output_BLMS(i) = w_BLMS * x1'; % Convolve (multiply)weights with input
    e_LMS(i) = desired_input(i)-filter_output_LMS(i); % Calculate error
    e_KLLMS(i) = desired_input(i)-filter_output_KLLMS(i); % Calculate error
    e_ISLMS(i) = desired_input(i)-filter_output_ISLMS(i); % Calculate error
    e_AALMS(i) = desired_input(i)-filter_output_AALMS(i); % Calculate error
    e_BLMS(i) = desired_input(i)-filter_output_BLMS(i); % Calculate error
    w_LMS = w_LMS + delta*e_LMS(i)*x1;% Adjust weights
    w_KLLMS = w_LMS + delta*x1*(desired_input(i)/filter_output_KLLMS(i) - 1);
    w_ISLMS = w_ISLMS + delta * x1 * ((desired_input(i)/(filter_output_ISLMS(i)+10e-6)^2)-(1/(filter_output_ISLMS(i)+10e-6)));
    w_AALMS = w_AALMS + (delta/alpha) * x1*((desired_input(i)^alpha/(filter_output_AALMS(i))^2)-1);
    w_BLMS = w_BLMS + delta * x1 * (filter_output_BLMS(i)^(beta-1)*e_BLMS(i));
end
subplot(5,2,3)
plot(1:1000,filter_output_KLLMS(1:1000))       % Plot of filtered output
title('Filtered Signal');
xlabel('Frequency');
ylabel('Amplitude');

subplot(5,2,4)
plot(1:1000,e_KLLMS(1:1000))%plot of error signal
title('error signal'); 
xlabel('Frequency');
ylabel('Amplitude');

subplot(5,2,5)
plot(1:1000,filter_output_ISLMS(1:1000)/100+0.2)       % Plot of filtered output
title('Filtered Signal');
xlabel('Frequency');
ylabel('Amplitude');

subplot(5,2,6)
plot(1:1000,e_ISLMS(1:1000))%plot of error signal
title('error signal'); 
xlabel('Frequency');
ylabel('Amplitude');

subplot(5,2,7)
plot(1:1000,filter_output_AALMS(1:1000))       % Plot of filtered output
title('Filtered Signal');
xlabel('Frequency');
ylabel('Amplitude');

subplot(5,2,8)
plot(1:1000,e_AALMS(1:1000))%plot of error signal
title('error signal'); 
xlabel('Frequency');
ylabel('Amplitude');

subplot(5,2,9)
plot(1:1000,filter_output_BLMS(1:1000))       % Plot of filtered output
% title('Filtered Signal');
xlabel('Frequency');
ylabel('Amplitude');

subplot(5,2,10)
plot(1:1000,e_BLMS(1:1000))%plot of error signal
title('error signal'); 
xlabel('Frequency');
ylabel('Amplitude');