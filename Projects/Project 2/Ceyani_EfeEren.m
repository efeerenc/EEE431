%% Part I: Signal Spaces
signal_data = load("signaldata2.mat").x;

Fs = 1200;

time1 = [0:1/Fs:1];
a = zeros([1, 10]);
f = zeros([1, 10]);

for k = 1:10
    energy1 = sum(sin(2*pi*(30*k)).^2);
    energy2 = sum(sin(2*pi*(30*k-10)).^2);
    energy3 = sum(sin(2*pi*(30*k-20)).^2);
    
    basis1 = sqrt(1/energy1)*sin(2*pi*(30*k).*time1);
    basis2 = sqrt(1/energy2)*sin(2*pi*(30*k-10).*time1);
    basis3 = sqrt(1/energy3)*sin(2*pi*(30*k-20).*time1);
    
    receiver1 = sum(basis1.*signal_data);
    receiver2 = sum(basis2.*signal_data);
    receiver3 = sum(basis3.*signal_data);
    
    if max([receiver1, receiver2, receiver3]) == receiver1
        a(k) = sqrt(energy1)*receiver1/600;
        f(k) = 30*k;
    elseif max([receiver1, receiver2, receiver3]) == receiver2
        a(k) = sqrt(energy2)*receiver2/600;
        f(k) = 30*k-10;
    else
        a(k) = sqrt(energy3)*receiver3/600;
        f(k) = 30*k-20;
    end 
end
%Frequencies are as the following: 30, 50, 80, 120, 150, 160, 190, 230,
%250, 300
%recon
% y = zeros([1,1201]);
% for l = 1:10
%     y = y + a(l)*sin(2*pi*f(l).*time1);
% end
% figure;
% plot(time1, signal_data);
% figure;
% plot(time1, y);

%d
gauss_noise1 = randn(size(signal_data))*1; % sigma^2 = 1
gauss_noise2 = randn(size(signal_data))*5; % sigma^2 = 25
gauss_noise3 = randn(size(signal_data))*10; % sigma^2 = 100

figure;
plot(time1, signal_data);
hold on;
plot(time1, signal_data + gauss_noise1);
xlabel("Time (s)");
title("Original signal with zero-mean Gaussian noise of \sigma^2 = 1");
legend(["x(t)","\sigma^2=1"]);

figure;
plot(time1, signal_data);
hold on;
plot(time1, signal_data + gauss_noise2);
xlabel("Time (s)");
title("Original signal with zero-mean Gaussian noise of \sigma^2 = 25");
legend(["x(t)","\sigma^2=25"]);

figure;
plot(time1, signal_data);
hold on;
plot(time1, signal_data + gauss_noise3);
xlabel("Time (s)");
title("Original signal with zero-mean Gaussian noise of \sigma^2 = 100");
legend(["x(t)","\sigma^2=100"]);

%e
a_1 = zeros([1,10]);
a_25 = zeros([1,10]);
a_100 = zeros([1,10]);

f_1 = zeros([1,10]);
f_25 = zeros([1,10]);
f_100 = zeros([1,10]);

for k = 1:10
    energy1 = sum(sin(2*pi*(30*k)).^2);
    energy2 = sum(sin(2*pi*(30*k-10)).^2);
    energy3 = sum(sin(2*pi*(30*k-20)).^2);
    
    basis1 = sqrt(1/energy1)*sin(2*pi*(30*k).*time1);
    basis2 = sqrt(1/energy2)*sin(2*pi*(30*k-10).*time1);
    basis3 = sqrt(1/energy3)*sin(2*pi*(30*k-20).*time1);
    
    receiver1 = sum(basis1.*(signal_data + gauss_noise1));
    receiver2 = sum(basis2.*(signal_data + gauss_noise1));
    receiver3 = sum(basis3.*(signal_data + gauss_noise1));
    
    if max([receiver1, receiver2, receiver3]) == receiver1
        a_1(k) = sqrt(energy1)*receiver1/600;
        f_1(k) = 30*k;
    elseif max([receiver1, receiver2, receiver3]) == receiver2
        a_1(k) = sqrt(energy2)*receiver2/600;
        f_1(k) = 30*k-10;
    else
        a_1(k) = sqrt(energy3)*receiver3/600;
        f_1(k) = 30*k-20;
    end 
end

for k = 1:10
    energy1 = sum(sin(2*pi*(30*k)).^2);
    energy2 = sum(sin(2*pi*(30*k-10)).^2);
    energy3 = sum(sin(2*pi*(30*k-20)).^2);
    
    basis1 = sqrt(1/energy1)*sin(2*pi*(30*k).*time1);
    basis2 = sqrt(1/energy2)*sin(2*pi*(30*k-10).*time1);
    basis3 = sqrt(1/energy3)*sin(2*pi*(30*k-20).*time1);
    
    receiver1 = sum(basis1.*(signal_data + gauss_noise2));
    receiver2 = sum(basis2.*(signal_data + gauss_noise2));
    receiver3 = sum(basis3.*(signal_data + gauss_noise2));
    
    if max([receiver1, receiver2, receiver3]) == receiver1
        a_25(k) = sqrt(energy1)*receiver1/600;
        f_25(k) = 30*k;
    elseif max([receiver1, receiver2, receiver3]) == receiver2
        a_25(k) = sqrt(energy2)*receiver2/600;
        f_25(k) = 30*k-10;
    else
        a_25(k) = sqrt(energy3)*receiver3/600;
        f_25(k) = 30*k-20;
    end 
end

for k = 1:10
    energy1 = sum(sin(2*pi*(30*k)).^2);
    energy2 = sum(sin(2*pi*(30*k-10)).^2);
    energy3 = sum(sin(2*pi*(30*k-20)).^2);
    
    basis1 = sqrt(1/energy1)*sin(2*pi*(30*k).*time1);
    basis2 = sqrt(1/energy2)*sin(2*pi*(30*k-10).*time1);
    basis3 = sqrt(1/energy3)*sin(2*pi*(30*k-20).*time1);
    
    receiver1 = sum(basis1.*(signal_data + gauss_noise3));
    receiver2 = sum(basis2.*(signal_data + gauss_noise3));
    receiver3 = sum(basis3.*(signal_data + gauss_noise3));
    
    if max([receiver1, receiver2, receiver3]) == receiver1
        a_100(k) = sqrt(energy1)*receiver1/600;
        f_100(k) = 30*k;
    elseif max([receiver1, receiver2, receiver3]) == receiver2
        a_100(k) = sqrt(energy2)*receiver2/600;
        f_100(k) = 30*k-10;
    else
        a_100(k) = sqrt(energy3)*receiver3/600;
        f_100(k) = 30*k-20;
    end 
end


%% Part II: Binary Modulation
%i=1
symbols_binary = randi([0, 1], 1, 5);
x_a = [];
F_sampling = 1000;
time_binary = [0:1/F_sampling:0.1-(1/F_sampling)];

signal_vector = (time_binary <= 0.025).*time_binary + ((0.025 < time_binary)&(time_binary <= 0.075)).*(0.05-time_binary) + ((0.075 < time_binary)&(time_binary <= 0.1)).*(time_binary-0.1);

for i = 1:5
    if symbols_binary(i) == 0
        x_a = [x_a, -signal_vector];
    else
        x_a = [x_a, signal_vector];
    end
end
time_elapsed = [0:1/(F_sampling):0.5-1/(F_sampling)];
figure;
plot(time_elapsed, x_a);
xlabel("Time (sec)");
ylabel("x(t)");
title("Signal generated from the bits");

%b
energy_binary = sum(signal_vector.^2);

basis_binary = signal_vector/sqrt(energy_binary);
figure;
plot(time_binary, basis_binary);
xlabel("Time (sec)");
ylabel("\delta(t)");
title("The basis vector of the signal space");

repr1 = sum(basis_binary.*signal_vector);
repr2 = sum(basis_binary.*(-signal_vector));

figure;
plot(time_binary, repr1*basis_binary);
xlabel("Time (sec)");
ylabel("-s(t)");
title("Negative signal represented in the signal space");

figure;
plot(time_binary, repr2*basis_binary);
xlabel("Time (sec)");
ylabel("s(t)");
title("Positive signal represented in the signal space");

%c

gauss_binary1 = randn(size(signal_vector))*0.01;
gauss_binary2 = randn(size(signal_vector))*0.1;
gauss_binary3 = randn(size(signal_vector))*1;

figure;
plot(time_binary, signal_vector);
hold on;
plot(time_binary, signal_vector + gauss_binary1);
xlabel("Time (sec)");
title("Original signal with zero-mean Gaussian noise with \sigma^2 = 10^{-4}");

figure;
plot(time_binary, signal_vector);
hold on;
plot(time_binary, signal_vector + gauss_binary2);
xlabel("Time (sec)");
title("Original signal with zero-mean Gaussian noise with \sigma^2 = 10^{-2}");

figure;
plot(time_binary, signal_vector);
hold on;
plot(time_binary, signal_vector + gauss_binary3);
xlabel("Time (sec)");
title("Original signal with zero-mean Gaussian noise with \sigma^2 = 1");

%e
%Iteration 1
bits_binary1 = string(randi([0,1],1,10^5));

%Generate signal
x_binary1 = zeros(1, 100*10^5);
loc_binary0 = find(bits_binary1 == "0");
loc_binary1 = find(bits_binary1 == "1");

for i = 1:size(loc_binary0, 2)
    x_binary1(100*(loc_binary0(i)-1)+1:100*loc_binary0(i)) = -signal_vector;
end

for i = 1:size(loc_binary1, 2)
    x_binary1(100*(loc_binary1(i)-1)+1:100*loc_binary1(i)) = signal_vector;
end

sigma1 = 1;
gauss_binary1 = randn(size(x_binary1))*sigma1;
perturbed_binary1 = gauss_binary1 + x_binary1;
snrdb_pb1 = 20*log10(energy_binary/(2*(sigma1^2)));

%Predict message
r1_binary = zeros(1, 10^5);
for i = 1:10^5
    r1_binary(i) = sum(perturbed_binary1(100*(i-1)+1:100*i).*basis_binary);
end

bit0 = r1_binary < 0;
bit1 = r1_binary >= 0;

loc_bit0 = find(bit0 == 1);
loc_bit1 = find(bit1 == 1);

predictions_binary1 = string(zeros(1,10^5));

for i = 1:size(loc_bit0, 2)
    predictions_binary1(loc_bit0(i)) = "0";
end

for i = 1:size(loc_bit1, 2)
    predictions_binary1(loc_bit1(i)) = "1";
end

prob_error_binary1 = 1 - sum(predictions_binary1 == bits_binary1)/(10^5);

%Iteration 2
bits_binary2 = string(randi([0,1],1,10^5));

%Generate signal
x_binary2 = zeros(1, 100*10^5);
loc_binary0 = find(bits_binary2 == "0");
loc_binary1 = find(bits_binary2 == "1");

for i = 1:size(loc_binary0, 2)
    x_binary2(100*(loc_binary0(i)-1)+1:100*loc_binary0(i)) = -signal_vector;
end

for i = 1:size(loc_binary1, 2)
    x_binary2(100*(loc_binary1(i)-1)+1:100*loc_binary1(i)) = signal_vector;
end

sigma2 = 0.1;
gauss_binary2 = randn(size(x_binary2))*sigma2;
perturbed_binary2 = gauss_binary2 + x_binary2;
snrdb_pb2 = 20*log10(energy_binary/(2*(sigma2^2)));

%Predict message
r1_binary = zeros(1, 10^5);
for i = 1:10^5
    r1_binary(i) = sum(perturbed_binary2(100*(i-1)+1:100*i).*basis_binary);
end

bit0 = r1_binary < 0;
bit1 = r1_binary >= 0;

loc_bit0 = find(bit0 == 1);
loc_bit1 = find(bit1 == 1);

predictions_binary2 = string(zeros(1,10^5));

for i = 1:size(loc_bit0, 2)
    predictions_binary2(loc_bit0(i)) = "0";
end

for i = 1:size(loc_bit1, 2)
    predictions_binary2(loc_bit1(i)) = "1";
end

prob_error_binary2 = 1 - sum(predictions_binary2 == bits_binary2)/(10^5);

%Iteration 3
bits_binary3 = string(randi([0,1],1,10^5));

%Generate signal
x_binary3 = zeros(1, 100*10^5);
loc_binary0 = find(bits_binary3 == "0");
loc_binary1 = find(bits_binary3 == "1");

for i = 1:size(loc_binary0, 2)
    x_binary3(100*(loc_binary0(i)-1)+1:100*loc_binary0(i)) = -signal_vector;
end

for i = 1:size(loc_binary1, 2)
    x_binary3(100*(loc_binary1(i)-1)+1:100*loc_binary1(i)) = signal_vector;
end

sigma3 = 0.04;
gauss_binary3 = randn(size(x_binary3))*sigma3;
perturbed_binary3 = gauss_binary3 + x_binary3;
snrdb_pb3 = 20*log10(energy_binary/(2*(sigma3^2)));

%Predict message
r1_binary = zeros(1, 10^5);
for i = 1:10^5
    r1_binary(i) = sum(perturbed_binary3(100*(i-1)+1:100*i).*basis_binary);
end

bit0 = r1_binary < 0;
bit1 = r1_binary >= 0;

loc_bit0 = find(bit0 == 1);
loc_bit1 = find(bit1 == 1);

predictions_binary3 = string(zeros(1,10^5));

for i = 1:size(loc_bit0, 2)
    predictions_binary3(loc_bit0(i)) = "0";
end

for i = 1:size(loc_bit1, 2)
    predictions_binary3(loc_bit1(i)) = "1";
end

prob_error_binary3 = 1 - sum(predictions_binary3 == bits_binary3)/(10^5);

%Theoretical part
sigmas_binary = [2:-0.0001:0.0015];
snrs = energy_binary./((F_sampling/2)*2*(sigmas_binary.^2));
q_error = 0.5*erfc(((2*snrs).^(0.5))./(sqrt(2)));

figure;
semilogy(20*log10(snrs), q_error);
xlabel("E_s/2\sigma^2 (dB)");
ylabel("P_e");
title("Theoretical probability of symbol error vs. SNR");
hold on;
plot(snrdb_pb1, prob_error_binary1, "r*");
plot(snrdb_pb2, prob_error_binary2, "k*");
plot(snrdb_pb3, prob_error_binary3, "g*");
legend(["Theoretical", "\sigma = 1.5", "\sigma = 0.1", "\sigma = 0.04"]);

%g
%Iteration 1
bits_unequal1 = string(ones(1,10^5).*(rand(1,10^5)>0.9));

%Generate signal
x_binary_unequal1 = zeros(1, 100*10^5);
loc_binary0 = find(bits_unequal1 == "0");
loc_binary1 = find(bits_unequal1 == "1");

for i = 1:size(loc_binary0, 2)
    x_binary_unequal1(100*(loc_binary0(i)-1)+1:100*loc_binary0(i)) = -signal_vector;
end

for i = 1:size(loc_binary1, 2)
    x_binary_unequal1(100*(loc_binary1(i)-1)+1:100*loc_binary1(i)) = signal_vector;
end

sigma1_unequal = 1;
gauss_binary1_unequal = randn(size(x_binary_unequal1))*sigma1_unequal;
perturbed_binary1_unequal = gauss_binary1_unequal + x_binary_unequal1;
snrdb_pb1_unequal = 20*log10(energy_binary/(2*(sigma1_unequal^2)));

%Predict message
r1_binary_unequal = zeros(1, 10^5);
for i = 1:10^5
    r1_binary_unequal(i) = sum(perturbed_binary1_unequal(100*(i-1)+1:100*i).*basis_binary);
end

bit0 = r1_binary_unequal < 0;
bit1 = r1_binary_unequal >= 0;

bit0map = r1_binary < log(9)*(2*sigma1_unequal^2)/(4*sqrt(energy_binary));
bit1map = r1_binary >= log(9)*(2*(sigma1_unequal^2))/(4*sqrt(energy_binary));

loc_bit0 = find(bit0 == 1);
loc_bit1 = find(bit1 == 1);
loc_bit0map = find(bit0map == 1);
loc_bit1map = find(bit1map == 1);

predictions_binary_unequal1 = string(zeros(1,10^5));
predictions_binary_unequal1_map = string(zeros(1,10^5));

for i = 1:size(loc_bit0, 2)
    predictions_binary_unequal1(loc_bit0(i)) = "0";
end

for i = 1:size(loc_bit1, 2)
    predictions_binary_unequal1(loc_bit1(i)) = "1";
end

for i = 1:size(loc_bit0map, 2)
    predictions_binary_unequal1_map(loc_bit0map(i)) = "0";
end

for i = 1:size(loc_bit1map, 2)
    predictions_binary_unequal1_map(loc_bit1map(i)) = "1";
end

prob_error_binary1_unequal = 1 - sum(predictions_binary_unequal1 == bits_unequal1)/(10^5);
prob_error_binary1_unequal_map = 1 - sum(predictions_binary_unequal1_map == bits_unequal1)/(10^5);

%Iteration 2
bits_unequal2 = string(ones(1,10^5).*(rand(1,10^5)>0.9));

%Generate signal
x_binary_unequal2 = zeros(1, 100*10^5);
loc_binary0 = find(bits_unequal2 == "0");
loc_binary1 = find(bits_unequal2 == "1");

for i = 1:size(loc_binary0, 2)
    x_binary_unequal2(100*(loc_binary0(i)-1)+1:100*loc_binary0(i)) = -signal_vector;
end

for i = 1:size(loc_binary1, 2)
    x_binary_unequal2(100*(loc_binary1(i)-1)+1:100*loc_binary1(i)) = signal_vector;
end

sigma2_unequal = 0.1;
gauss_binary2_unequal = randn(size(x_binary_unequal2))*sigma2;
perturbed_binary2_unequal = gauss_binary2_unequal + x_binary_unequal2;
snrdb_pb2_unequal = 20*log10(energy_binary/(2*(sigma2_unequal^2)));

%Predict message
r1_binary = zeros(1, 10^5);
for i = 1:10^5
    r1_binary(i) = sum(perturbed_binary2_unequal(100*(i-1)+1:100*i).*basis_binary);
end

bit0 = r1_binary < 0;
bit1 = r1_binary >= 0;

bit0map = r1_binary < log(9)*(2*sigma2_unequal^2)/(4*sqrt(energy_binary));
bit1map = r1_binary >= log(9)*(2*(sigma2_unequal^2))/(4*sqrt(energy_binary));


loc_bit0 = find(bit0 == 1);
loc_bit1 = find(bit1 == 1);
loc_bit0map = find(bit0map == 1);
loc_bit1map = find(bit1map == 1);

predictions_binary_unequal2 = string(zeros(1,10^5));
predictions_binary_unequal2_map = string(zeros(1,10^5));

for i = 1:size(loc_bit0, 2)
    predictions_binary_unequal2(loc_bit0(i)) = "0";
end

for i = 1:size(loc_bit1, 2)
    predictions_binary_unequal2(loc_bit1(i)) = "1";
end

for i = 1:size(loc_bit0map, 2)
    predictions_binary_unequal2_map(loc_bit0map(i)) = "0";
end

for i = 1:size(loc_bit1map, 2)
    predictions_binary_unequal2_map(loc_bit1map(i)) = "1";
end

prob_error_binary2_unequal = 1 - sum(predictions_binary_unequal2 == bits_unequal2)/(10^5);
prob_error_binary2_unequal_map = 1 - sum(predictions_binary_unequal2_map == bits_unequal2)/(10^5);

%Iteration 3
zero_dist = 0.9;
bits_unequal3 = string(ones(1,10^5).*(rand(1,10^5)>0.9));

%Generate signal
x_binary_unequal3 = zeros(1, 100*10^5);
loc_binary0 = find(bits_unequal3 == "0");
loc_binary1 = find(bits_unequal3 == "1");

for i = 1:size(loc_binary0, 2)
    x_binary_unequal3(100*(loc_binary0(i)-1)+1:100*loc_binary0(i)) = -signal_vector;
end

for i = 1:size(loc_binary1, 2)
    x_binary_unequal3(100*(loc_binary1(i)-1)+1:100*loc_binary1(i)) = signal_vector;
end

sigma3_unequal = 0.04;
gauss_binary3 = randn(size(x_binary_unequal3))*sigma3_unequal;
perturbed_binary3_unequal = gauss_binary3 + x_binary_unequal3;
snrdb_pb3_unequal = 20*log10(energy_binary/(2*(sigma3_unequal^2)));

%Predict message
r1_binary = zeros(1, 10^5);
for i = 1:10^5
    r1_binary(i) = sum(perturbed_binary3_unequal(100*(i-1)+1:100*i).*basis_binary);
end

bit0 = r1_binary < 0;
bit1 = r1_binary >= 0;

bit0map = r1_binary < log(9)*(2*(sigma3_unequal^2))/(4*sqrt(energy_binary));
bit1map = r1_binary >= log(9)*(2*(sigma3_unequal^2))/(4*sqrt(energy_binary));

loc_bit0 = find(bit0 == 1);
loc_bit1 = find(bit1 == 1);
loc_bit0map = find(bit0map == 1);
loc_bit1map = find(bit1map == 1);

predictions_binary_unequal3 = string(zeros(1,10^5));
predictions_binary_unequal3_map = string(zeros(1,10^5));

for i = 1:size(loc_bit0, 2)
    predictions_binary_unequal3(loc_bit0(i)) = "0";
end

for i = 1:size(loc_bit1, 2)
    predictions_binary_unequal3(loc_bit1(i)) = "1";
end

for i = 1:size(loc_bit0map, 2)
    predictions_binary_unequal3_map(loc_bit0map(i)) = "0";
end

for i = 1:size(loc_bit1map, 2)
    predictions_binary_unequal3_map(loc_bit1map(i)) = "1";
end

prob_error_binary3_unequal = 1 - sum(predictions_binary_unequal3 == bits_unequal3)/(10^5);
prob_error_binary3_unequal_map = 1 - sum(predictions_binary_unequal3_map == bits_unequal3)/(10^5);

%Theoretical part
sigmas_binary_unequal = [2:-0.0001:0.0015];
snrs_unequal = energy_binary./((F_sampling/2)*2*(sigmas_binary_unequal.^2));
q_error_unequal = 0.5*erfc(((2*snrs_unequal).^(0.5))./(sqrt(2)));

figure;
semilogy(20*log10(snrs_unequal), q_error_unequal);
xlabel("E_s/2\sigma^2 (dB)");
ylabel("P_e");
title("Theoretical probability of symbol error vs. SNR (part g)");
hold on;
plot(snrdb_pb1_unequal, prob_error_binary1_unequal, "r*");
plot(snrdb_pb2_unequal, prob_error_binary2_unequal, "k*");
plot(snrdb_pb3_unequal, prob_error_binary3_unequal, "g*");
plot(snrdb_pb1_unequal, prob_error_binary1_unequal_map, "ro");
plot(snrdb_pb2_unequal, prob_error_binary2_unequal_map, "ko");
plot(snrdb_pb3_unequal, prob_error_binary3_unequal_map, "go");
legend(["Theoretical (ML)", "\sigma = 1.5", "\sigma = 0.1", "\sigma = 0.04"]);

%h
alphas = [0.1:0.005:0.5];
snrs_parth = zeros(1,size(alphas, 2));
prob_error_ml = zeros(1,size(alphas, 2));
prob_error_map = zeros(1,size(alphas, 2));

for i = 1:size(alphas, 2)

    bits_part_h1 = string(ones(1,10^5).*(rand(1,10^5)>(1-alphas(i))));

    %Generate signal
    x_part_h1 = zeros(1, 100*10^5);
    loc_binary0 = find(bits_part_h1 == "0");
    loc_binary1 = find(bits_part_h1 == "1");

    for j = 1:size(loc_binary0, 2)
        x_part_h1(100*(loc_binary0(j)-1)+1:100*loc_binary0(j)) = -signal_vector;
    end

    for j = 1:size(loc_binary1, 2)
        x_part_h1(100*(loc_binary1(j)-1)+1:100*loc_binary1(j)) = signal_vector;
    end

    sigma_parth = 0.1;
    gauss_part_h1 = randn(size(x_part_h1))*sigma_parth;
    perturbed_part_h1 = gauss_part_h1 + x_part_h1;
    snrdb_part_h1 = 20*log10(energy_binary/(2*(sigma_parth^2)));
    snrs_parth(i) = snrdb_part_h1;

    %Predict message
    r1_binary = zeros(1, 10^5);
    for j = 1:10^5
        r1_binary(j) = sum(perturbed_part_h1(100*(j-1)+1:100*j).*basis_binary);
    end

    bit0 = r1_binary < 0;
    bit1 = r1_binary >= 0;

    bit0map = r1_binary < log((1-alphas(i))/alphas(i))*(2*(sigma_parth^2))/(4*sqrt(energy_binary));
    bit1map = r1_binary >= log((1-alphas(i))/alphas(i))*(2*(sigma_parth^2))/(4*sqrt(energy_binary));

    loc_bit0 = find(bit0 == 1);
    loc_bit1 = find(bit1 == 1);
    loc_bit0map = find(bit0map == 1);
    loc_bit1map = find(bit1map == 1);

    predictions_part_h1 = string(zeros(1,10^5));
    predictions_part_h1_map = string(zeros(1,10^5));

    for j = 1:size(loc_bit0, 2)
        predictions_part_h1(loc_bit0(j)) = "0";
    end

    for j = 1:size(loc_bit1, 2)
        predictions_part_h1(loc_bit1(j)) = "1";
    end

    for j = 1:size(loc_bit0map, 2)
        predictions_part_h1_map(loc_bit0map(j)) = "0";
    end

    for j = 1:size(loc_bit1map, 2)
        predictions_part_h1_map(loc_bit1map(j)) = "1";
    end

    prob_error_part_h1 = 1 - sum(predictions_part_h1 == bits_part_h1)/(10^5);
    prob_error_part_h1_map = 1 - sum(predictions_part_h1_map == bits_part_h1)/(10^5);
    
    prob_error_ml(i) = prob_error_part_h1;
    prob_error_map(i) = prob_error_part_h1_map;
    
end

figure;
xlabel("\alpha");
ylabel("P_e");
title("Probability of error comparison between ML and MAP rule for fixed SNR with ranging \alpha");
hold on;
plot(alphas, prob_error_ml);
plot(alphas, prob_error_map);
legend(["ML Rule","MAP Rule"]);

%% Part III: Frequency Shift Keying (FSK)
%i = rem(21903359, 2) = 1
%f_i = f_1 = 250Hz

%Parameters
T = 0.1;
F_s = 5000;
time = [0:1/F_s:T-(1/F_s)];
f1 = 250;

%Function g(t)
g = time == time;

%Signals
s0 = cos(2*pi*f1.*time).*g;
s1 = -s0;
s2 = cos(2*pi*2*f1.*time).*g;
s3 = -s2;
energy = sum(s0.*s0);

%Part a
%Random symbols
symbols = [];

for i = 1:10
    if rem(i,2) == 0
        bit2 = randi([0, 1]);
        symbols = [symbols, strcat(string(bit), string(bit2))];
    end
    bit = randi([0,1]);
end

%Generate signal from symbols
x = [];

for i = 1:length(symbols)
    if symbols(i) == "00"
        x = [x, s0];
    elseif symbols(i) == "01"
        x = [x, s1];
    elseif symbols(i) == "10"
        x = [x, s2];
    else
        x = [x, s3];
    end
end

extended_time = [0:1/(F_s):5*T-(1/F_s)];
figure;
plot(extended_time, x);
xlabel("Time (s)");
ylabel("x(t)");
title("Signal generated from the symbols");

%Part b
psi1 = sqrt(1/energy)*cos(500*pi.*time);
psi2 = sqrt(1/energy)*cos(1000*pi.*time);

figure;
plot(time, psi1);
xlabel("Time (s)");
ylabel("\Psi_1(t)");
title("First basis function");
figure;
plot(time, psi2);
xlabel("Time (s)");
ylabel("\Psi_2(t)");
title("Second basis function");

projection0 = [sum(psi1.*s0), sum(psi2.*s0)];
projection1 = [sum(psi1.*s1), sum(psi2.*s1)];
projection2 = [sum(psi1.*s2), sum(psi2.*s2)];
projection3 = [sum(psi1.*s3), sum(psi2.*s3)];

representation0 = projection0(1)*psi1 + projection0(2)*psi2;
representation1 = projection1(1)*psi1 + projection1(2)*psi2;
representation2 = projection2(1)*psi1 + projection2(2)*psi2;
representation3 = projection3(1)*psi1 + projection3(2)*psi2;

figure;
plot(time, representation0);
xlabel("Time (s)");
ylabel("s_1(t)");
title("First signal's representation on the signal space");
figure;
plot(time, representation1);
xlabel("Time (s)");
ylabel("s_2(t)");
title("Second signal's representation on the signal space");
figure;
plot(time, representation2);
xlabel("Time (s)");
ylabel("s_3(t)");
title("Third signal's representation on the signal space");
figure;
plot(time, representation3);
xlabel("Time (s)");
ylabel("s_4(t)");
title("Fourth signal's representation on the signal space");

%Part c
gaussian1 = randn(size(x))*(1/10);
gaussian2 = randn(size(x));
gaussian3 = randn(size(x))*(10);

perturbed1 = gaussian1 + x;
perturbed2 = gaussian2 + x;
perturbed3 = gaussian3 + x;

figure;
plot(extended_time, x);
hold on;
plot(extended_time, perturbed1);
xlabel("Time (s)");
title("x(t) with zero-mean white Gaussian of \sigma^2 = 10^{-2}");
legend(["x(t)","\sigma^2 = 10^{-2}"]);

figure;
plot(extended_time, x);
hold on;
plot(extended_time, perturbed2);
xlabel("Time (s)");
title("x(t) with zero-mean white Gaussian of \sigma^2 = 1");
legend(["x(t)","\sigma^2 = 1"]);

figure;
plot(extended_time, x);
hold on;
plot(extended_time, perturbed3);
xlabel("Time (s)");
title("x(t) with zero-mean white Gaussian of \sigma^2 = 10^{2}");
legend(["x(t)","\sigma^2 = 10^2"]);

%Part e


%Iteration 1
bits1 = randi([0,1],1,2*10^5);
bits1 = string(reshape(bits1, [], 2));
symbols1 = reshape(strcat(bits1(:,1), bits1(:,2)), [1, 100000]);

%Generate signal
x1 = zeros(1,500*10^5);
loc0 = find(symbols1 == "00");
loc1 = find(symbols1 == "01");
loc2 = find(symbols1 == "10");
loc3 = find(symbols1 == "11");

for i = 1:size(loc0, 2)
    x1(500*(loc0(i)-1)+1:500*loc0(i)) = s0;
end

for i = 1:size(loc1, 2)
    x1(500*(loc1(i)-1)+1:500*loc1(i)) = s1;
end

for i = 1:size(loc2, 2)
    x1(500*(loc2(i)-1)+1:500*loc2(i)) = s2;
end

for i = 1:size(loc3, 2)
    x1(500*(loc3(i)-1)+1:500*loc3(i)) = s3;
end

sigma_x1 = 10;
noise_x1 = randn(size(x1))*(sigma_x1);
perturbed_x1 = noise_x1 + x1;
%snrdb_x1 = 20*log10(sum(x1.^2)/sum(noise_x1.^2));
snrdb_x1 = 20*log10(energy/(2*(sigma_x1^2)));

%Predict message
r1 = zeros(1, 10^5);
r2 = zeros(1, 10^5);
for i = 1:10^5
    r1(i) = sum(perturbed_x1(500*(i-1)+1:500*i).*psi1);
    r2(i) = sum(perturbed_x1(500*(i-1)+1:500*i).*psi2);      
end

m0 = r1 >= r2 & r1 >= -r2;
m1 = r1 < r2 & r1 < -r2;
m2 = r1 < r2 & r1 >= -r2;
m3 = r1 >= r2 & r1 < -r2;

locm0 = find(m0 == 1);
locm1 = find(m1 == 1);
locm2 = find(m2 == 1);
locm3 = find(m3 == 1);

predictions1 = string(zeros(1,10^5));

for i = 1:size(locm0, 2)
    predictions1(locm0(i)) = "00";
end

for i = 1:size(locm1, 2)
    predictions1(locm1(i)) = "01";
end

for i = 1:size(locm2, 2)
    predictions1(locm2(i)) = "10";
end

for i = 1:size(locm3, 2)
    predictions1(locm3(i)) = "11";
end

prob_error_x1 = 1 - sum(predictions1 == symbols1)/(10^5);

%Iteration 2
bits2 = randi([0,1],1,2*10^5);
bits2 = string(reshape(bits2, [], 2));
symbols2 = reshape(strcat(bits2(:,1), bits2(:,2)), [1, 100000]);

%Generate signal
x2 = zeros(1,500*10^5);
loc0 = find(symbols2 == "00");
loc1 = find(symbols2 == "01");
loc2 = find(symbols2 == "10");
loc3 = find(symbols2 == "11");

for i = 1:size(loc0, 2)
    x2(500*(loc0(i)-1)+1:500*loc0(i)) = s0;
end

for i = 1:size(loc1, 2)
    x2(500*(loc1(i)-1)+1:500*loc1(i)) = s1;
end

for i = 1:size(loc2, 2)
    x2(500*(loc2(i)-1)+1:500*loc2(i)) = s2;
end

for i = 1:size(loc3, 2)
    x2(500*(loc3(i)-1)+1:500*loc3(i)) = s3;
end

sigma_x2 = 5;
noise_x2 = randn(size(x1))*(sigma_x2);
perturbed_x2 = noise_x2 + x2;
snrdb_x2 = 20*log10(energy/(2*(sigma_x2^2)));

%Predict message
r1 = zeros(1, 10^5);
r2 = zeros(1, 10^5);
for i = 1:10^5
    r1(i) = sum(perturbed_x2(500*(i-1)+1:500*i).*psi1);
    r2(i) = sum(perturbed_x2(500*(i-1)+1:500*i).*psi2);  
end

m0 = r1 >= r2 & r1 >= -r2;
m1 = r1 < r2 & r1 < -r2;
m2 = r1 < r2 & r1 >= -r2;
m3 = r1 >= r2 & r1 < -r2;

locm0 = find(m0 == 1);
locm1 = find(m1 == 1);
locm2 = find(m2 == 1);
locm3 = find(m3 == 1);

predictions2 = string(zeros(1,10^5));

for i = 1:size(locm0, 2)
    predictions2(locm0(i)) = "00";
end

for i = 1:size(locm1, 2)
    predictions2(locm1(i)) = "01";
end

for i = 1:size(locm2, 2)
    predictions2(locm2(i)) = "10";
end

for i = 1:size(locm3, 2)
    predictions2(locm3(i)) = "11";
end

prob_error_x2 = 1 - sum(predictions2 == symbols2)/(10^5);

%Iteration 3
bits3 = randi([0,1],1,2*10^5);
bits3 = string(reshape(bits3, [], 2));
symbols3 = reshape(strcat(bits3(:,1), bits3(:,2)), [1, 100000]);

%Generate signal
x3 = zeros(1,500*10^5);
loc0 = find(symbols3 == "00");
loc1 = find(symbols3 == "01");
loc2 = find(symbols3 == "10");
loc3 = find(symbols3 == "11");

for i = 1:size(loc0, 2)
    x3(500*(loc0(i)-1)+1:500*loc0(i)) = s0;
end

for i = 1:size(loc1, 2)
    x3(500*(loc1(i)-1)+1:500*loc1(i)) = s1;
end

for i = 1:size(loc2, 2)
    x3(500*(loc2(i)-1)+1:500*loc2(i)) = s2;
end

for i = 1:size(loc3, 2)
    x3(500*(loc3(i)-1)+1:500*loc3(i)) = s3;
end

sigma_x3 = 3;
noise_x3 = randn(size(x3))*(sigma_x3);
perturbed_x3 = noise_x3 + x3;
snrdb_x3 = 20*log10(energy/(2*(sigma_x3^2)));

%Predict message
r1 = zeros(1, 10^5);
r2 = zeros(1, 10^5);
for i = 1:10^5
    r1(i) = sum(perturbed_x3(500*(i-1)+1:500*i).*psi1);
    r2(i) = sum(perturbed_x3(500*(i-1)+1:500*i).*psi2);  
end

m0 = r1 >= r2 & r1 >= -r2;
m1 = r1 < r2 & r1 < -r2;
m2 = r1 < r2 & r1 >= -r2;
m3 = r1 >= r2 & r1 < -r2;

locm0 = find(m0 == 1);
locm1 = find(m1 == 1);
locm2 = find(m2 == 1);
locm3 = find(m3 == 1);

predictions3 = string(zeros(1,10^5));

for i = 1:size(locm0, 2)
    predictions3(locm0(i)) = "00";
end

for i = 1:size(locm1, 2)
    predictions3(locm1(i)) = "01";
end

for i = 1:size(locm2, 2)
    predictions3(locm2(i)) = "10";
end

for i = 1:size(locm3, 2)
    predictions3(locm3(i)) = "11";
end

prob_error_x3 = 1 - sum(predictions3 == symbols3)/(10^5);

%Theoretical result
sigmaa = [0.5:-0.001:0.055];
snrs = energy./((F_s/2)*2*(sigmaa.^2)); % edited snrs = energy./(2*(sigmaa.^2));
qs = 0.5*erfc((snrs.^(0.5))./(sqrt(2)));
y = 1-((1-qs).^2);
figure;
semilogy(20*log10(snrs), y);
xlabel("E_s/2\sigma^2 (dB)");
ylabel("P_e");
title("Theoretical probability of symbol error vs. SNR");
hold on;
plot(snrdb_x1, prob_error_x1, "r*");
plot(snrdb_x2, prob_error_x2, "k*");
plot(snrdb_x3, prob_error_x3, "g*");
legend(["Theoretical","\sigma = 10","\sigma = 5","\sigma = 3"]);