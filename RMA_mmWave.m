clear; clc

offsetData = 10; % azimuth sample offset
antenna_num = 2; % antenna ID
numADCBits = 16;
numRX = 4;
isReal = false;
valid_frame = 250;
resolution_ratio = 8;

%% mmWave Sensor Parameters
c = 299792458;                      % Free Space Wave Velocity  [m/sec]
F0 = 77e9;                          % Start Frequency           [Hz]
K = 70.295e+12;                     % Frequency Slope           [Hz/sec]
T_idle = 7e-6;                      % Idle Time                 [sec]
T_ADC = 4.66e-6;                    % ADC Start Time            [sec]
ADC_Samples = 256;                  % ADC Samples               [개수]
Samples_r = 5000e+3;                % Sample Rate               [Hz]
T_ramp = 56.9e-6;                   % Ramp End Time             [sec]
No_Chirp = 128;                     % No of Chirp Loops         [개수]
No_Frame = 300;                     % No of Frames              [개수]
T_frame = 40e-3;                    % Frame Periodicity         [sec]

%% Calculated Parameters
B_ideal = K * T_ramp;               % Ideal Bandwidth
T_valid = ADC_Samples / Samples_r;  % Real Chirp Time
B_valid = K * T_valid;              % Valid Bandwidth
PRI = T_idle + T_ramp;              % Chrip Repetition Period
T_compliance = PRI * No_Chirp;      % Compliance Chirp Time
R_ideal = T_compliance / T_frame;   % Duty Cycle
T_active = T_ramp * No_Chirp;       % Active Chirping Time
R_active = T_active / T_frame;      % Active-Ramp Duty Cycle
F_end = F0 + B_ideal;               % End Frequency
F0v = F0 + T_ADC * K;               % Valid Start Frquency
F_cv = F0v + B_valid / 2;           % Valid Carrier Frequency
F_endv = F0v + B_valid;             % Valid End Freqeuncy
dx = 1 / valid_frame;               % interelement space in synthetic array (lenth L) [m], rail 1m 초당 5cm

nFFTx = resolution_ratio * valid_frame;
nFFTy = resolution_ratio * ADC_Samples;

X = linspace(-0.5, 0.5, valid_frame).'; % Position of Array Elements
Kx = linspace(-pi/dx, pi/dx, nFFTx).';  % Azimuthal Spatial Wavenumber

freq = linspace(F0v, F_endv, ADC_Samples);  % Wideband Frequency
Kr = 2 * (2 * pi * freq / c);               % Wavenumber

path = "drone/adc_data_Raw_0.bin";
mndata = readDCA1000(path, ADC_Samples, numADCBits, numRX, isReal);

raw = reshape(mndata(antenna_num,:), [ADC_Samples, No_Chirp, No_Frame]);
raw = permute(raw, [3, 2, 1]);

raw = raw(offsetData:offsetData-1+valid_frame,:,:);
raw_r = real(raw);
raw_i = imag(raw);
mean_r = mean(raw_r, 2);
mean_i = mean(raw_i, 2);
raw = mean_r + mean_i * 1i;
raw = squeeze(raw);
figure; imagesc(abs(raw)); axis xy
% raw = raw(:,1:ADC_Samples/2);

%% Data Padding for Azimuth FFT
sarDataPadded = single(raw);
sarDataPadded = padarray(sarDataPadded,[floor((nFFTx-valid_frame)/2)],0,'pre');
sarDataPadded = padarray(sarDataPadded,[ceil((nFFTx-valid_frame)/2)],0,'post');

%% Matched Filter
Ky = single(sqrt(Kr.^2 - Kx.^2)); % x : azimuth axis, r, y : range axis
sarDataFFT2D = fftshift(fft(sarDataPadded, [], 1), 1); % Spatial Frequency Domain
Rs = 1; % Target Range
% phaseFactor = exp(-1i * Rs * Ky); % Compensation Function % 식이랑 다르다.
phaseFactor = exp(-1i * Rs * (Kr - Ky)); % 수정한 식
phaseFactor((Kx.^2) > (2*Kr).^2) = 0; % Consider Positive Wavenumber

sarDataFFT = Ky .* sarDataFFT2D; % Amplitude Correction
sarDataFFT = sarDataFFT .* phaseFactor; % Phase Correction

X2 = [-nFFTx / 2 + 1:nFFTx / 2] * dx;
Ky_int = Ky(nFFTx / 2, :); % Data Interpolation to Ky_int Cross-Range Wavenumber
dy = 2 * pi / (Ky_int(end) - Ky_int(1));
Y2 = 0 + dy:dy:ADC_Samples * dy;

% Visualization
SAR_db_Vis(fft(raw, [], 2), X2, Y2, 'Range FFT from raw');

SAR_db_Vis(fft2(sarDataFFT), X2, Y2, 'After Matched Filter');
SAR_db_Vis(fft(sarDataFFT, [], 2), X2, Y2, 'aa');
SAR_db_Vis(sarDataFFT, X2, Y2, 'aaa');

%% Stolt Interpolation
sarDataFFTinterp = my_spline(Ky, sarDataFFT, Ky_int); % 1D Data Interpolation by Spline Function
sarDataFFTinterp(find(isnan(sarDataFFTinterp))) = 1E-30; % set all Nan values to 0

figure; imagesc(abs(sarDataFFT)); colorbar
figure; imagesc(abs(sarDataFFTinterp)); title('Proposed Kx-Ky Map'); colorbar
SAR_db_Vis(fft2(sarDataFFTinterp), X2, Y2, 'After Original Stolt Interpolation');

%% Last Step
sarImg = fft2(sarDataFFTinterp, nFFTx, nFFTy);
SAR_db_img = SAR_db_Vis(sarImg, X2, Y2, 'Final SAR Image (dB)');