clear;
clc;
close all;

%% 1. 讀取 EEG data
filename = 'EEG資料3.csv';
data = readtable(filename);

raw_EEG = data.EEG(:);

% 讀取 reference
ref.alpha = data.alpha(:);
ref.beta  = data.beta(:);
ref.delta = data.delta(:);
ref.theta = data.theta(:);

fs = 200;

N = length(raw_EEG);
time = (0:N-1)' / fs;

% 檢查資料正常
if height(data) ~= N || any(~isfinite(raw_EEG))
    error('Data length mismatch or invalid EEG values.');
end

if N < 2
    error('Not enough EEG samples.');
end

fprintf('Number of samples: %d\n', N);
fprintf('Sampling frequency: %.2f Hz\n', fs);
fprintf('Duration: %.3f seconds\n', (N-1)/fs);

%% 2. 定義閉眼、張眼區間
closed1 = time < 20.42;
openEye = time >= 20.42 & time < 41.11;
closed2 = time >= 41.11;

closedMask = closed1 | closed2;

fprintf('\nSegment durations:\n');
fprintf('Eyes closed: %.2f seconds total\n', sum(closedMask)/fs);
fprintf('Eyes open: %.2f seconds\n', sum(openEye)/fs);

%% 3. 移除DC offset
raw_EEG = raw_EEG - mean(raw_EEG);

%% 4. 定義頻帶
bands = struct( ...
    'name',  {'Alpha', 'Beta', 'Delta', 'Theta'}, ...
    'field', {'alpha', 'beta', 'delta', 'theta'}, ...
    'range', {[8 13], [13 30], [1 4], [4 8]});

filter_order = 4;
designed = struct();

%% 5. Butterworth bandpass filters
for k = 1:length(bands)

    field = bands(k).field;

    % 正規化截止頻率
    wn = bands(k).range / (fs/2);

    % Butterworth bandpass filter
    [b, a] = butter(filter_order, wn, 'bandpass');
    % 進行濾波
    designed.(field) = filtfilt(b, a, raw_EEG);

    fprintf('%s filter: %.1f-%.1f Hz\n', ...
        bands(k).name, ...
        bands(k).range(1), bands(k).range(2));
end

%% 6. 畫圖
for k = 1:length(bands)

    field = bands(k).field;

    % 取出濾波後的訊號
    sig = designed.(field);
    x_designed_closed = sig(closedMask);
    x_designed_open = sig(openEye);

    % 取出ref訊號
    sig = ref.(field);
    x_ref_closed = sig(closedMask);
    x_ref_open = sig(openEye);

    % 設定視窗
    figure('Color', 'w', ...
        'Name', bands(k).name, ...
        'Position', [100 100 1100 750]);

    sgtitle(sprintf('%s Band (%.0f-%.0f Hz)', ...
        bands(k).name, ...
        bands(k).range(1), bands(k).range(2)), ...
        'FontWeight', 'bold', ...
        'FontSize', 15);

    % 畫Designed, eyes closed的圖
    subplot(2,2,1);
    plotSpectrum(x_designed_closed, fs);
    title('Designed - Eyes Closed');

    % 畫Reference, eyes closed的圖
    subplot(2,2,2);
    plotSpectrum(x_ref_closed, fs);
    title('Reference - Eyes Closed');

    % 畫Designed, eyes open的圖
    subplot(2,2,3);
    plotSpectrum(x_designed_open, fs);
    title('Designed - Eyes Open');

    % 畫Reference, eyes open的圖
    subplot(2,2,4);
    plotSpectrum(x_ref_open, fs);
    title('Reference - Eyes Open');

end

%% 7.function plotSpectrum
function plotSpectrum(x, fs)

    x = x(:);
    x = x - mean(x);
    N = length(x);
    if N < 2 || max(abs(x)) == 0
        plot(NaN, NaN);
        xlim([-100 100]);
        ylim([0 1]);
        xticks(-100:25:100);
        yticks(0:0.2:1);
        xlabel('Frequency (Hz)');
        ylabel('Normalized Amplitude');
        grid on;
        box on;
        return;
    end

    % Hann window 降低頻譜洩漏
    w = hann(N, 'periodic');
    x_windowed = x .* w;

    % 傅立葉轉換並將輸出重新排列
    Y = fftshift(fft(x_windowed));
    amplitude = abs(Y);

    % 正規化
    peak = max(amplitude);
    if peak > 0
        amplitude = amplitude / peak;
    end

    % 建立FFT每一個資料的頻率
    f = (-floor(N/2):ceil(N/2)-1)' * fs/N;

    % 畫圖
    plot(f, amplitude, 'LineWidth', 1);
    xlim([-100 100]);
    ylim([0 1]);

    xticks(-100:25:100);
    yticks(0:0.2:1);

    xlabel('Frequency (Hz)');
    ylabel('Normalized Amplitude');

    grid on;
    box on;
end
