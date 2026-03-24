% =========================================================
% Processamento de Sinais — Trabalho 1
% Análise e Síntese de Tons
% Afonso Pires nº2021032
% =========================================================


% ---------------------------------------------------------
% Carregar o ficheiro de áudio

% Define o caminho para a pasta de áudio
audio_dir = '/Users/afonsopiresmacmini/ESML_OFFLINE/PS_T1/audio/'
figures_dir = '/Users/afonsopiresmacmini/ESML_OFFLINE/PS_T1/figures/';

[x, fs] = audioread('/Users/afonsopiresmacmini/ESML_OFFLINE/PS_T1/audio/cmajorscale-vibrophone.wav');

signal_duration = length(x) / fs; % Duração do sinal em segundos

N = length(x); % Número de amostras

time = linspace(0, signal_duration, N); % Vetor de tempo correspondente a cada amostra

figure('Visible', 'off') % Criar figura sem a mostrar no ecrã
plot(time, x, 'b');
grid on;
xlabel('Tempo (s)', 'FontSize', 10);
ylabel('Amplitude', 'FontSize', 10);
title('Escala de Dó Maior - Vibrafone', 'FontSize', 12);
xlim([0, signal_duration]);
ylim([-1, 1]);

figure_name = [figures_dir 'Escala_Completa.png'];
print(figure_name, '-dpng', '-r150');

% ---------------------------------------------------------
% Separação das notas individuais

% Nomes das notas - usados para nomear ficheiros e gráficos:
notes = {'Do3', 'Re3', 'Mi3', 'Fa3', 'Sol3', 'La3', 'Si3', 'Do4'};

% Os limites temporais de cada nota foram determinados por
% análise visual da forma de onda exportada no Passo 1.
% Esta abordagem é adequada para um ficheiro de áudio conhecido
% e com silêncios bem definidos entre notas.

note_times = [0.05, 7.9;   % Do3
              7.9,  15.3;  % Re3
              15.3, 21.9;  % Mi3
              21.9, 28.7;  % Fa3
              28.7, 34.6;  % Sol3
              34.6, 40.7;  % La3
              40.7, 46.8;  % Si3
              46.8, 55.4]; % Do4

fade_length = round(0.1 * fs); % Comprimento do fade em amostras (100 ms)

fade_in = linspace(0, 1, fade_length)'; % Curva de fade-in
fade_out = linspace(1, 0, fade_length)'; % Curva de fade-out

x_notes = cell(1, length(notes)); % Célula para armazenar as notas individuais

for i = 1:length(notes)

    idx_start = round(note_times(i, 1) * fs) + 1; % Índice de início da nota
    idx_end = round(note_times(i, 2) * fs);       % Índice de fim da nota

    segment = x(idx_start:idx_end); % Extrair segmento da nota

    segment(1:fade_length) = segment(1:fade_length) .* fade_in; % Aplicar fade-in

    segment(end-fade_length+1:end) = segment(end-fade_length+1:end) .* fade_out; % Aplicar fade-out

    x_notes{i} = segment; % Armazenar a nota processada

    % Exportar cada nota como ficheiro .wav individual

    filename = [audio_dir notes{i} '.wav'];
    audiowrite(filename, segment, fs, 'BitsPerSample', 24);

end

% ---------------------------------------------------------
% Análise FFT de cada nota

n_parciais = 30; % Número de parciais a identificar por nota

peaks_mag = zeros(n_parciais, length(notes)); % Matriz para armazenar magnitudes dos picos
peaks_freq = zeros(n_parciais, length(notes)); % Matriz para armazenar frequências dos picos

for i = 1:length(notes)

    segment = x_notes{i}; % Carregar a nota processada
    N_fft = length(segment); % Número de pontos para a FFT (igual ao comprimento do segmento)

    % Janela de Blackman-Harris — reduz o spectral leakage
    window = blackman(N_fft); % Janela de Blackman
    segment_windowed = segment .* window; % Aplicar a janela ao segmento

    X = fft(segment_windowed, N_fft) / N_fft; % Calcular a FFT e normalizar pela quantidade de amostras

    magnitude = 2 * abs(X(1:N_fft/2+1)); % Magnitude da FFT (apenas a metade positiva)

    freq_axis = linspace(0, fs/2, N_fft/2+1); % Eixo de frequências correspondente à FFT

    % Deteção dos picos espectrais
    % MinPeakDistance — distância mínima entre picos em bins - Evita que o mesmo pico seja detetado múltiplas vezes

    min_dist = round(0.003 * fs);
    [pks, locs] = findpeaks(magnitude, 'MinPeakDistance', min_dist, 'NPeaks', n_parciais, 'SortStr', 'descend');

    % Guarda as magnitudes e frequências dos picos detetados
    n_found = length(pks);
    peaks_mag(1:n_found, i)  = pks;
    peaks_freq(1:n_found, i) = freq_axis(locs);

    % Gráfico do espectro de magnitude para cada nota
    figure('Visible', 'off');
    plot(freq_axis, magnitude, 'r');
    hold on;
    scatter(freq_axis(locs), pks, 'b', 'filled');
    xlabel('Frequência [Hz]', 'FontSize', 10);
    ylabel('Magnitude', 'FontSize', 10);
    title(['Espectro FFT — ' notes{i}], 'FontSize', 11);
    xlim([0 fs/2]);
    figure_name = [figures_dir 'FFT_' notes{i} '.png'];
    print(figure_name, '-dpng', '-r150');

end
