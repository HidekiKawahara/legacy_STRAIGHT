%% Regression test data generator
% This program should be executed at the very beginning of refactoring a
% major revision. This is for making legacy STRAIGHT v40_007d to be
% compatible with MATLAB R2015b and Octave

% Copyright(c) 2016, Hideki Kawahara, (kawahara@sys.wakayama-u.ac.jp)

clear all
close all

original_speech_dir = '~/Music/VCTK_CORPUS/VCTK-Corpus/wav48/';
target_analysis_dir = '~/m-file/STRAIGHTV40_007e/analysisData/';
target_wave_dir = '~/m-file/STRAIGHTV40_007e/waveData/';

mkdir(target_analysis_dir);
mkdir(target_wave_dir);

dir_list = dir([original_speech_dir 'p*']);

%%

n_dirs = length(dir_list);
n_files = 0;
for ii = 1:n_dirs
  tmp_files = dir([original_speech_dir dir_list(ii).name '/*.wav']);
  n_files = n_files + length(tmp_files);
end;
%%
n_test = 2; % number of files tested for each speaker
l_segment = 0.1; % 100 ms segment
for ii = 1:n_dirs
  rng(12345); % initialize frozen random number
  seq_id = 0;
  basic_stat_table = zeros(n_files, 4);
  tmp_files = dir([original_speech_dir dir_list(ii).name '/*.wav']);
  for jj = 1:length(tmp_files)
    [x, fs] = audioread([original_speech_dir dir_list(ii).name '/' ...
      tmp_files(jj).name ]);
    seq_id = seq_id + 1;
    l_in_sample_segment = min(length(x), round(fs * l_segment));
    n_segment = floor(length(x) / l_in_sample_segment);
    rms_level = zeros(n_segment, 1);
    for kk = 1:n_segment
      rms_level(kk) = 20 * ...
        log10(std(x((kk - 1) * l_in_sample_segment + ...
        (1:l_in_sample_segment))));
    end;
    sorted_level = sort(rms_level);
    basic_stat_table(seq_id, 1) = length(x) / fs;
    basic_stat_table(seq_id, 2) = max(rms_level) - min(rms_level);
    basic_stat_table(seq_id, 3) = ...
      sorted_level(round(length(sorted_level) * 0.85));
    basic_stat_table(seq_id, 4) = max(abs(x));
  end;
  basic_stat_table = basic_stat_table(1:seq_id, :);
  % select safe region
  sorted_length = sort(basic_stat_table(:, 1));
  sorted_dynamic_range = sort(basic_stat_table(:, 2));
  sorted_85percent = sort(basic_stat_table(:, 3));
  l_10 = sorted_length(round(seq_id * 0.1));
  l_90 = sorted_length(round(seq_id * 0.9));
  d_10 = sorted_dynamic_range(round(seq_id * 0.1));
  d_90 = sorted_dynamic_range(round(seq_id * 0.9));
  v_10 = sorted_85percent(round(seq_id * 0.1));
  v_90 = sorted_85percent(round(seq_id * 0.9));
  index_list = 1:seq_id;
  safe_index = index_list( ...
    l_10 < basic_stat_table(:, 1) & l_90 > basic_stat_table(:, 1) & ...
    d_10 < basic_stat_table(:, 2) & d_90 > basic_stat_table(:, 2) & ...
    v_10 < basic_stat_table(:, 3) & v_90 > basic_stat_table(:, 3) & ...
    basic_stat_table(:, 4) < 0.95);
  selection_index = 1:length(safe_index);
  [~, tmp_index] = sort(rand(n_test, 1));
  selection_index = selection_index(tmp_index(1:n_test));
  for kk = 1:n_test
    id = safe_index(selection_index(kk));
    [x, fs] = audioread([original_speech_dir dir_list(ii).name '/' ...
      tmp_files(id).name ]);
    rng(12345); % initialize frozen random number
    f0raw = MulticueF0v14(x,fs);
    ap = exstraightAPind(x,fs,f0raw);
    n3sgram=exstraightspec(x,f0raw,fs);
    rng(12345); % initialize frozen random number
    y = exstraightsynth(f0raw,n3sgram,ap,fs);
    disp([num2str(kk) ': ' tmp_files(id).name  ' at:' datestr(now)]);
    audiowrite([target_wave_dir '/' tmp_files(id).name], ...
      y / max(abs(y)) * 0.9, fs);
    path_name_f0 = [target_analysis_dir '/' tmp_files(id).name(1:end-4) 'f0.bin'];
    path_name_ap = [target_analysis_dir '/' tmp_files(id).name(1:end-4) 'ap.bin'];
    path_name_sp = [target_analysis_dir '/' tmp_files(id).name(1:end-4) 'sp.bin'];
    WriteBinaryData(path_name_f0, f0raw)
    WriteBinaryData(path_name_ap, ap)
    WriteBinaryData(path_name_sp, n3sgram)
  end;
end;
