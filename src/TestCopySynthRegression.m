function output = TestCopySynthRegression(n_test)
output = false;
original_speech_dir = '~/Music/VCTK_CORPUS/VCTK-Corpus/wav48/';
target_analysis_dir = '~/m-file/STRAIGHTV40_007e/analysisData/';
target_wave_dir = '~/m-file/STRAIGHTV40_007e/waveData/';
target_files = dir([target_wave_dir '*.wav']);
n_files = length(target_files);
selected_id = randi(n_files, n_test);
for ii = 1:n_test
  tmp_name = target_files(selected_id(ii)).name;
  [x, fs] = audioread([original_speech_dir tmp_name(1:4) '/' tmp_name]);
  disp([num2str(ii) ': ' tmp_name ' ' datestr(now)]);
  rng(12345); % initialize frozen random number
  f0raw = MulticueF0v14(x, fs);
  ap = exstraightAPind(x, fs, f0raw);
  n3sgram=exstraightspec(x, f0raw, fs);
  tmp_name_root = tmp_name(1:end - 4);
  if ~CheckAnalysisData(f0raw, ap, n3sgram, target_analysis_dir, tmp_name_root)
    disp(['Failed: ' tmp_name ' data is not similar.']);
  end;
  wave_pathname = [target_wave_dir tmp_name];
  [sy, fs] = audioread(wave_pathname);
  y = exstraightsynth(f0raw,n3sgram,ap,fs);
  if std(sy - y / max(abs(y)) * 0.9) / std(sy) > 10 ^ (-3)
    disp(['Failed! ' tmp_name ' copy synthesis test.']);
    keyboard
    return;
  end;
end;
disp(['Success! ' num2str(n_test) ' files are passed copy-synth regression.']);
output = true;
end