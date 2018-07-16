function output = TestCopySynthRegressionR(n_test)
output = false;
original_speech_dir = '/Users/kawahara/Music/VCTK_CORPUS/VCTK-Corpus/wav48/';
if isOctave
  target_analysis_dir = '/Users/kawahara/m-file/STRAIGHTV40_007e/analysisDataO/';
  target_wave_dir = '/Users/kawahara/m-file/STRAIGHTV40_007e/waveDataO/';
else
  target_analysis_dir = '/Users/kawahara/m-file/STRAIGHTV40_007e/analysisDataR/';
  target_wave_dir = '/Users/kawahara/m-file/STRAIGHTV40_007e/waveDataR/';
end;
target_files = dir([target_wave_dir '*.wav']);
n_files = length(target_files);
selected_id = randi(n_files, n_test);
command1 = 'rand("seed", 12345);';
command2 = 'randn("seed", 12345);';
for ii = 1:n_test
  tmp_name = target_files(selected_id(ii)).name;
  [x, fs] = audioread([original_speech_dir tmp_name(1:4) '/' tmp_name]);
  disp([num2str(ii) ': ' tmp_name ' ' datestr(now)]);
  if isOctave
    eval(command1);
    eval(command2);
  else
    rng(12345); % initialize frozen random number
  end;
  f0raw = MulticueF0v14(x, fs);
  ap = exstraightAPind(x, fs, f0raw);
  n3sgram=exstraightspec(x, f0raw, fs);
  tmp_name_root = tmp_name(1:end - 4);
  if ~CheckAnalysisData(f0raw, ap, n3sgram, target_analysis_dir, tmp_name_root)
    disp(['Failed: ' tmp_name ' data is not similar.']);
  end;
  wave_pathname = [target_wave_dir tmp_name];
  [sy, fs] = audioread(wave_pathname);
  if isOctave
    eval(command1);
    eval(command2);
  else
    rng(12345); % initialize frozen random number
  end;
  y = exstraightsynth(f0raw,n3sgram,ap,fs);
  disp(['Relative error SD: ' num2str(100 * std(sy - y / max(abs(y)) * 0.9) / std(sy)) ' %']);
  if std(sy - y / max(abs(y)) * 0.9) / std(sy) > 10 ^ (-3)
    disp(['Failed! ' tmp_name ' copy synthesis test.']);
    keyboard
    return;
  end;
end;
disp(['Success! ' num2str(n_test) ' files are passed copy-synth regression.']);
output = true;
end