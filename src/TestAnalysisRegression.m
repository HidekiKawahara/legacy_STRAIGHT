function output = TestAnalysisRegression(n_test)
if ~isOctave
  rng('shuffle'); % initialize frozen random number
end;
output = false;
original_speech_dir = '/Users/kawahara/Music/VCTK_CORPUS/VCTK-Corpus/wav48/';
target_analysis_dir = '/Users/kawahara/m-file/STRAIGHTV40_007e/analysisData/';
target_wave_dir = '/Users/kawahara/m-file/STRAIGHTV40_007e/waveData/';
target_files = dir([target_wave_dir '*.wav']);
n_files = length(target_files);
selected_id = randi(n_files, n_test);
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
    return;
  end;
  %y = exstraightsynth(f0raw,n3sgram,ap,fs);
end;
disp(['Success! ' num2str(n_test) ' files are passed analysis regression.']);
output = true;
end