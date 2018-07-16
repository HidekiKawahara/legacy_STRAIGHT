function output = ...
  CheckAnalysisData(f0raw, ap, n3sgram, target_analysis_dir, tmp_name_root)
output = true;
tolerance = 10 ^ (-6);
f0_file_path = ...
  [target_analysis_dir tmp_name_root 'f0.bin'];
ap_file_path = ...
  [target_analysis_dir tmp_name_root 'ap.bin'];
sp_file_path = ...
  [target_analysis_dir tmp_name_root 'sp.bin'];
f0_ref = ReadBinaryData(f0_file_path);
ap_ref = ReadBinaryData(ap_file_path);
sp_ref = ReadBinaryData(sp_file_path);
f0_median = median(f0_ref(f0_ref > 30 & f0_ref < 1000));
ap_std = std(ap_ref(:));
sp_std = std(sp_ref(:));
if std(f0raw(:) - f0_ref(:)) / f0_median > tolerance
  return;
end;
if std(ap(:) - ap_ref(:)) / ap_std > tolerance
  return;
end;
if std(n3sgram(:) - sp_ref(:)) / sp_std > tolerance
  return;
end;
end