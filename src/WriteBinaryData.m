function WriteBinaryData(path_name, data)
[n_row, n_column] = size(data);
fid = fopen(path_name, 'w');
magic = int8('magic');
fwrite(fid, magic, 'int8');
fwrite(fid, int32(n_row), 'int32');
fwrite(fid, int32(n_column), 'int32');
for ii = 1:n_row
  fwrite(fid, single(data(ii, :)), 'float32'); 
end;
fclose(fid);
end