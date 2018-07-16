function data = ReadBinaryData(path_name)
data = [];
fid = fopen(path_name);
magic = int8('magic');
read_magic = fread(fid, 5, 'int8');
for ii = 1:5
    if magic(ii) ~= read_magic(ii)
        return;
    end;
end;
n_row = fread(fid, 1, 'int32');
n_column = fread(fid, 1, 'int32');
data = zeros(n_row, n_column);
for ii = 1:n_row
    data(ii, :) = double(fread(fid, n_column, 'float32'));
end;
fclose(fid);
end