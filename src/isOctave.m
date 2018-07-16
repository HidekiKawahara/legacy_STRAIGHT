function output = isOctave
v = ver;
output = strcmp('Octave', v(1).Name);
end