function write_csv(fname, header, data)
%WRITE_CSV  Minimal numeric CSV writer that behaves identically in MATLAB and Octave.
fid = fopen(fname, 'w');
fprintf(fid, '%s', header{1}); fprintf(fid, ',%s', header{2:end}); fprintf(fid, '\n');
fprintf(fid, [repmat('%.12g,', 1, size(data,2)-1) '%.12g\n'], data');
fclose(fid);
end
