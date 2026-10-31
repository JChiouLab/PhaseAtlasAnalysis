function d = sa_m3_testdir(E, name)
%SA_M3_TESTDIR  Fresh scratch folder for one milestone-3 test.
d = fullfile(E.scratch, sprintf('%s_%s', name, datestr(now, 'yyyymmdd_HHMMSSFFF'))); %#ok<DATST,TNOW1>
mkdir(d);
end
