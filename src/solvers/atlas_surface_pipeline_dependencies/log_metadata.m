function log_metadata(id, filename, purpose, par, L, M_BG, M_P, rho_GEF, rho_PAK, notes)
    % LOG_METADATA Appends an entry to metadata_log.csv in the current folder
    log_file = 'metadata_log.csv';
    
    % Ensure directory exists
    if ~exist(log_file, 'file')
        fid = fopen(log_file, 'w');
        fprintf(fid, 'timestamp,computation_id,filename,purpose,a,b,c,d,e,f,Dm,Dc,L,M_BG,M_P,rho_GEF,rho_PAK,notes\n');
        fclose(fid);
    end
    
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    
    % Prepare data row
    row = sprintf('"%s","%s","%s","%s",%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f,%.2f,%.2f,%.2f,%.2f,"%s"\n', ...
        timestamp, id, filename, purpose, ...
        par.a, par.b, par.c, par.d, par.e, par.f, ...
        par.Dm, par.Dc, L, M_BG, M_P, rho_GEF, rho_PAK, notes);
    
    % Append to file
    fid = fopen(log_file, 'a');
    fprintf(fid, '%s', row);
    fclose(fid);
end
