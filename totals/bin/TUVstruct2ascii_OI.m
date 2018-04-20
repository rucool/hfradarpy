function [fname] = TUVstruct2ascii_OI(T,outDir)

try, outDir;
catch
    outDir = pwd;
end


uvFlag = -999.0;

% See how many valid times exist
index = find(isfinite(T.TimeStamp));
if isempty(index)
    return;
end

for i = 1:length(index)
    
    [yy,mm,dd,hh,MM,ss]=datevec(T.TimeStamp(index(i)));
    % Create the file name
    fname = sprintf('%s_%s_%.4d_%.2d_%.2d_%.2d%.2d', ...
            T.Type,T.DomainName,yy,mm,dd,hh,MM);
   
%         [tdn,tfn] = datenum_to_directory_filename( outDir, T.TimeStamp((i)),'OI_BPU_','' );%, conf.MonthFlag );%T.Type,T.DomainName

    fname = fullfile(outDir,fname);
    fout = fopen(fname,'w');
    
    % Create the header information
    fprintf(fout,'%%TimeStamp: %.4d %.2d %.2d %.2d %.2d\n', yy,mm,dd,hh,MM);
    fprintf(fout,'%%TimeZone: %s\n',T.TimeZone);
    fprintf(fout,'%%Domain: %s\n',T.DomainName);
    fprintf(fout,'%%Type: %s\n',T.Type);
    fprintf(fout,'%%DataCreationInfo: %s\n',T.CreationInfo);
    fprintf(fout,'%%DataCreationTimeStamp: %s\n', ...
                  datestr(T.CreateTimeStamp,0));
    fprintf(fout,'%%DataCreationTimeZone: %s\n',T.CreateTimeZone);
    fprintf(fout,'%%ProcessingProgram: %s %s\n',mfilename,datestr(now));
    fprintf(fout,'%%TUV_structVersion: %s\n',T.TUV_struct_version);

    fprintf(fout,'%%MinNumSites: %4.0f\n',T.OtherMetadata.makeTotalsOI.parameters.MinNumSites);
    fprintf(fout,'%%MinNumRads: %3.0f\n',T.OtherMetadata.makeTotalsOI.parameters.MinNumRads);
    fprintf(fout,'%%mdlvar: %6.2f\n',T.OtherMetadata.makeTotalsOI.parameters.mdlvar);
    fprintf(fout,'%%errvar: %6.2f\n',T.OtherMetadata.makeTotalsOI.parameters.errvar);
    fprintf(fout,'%%sx: %6.2f\n',T.OtherMetadata.makeTotalsOI.parameters.sx);
    fprintf(fout,'%%sy: %6.2f\n',T.OtherMetadata.makeTotalsOI.parameters.sy);
    fprintf(fout,'%%tempthresh: %7.6f\n',T.OtherMetadata.makeTotalsOI.parameters.tempthresh);

    fprintf(fout,'%%Longitude  Latitude  U comp  V comp    Uerr   Verr   NumRad   Site\n');
    fprintf(fout,'%% (deg)      (deg)    (cm/s)  (cm/s)  (norm) (norm)            Code\n');
    
    % Create the data - flag any NaN's with uvFlag for ascii output
    U = T.U(:,index(i));
    V = T.V(:,index(i));
    
    ind = isnan(U);
    U(ind) = uvFlag;
    V(ind) = uvFlag;
    ind = isnan(V);
    U(ind) = uvFlag;
    V(ind) = uvFlag;    

    Ue = T.ErrorEstimates.Uerr(:,index(i));
    Ve = T.ErrorEstimates.Verr(:,index(i));

    ind = isnan(Ue);
    Ue(ind) = uvFlag;
    Ve(ind) = uvFlag;
    ind = isnan(Ve);
    Ue(ind) = uvFlag;
    Ve(ind) = uvFlag;    

    nRad = T.OtherMatrixVars.makeTotalsOI_TotalsNumRads(:,index(i));
    SC   = T.OtherMatrixVars.makeTotalsOI_TotalsSiteCode(:,index(i));
    
    out = [T.LonLat U V Ue Ve nRad SC];
    fprintf(fout,'%10.5f %9.5f %7.2f %7.2f %7.2f %7.2f %6.0f %7.0f\n',out');
    fclose(fout);
end
