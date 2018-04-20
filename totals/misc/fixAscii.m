% Setup paths
javaaddpath('mysql-connector-java-5.1.6-bin.jar');
addpath(genpath('/home/codaradm/operational_scripts/totals_toolbox/'))
addpath(genpath('/home/codaradm/HFR_Progs-2_1_3beta/matlab/general'))
add_subdirectories_to_path('/home/codaradm/HFR_Progs-2_1_3beta/matlab/',{'CVS','private','@'});

% Add the snctools jar file
disp( 'Adding NetCDF libraries...' );
lasterr('');
try
    javaaddpath('/home/coolgroup/matlabToolboxes/netcdfAll-4.2.jar')
    addpath(genpath('/home/coolgroup/matlabToolboxes'));
catch
    disp(['NetCDF-4 Error: ' lasterr]);
end


tS = datenum(2009, 1, 1, 0, 0, 0);
tE = datenum(2009, 5, 20, 1, 0, 0);

tI = tE:-1/24:tS;
matDir = '/home/codaradm/data/totals/maracoos/oi/mat/5MHz/';
dataDir = '/home/codaradm/data/';
            
for m = 1:length(tI)
    tM = datevec(tI(m));
    tY = num2str(tM(1));
    
    if tM(2) < 10;
        tMo = ['0' num2str(tM(2))];
    else
        tMo = num2str(tM(2));
    end

    if tM(3) < 10;
        tD = ['0' num2str(tM(3))];
    else
        tD = num2str(tM(3));
    end
    
    if tM(4) == 0;
        tH = [num2str(tM(4)) '0'];
    elseif tM(4) < 10 & tM(4) > 0
        tH = ['0' num2str(tM(4))];
    else
        tH = [num2str(tM(4))];
    end
    
    fileName = [matDir tY '_' tMo '/tuv_oi_MARA_' tY '_' tMo '_' tD '_' tH '00.mat'];
    
    load(fileName);

    conf = CODAR_configuration(2, [], [], dataDir, dataDir, 7);

    % Clean totals
    [TUVfix,I] = cleanTotals( TUVorig, conf.Totals.MaxTotSpeed);%, ...
    fprintf(1, '%d totals removed by cleanTotals. \n',sum(I(:)>0));

    % Mask totals
    [TUVfix,I]=maskTotals(TUVfix,conf.Totals.MaskFile,false);
    fprintf(1, '%d totals masked out. \n',sum(~I(:)));

    % Check if Ascii Directory exists. If not, create it.
    if ~exist( conf.OI.AsciiDir, 'dir' )
      mkdir(conf.OI.AsciiDir);
    end

    % Save the totals to ASCII files
    procFname = TUVstruct2ascii_OI(TUVfix,conf.OI.AsciiDir);

    ingestCodar(2, procFname, 7);

end 
