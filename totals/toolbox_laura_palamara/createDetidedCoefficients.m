%% get detiding coefficients

starttime=datenum(2014,5,1);
endtime=datenum(2014,6,1);

% will have to switch period to the Palmer tidal components
period=[12.42, 12.00, 12.66, 23.93, 25.82];
inputdir='/Users/codar/Desktop/Data/totals/maracoos/25MHz/1km/oi/nc/ideal/';
coefficientfile=['detidingcoefficients_' datestr(starttime,'yyyy_mm_dd_HHMM') 'to' datestr(endtime,'yyyy_mm_dd_HHMM') '.nc'];

if exist([inputdir 'detidingcoefficients/Archive'], 'dir')
    disp('Coefficient directories exist')
else
    mkdir([inputdir 'detidingcoefficients/Archive']);
end

saveDetidedCoeffs(coefficientfile,starttime,endtime,period,'inputdir',inputdir, 'overwrite', true)

