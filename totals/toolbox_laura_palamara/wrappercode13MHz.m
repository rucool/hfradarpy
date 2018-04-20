% % % %% get detiding coefficients
% % % 
% % % starttime=datenum(2014,6,1);
% % % endtime=datenum(2014,7,1);
% % % % will have to switch period to the Palmer tidal components
% % % period=[12.42, 12.00, 12.66, 23.93, 25.82];
% % % inputdir='/Users/palamara/Desktop/antarctic_test/test_data/test13MHz/';
% % % coefficientfile=[inputdir 'detidingcoefficients/detidingcoefficients' datestr(starttime,'yyyy_mm_dd_HHMM') 'to' datestr(endtime,'yyyy_mm_dd_HHMM') '.nc'];
% % % 
% % % saveDetidedCoeffs(coefficientfile,starttime,endtime,period,'inputdir',inputdir)
% % % 


%% detide, filter, divergence, divergence trend (for each file)
inputdir = '/Users/codar/Desktop/Data/totals/pldp/25MHz/0.5km/oi/nc/';
inputfile=[inputdir 'OI_PLDP_2014_01_21_1600.totals.nc'];
earlierinputfile=[inputdir 'RU_13MHz_2014_06_10_0300.totals.nc'];

% can be done as files come in
[~]=detideCurrents(inputfile,coefficientfile,'addtide',true,'adddetided',true);

% use on earlier file (maybe 6 h previous?)
% timestep option will have to switch to .5
[~]=filterCurrents(earlierinputfile,'timestep',1,'addlowpass',true,'addhighpass',true,'currenttype','X_detided');

% get divergence from raw u and v (can be done as files come in on all but
% filtered)
[~]=divergence('netcdf','ncfile',inputfile,'adddivergence',true,'addvorticity',true,'currenttype','X');

% get 7-day previous raw divergence trend
[~]=trendData(inputfile,'div',.1,'addtrend',true,'overwrite',true);

% get 7-day previous raw vorticity trend
[~]=trendData(inputfile,'vor',.02,'addtrend',true,'overwrite',true);

