function [fname] = codar_driver_totals(dtime,conf, systemType, varargin)
%
% Rutgers HFRadar Processing Toolbox
%
% Usage: CODAR_driver_totals(dtime, conf, varargin)
%
% This script processes radial files for input into total generation. 
%
% First, it converts the conf variable called from CODAR_configuration.m
% and converts it into a format that is useable by the HFRProgs toolbox.
% Then, it loads the radials from all available sites. After loading, it
% cleans the radials of bad data and masks the radials that are lying over
% the land. After this process, it calls makeTotalsOI.m to create the total
% vectors. Once the totals are generated, this script also cleans and masks
% out any bad data that may have been generated during total creation.
% 
% Results are written to the data directory in both matlab and ascii formats
%
% See also CheckTotals, CODAR_configuration, CODAR_driver_totals,
% makeTotalsOI
% 

%% Check to see if only one time was input
% if numel(dtime)>1
%   error('This function can only process one hour at a time');
% end

% Convert conf variable into array that HFRProfs can handle
conf = HFRPdriver_default_conf( conf );

% Mandatory parameters required by user
mand_params = { 'Radials.Sites', 'Radials.Types', 'Radials.RangeLims', ...
                'Radials.BearLims', 'Totals.DomainName', 'Totals.GridFile' };
            
% Check if mandatory input arguments are satisfied
conf = checkParamValInputArgs( conf, {}, mand_params, varargin{:} );

% Load Radial Data
fprintf(1, 'Loading Radial Data. \n');
fprintf(1, '---------------------------------------------------------------- \n');

% Construct filenames
F = filenames_standard_filesystem(conf.Radials.BaseDir, conf.Radials.Sites(:), conf.Radials.Types(:), dtime, conf.Radials.MonthFlag, conf.Radials.TypeFlag);

Rorig = loadRDLFile(F(:));
fprintf(1, '---------------------------------------------------------------- \n');

% If files contain no data, filter them out. 
ii = false(size(Rorig));
  for j = 1:numel(Rorig)
    ii(j) = numel(Rorig(j).U) == 0;
  end

  missingRadials.FileNames = [ Rorig(ii).FileName ];
  [missingRadials.TimeStamps,missingRadials.Sites,missingRadials.Types] ...
      = parseRDLFileName( missingRadials.FileNames );

  Rorig(ii) = [];
  if isempty(Rorig)
    error( 'No data at this timestep.' );
  end

% Get rid of stuff for missing files - if maskfile, rangelims or bearlims
% are missing, just don't do anything
try, conf.Radials.MaskFiles(ii) = []; end
try, conf.Radials.RangeLims(ii,:) = []; end
try, conf.Radials.BearLims(ii,:) = []; end

% Remove radials above MaxRadSpeed threshold
fprintf(1, 'Cleaning Radials. \n');
Rclean = cleanRadials( Rorig, conf.Radials.MaxRadSpeed );

% Remove radials over land
fprintf('Masking Radials. \n');
Rmask = maskRadials( Rclean, conf.Radials.MaskDir, 0);

%-------------------------------------------------
% Interpolate missing radials
fprintf(1, 'Interpolating Radials. \n');
fprintf(1, '---------------------------------------------------------------- \n');
for n = 1:numel(Rmask)
    % Allow for possibilty of RangeLims and/or BearLims to be not defined.
    try
        RL = conf.Radials.RangeLims(n,:);
    catch
%         fprintf('Radials.RangeLims(%d,:) not set, will default to [ ]\n',n);
        RL = [];
    end
    try
        BL = conf.Radials.BearLims(n,:);
    catch
        BL = [];
%         fprintf('Radials.BearLims(%d,:) not set, will default to [ ]\n',n);
    end
        
    % If there is only one radial, or maybe some other conditons that I'm 
    % not thinking of, then the interpolation will fail.  Set up a
    % try/catch block to keep things from failing.
    try
        Rinterp(n) = interpRadials( Rmask(n), ...
                            'RangeLims', RL, ...
                            'BearLims', BL, ...
                            'RangeDelta', conf.Radials.RangeBearSlop(n,1), ...
                            'BearDelta', conf.Radials.RangeBearSlop(n,2), ...
                            'MaxRangeGap', conf.Radials.RangeGap, ...
                            'MaxBearGap', conf.Radials.BearGap, ...
                            'CombineMethod', 'average' );
    catch
        fprintf(1, 'Warning: ## interpRadials failed for Site: %s\n',Rmask(n).SiteName);
        res=lasterror;
        fprintf(1, '%s\n',res.message)
        Rinterp(n) = Rmask(n);
        Rinterp(n).ProcessingSteps{end+1} = 'Interpolation failed, revert to uninterpolated';
    end
    
    % Check for case of interpolation creating all NaN's.  Use 90 % as the 
    % threshold.  Replace with uninterpolated radials and warn the user.
    if (sum(~isnan(Rinterp(n).U)) < sum(~isnan(Rmask(n).U)) * 0.9)
        fprintf(1, '%s:\n',char(Rinterp(n).FileName));
        fprintf(1, 'probably not interpolated properly ... using uninterpolated data instead\n')
        tmp = Rinterp(n).ProcessingSteps;
        Rinterp(n) = Rmask(n);
        Rinterp(n).ProcessingSteps = tmp;
        Rinterp(n).ProcessingSteps{end+1} = 'Revert to uninterpolated';
    end
end
fprintf(1, '---------------------------------------------------------------- \n');

% Load the total grid
[grid,fn,c] = loadDataFileWithChecks( conf.Totals.GridFile );
if c >= 100
  error( 'Could not find totals grid.' );
end

% Create OI Totals
% Start with Rmask (pre interpolation)
RTUV = Rmask;

% Call makeTotalsOI to generate the totals.
fprintf(1, 'Generating OI totals. \n');
[TUVorig, RTUV]=makeTotalsOI(RTUV,'Grid',grid,'TimeStamp',dtime, ...
      'mdlvar', conf.OI.mdlvar, 'errvar', conf.OI.errvar, ...
      'sx', conf.OI.sx, 'sy', conf.OI.sy, ...
      'tempthresh',conf.OI.tempthresh, ...
      'DomainName',conf.Totals.DomainName, ...
      'CreationInfo',conf.Totals.CreationInfo);

% Clean totals
[TUV,I] = cleanTotals( TUVorig, conf.Totals.MaxTotSpeed); %, ...
%conf.OI.cleanTotalsVarargin{:} );
%
fprintf(1, '%d totals removed by cleanTotals. \n',sum(I(:)>0));

if systemType == 2
    % Mask totals
     [TUV,I]=maskTotals(TUV,conf.Totals.MaskFile,false);
     fprintf(1, '%d totals masked out. \n',sum(~I(:)));
end
% 
% %% ------------------------------------------------------------------------
% %% This section of code was a contribution from Erick Fredj to gap fill the data
% % Masking is a destructive process, any totals current inside the mask will be
% % set to zero.
% 
% fprintf(1, 'Starting Smooth Total Field. \n');
% fprintf(1, '---------------------------------------------------------------- \n');
% 
% mask=load(conf.OSN.BestCoverageFile);
% hfrcvrg = inpolygon(TUV.LonLat(:,1),TUV.LonLat(:,2),mask(:,1),mask(:,2));
% 
% % Robust Smooth
% TUVosn = TUV;
% TUVosn.CreationInfo= 'Erick Fredj';
% 
% U=TUVosn.U(hfrcvrg);
% V=TUVosn.V(hfrcvrg);
% 
% % set to reset TUVs.U to NaN
% TUVosn.U = NaN(size(TUVosn.U));
% % set to reset TUVs.V to NaN
% TUVosn.V = NaN(size(TUVosn.V));
% 
% %% this function smoothn is located in toolbox_eric_fredj
% Vs = smoothn({U,V},'robust');
% 
% TUVosn.U(hfrcvrg)=Vs{1};
% TUVosn.V(hfrcvrg)=Vs{2};

%%-------------------------------------------------------------------------

% Save results
[tdn,tfn] = datenum_to_directory_filename( conf.OI.BaseDir, dtime, ...
                                           conf.OI.FilePrefix, ...
                                           conf.Totals.FileSuffix, conf.MonthFlag );
tdn = tdn{1};%

if ~exist( tdn, 'dir' )
  mkdir(tdn);
end
save(fullfile(tdn,tfn{1}),'conf','missingRadials','RTUV','TUVorig','TUV')