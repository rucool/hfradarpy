function conf = codar_configuration(system_type, baseDir, saveDir)
%-------------------------------------------------
%  HF-Radar Processing Toolbox
%  Create Codar configuration file
%-------------------------------------------------

% Setup the configuration matrix directories, Sites, and Pattern Types
conf.Radials.BaseDir = baseDir;
conf.Radials.RangeLims = [];
conf.Radials.BearLims = [];
conf.Radials.MonthFlag = false;
conf.Radials.MonthSeperatorFlag = false;
conf.Radials.TypeFlag=false;

switch system_type
    case 4
        conf.Radials.Sites = {'PORT', 'SILD', 'MISQ', 'SUNS', 'SLTR', 'PORT', 'CBBT', 'GCAP', 'STLI', 'VIEW', 'BISL', 'CPHN'};
        conf.Radials.Types = {'RDLm', 'RDLi', 'RDLi', 'RDLm', 'RDLi', 'RDLi', 'RDLm', 'RDLi', 'RDLi', 'RDLm', 'RDLi', 'RDLm'};
        conf.Radials.MaskDir = '/Users/codar/Documents/MATLAB/totals_toolbox/mask_files/25MHz_1kmMask.txt';% MASK FILE
        conf.Radials.MaxRadSpeed = 100;

        % Total Configuration
        conf.Totals.DomainName = 'MARASR';
        conf.Totals.CreationInfo = '25MHz/MARACOOS Domain';
        conf.Totals.BaseDir = [saveDir 'totals/maracoos/lsq/25MHz/'];
        conf.Totals.FilePrefix = strcat('tuv_',conf.Totals.DomainName,'_');
        conf.Totals.GridFile = '/Users/codar/Documents/MATLAB/totals_toolbox/grid_files/old_grids/NY_1kmgrid.txt'; % GRID FILE
        conf.Totals.spatthresh = 1.6; %km 1.6
        conf.Totals.tempthresh = 1/24/2-eps; %
        conf.Totals.MaxTotSpeed = 150;
        conf.Totals.cleanTotalsVarargin = {{'GDOP','TotalErrors',1.05}};%1.25

        % OI Configuration
        conf.OI.BaseDir = [saveDir 'totals/maracoos/25MHz/1km/oi/mat/'];
        conf.OI.AsciiDir = [saveDir 'totals/maracoos/25MHz/1km/oi/ascii/'];
        conf.OI.FilePrefix = strcat('tuv_oi_',conf.Totals.DomainName,'_');
        conf.OI.FileSuffix = '.mat';
        conf.OI.mdlvar = 420;
        conf.OI.errvar = 66;
        conf.OI.sx = 1.6;% km 1.6
        conf.OI.sy = 1.6; % km 1.6
        conf.OI.tempthresh = 1/24/2-eps;
        conf.OI.cleanTotalsVarargin = {{'OIuncert','Uerr',0.7}, {'OIuncert','Verr',0.7}};%0.6
    case 3
        conf.Radials.Sites = {'BRAD', 'SPRK', 'BRNT','BRMR','RATH','WOOD'};
        conf.Radials.Types = {'RDLm', 'RDLm', 'RDLm','RDLm','RDLm','RDLm'}; %remove belm august 31, 2012
        conf.Radials.MaskDir = ['/home/codaradm/data/mask_files/BPU_2km_extendedMask.txt'];% MASK FILE
        conf.Radials.MaxRadSpeed = 150;

        % Total Configuration
        conf.Totals.DomainName = 'BPU';
        conf.Totals.CreationInfo = 'BPU/Domain';
        conf.Totals.BaseDir = [saveDir 'totals/maracoos/lsq/13MHz' suffix];
        conf.Totals.FilePrefix = strcat('tuv_',conf.Totals.DomainName,'_');
        conf.Totals.GridFile = '/home/codaradm/data/grid_files/OI_2km_Grid.txt';
        conf.Totals.MaskFile = '/home/codaradm/data/mask_files/BPU_2km_extendedMask.txt';
        conf.Totals.spatthresh = 5; %km 1.6
        conf.Totals.tempthresh = 1/24/2-eps; %
        conf.Totals.MaxTotSpeed = 150;
        conf.Totals.cleanTotalsVarargin = {{'GDOP','TotalErrors',1.25}};

        % OI Configuration
        conf.OI.BaseDir = [saveDir 'totals/maracoos/oi/mat/13MHz' suffix];
        conf.OI.AsciiDir = [saveDir 'totals/maracoos/oi/ascii/13MHz' suffix];
        conf.OI.FilePrefix = strcat('tuv_oi_',conf.Totals.DomainName,'_');
        conf.OI.FileSuffix = '.mat';
        conf.OI.mdlvar = 420;
        conf.OI.errvar = 66;
        conf.OI.sx = 5;% km 1.6
        conf.OI.sy = 8; % km 1.6
        conf.OI.tempthresh = 1/24/2-eps;
        conf.OI.cleanTotalsVarargin = {{'OIuncert','Uerr',0.7}, {'OIuncert','Verr',0.7}};
    case 2
        conf.Radials.Types = {'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi'}; 
        conf.Radials.Sites = {'NANT', 'BLCK', 'AMAG', 'MRCH', 'HEMP', 'HOOK', 'LOVE', 'BRIG', 'WILD', 'ASSA', 'CEDR', 'LISL', 'DUCK', 'HATY', 'CORE'}; %Current
        conf.Radials.MaskDir = './mask_files/MARACOOS_6kmMask.txt';
        conf.Radials.MaxRadSpeed = 300; %Please See notes

        % Total Configuration
        conf.Totals.DomainName = 'MARA';
        conf.Totals.CreationInfo = 'Rutgers/MARACOOS Domain';
        conf.Totals.BaseDir = [saveDir 'totals/maracoos/lsq/5MHz/'];
        conf.Totals.FilePrefix = strcat('tuv_',conf.Totals.DomainName,'_');
        conf.Totals.GridFile = './grid_files/maracoos_grid_6km_extended.txt';
        conf.Totals.MaskFile = './mask_files/MARACOOS_6kmMask.txt';
        conf.Totals.spatthresh = 10; %km
        conf.Totals.tempthresh = 1/24/2-eps;
        conf.Totals.MaxTotSpeed = 300;% Please See Notes
        conf.Totals.cleanTotalsVarargin = {{'GDOP','TotalErrors',1.25}};

        % OI Configuration
        conf.OI.BaseDir = saveDir;
        conf.OI.AsciiDir = saveDir;
        conf.OI.FilePrefix = strcat('tuv_oi_',conf.Totals.DomainName,'_');
        conf.OI.FileSuffix = '.mat';
        conf.OI.mdlvar = 420;
        conf.OI.errvar = 66;
        conf.OI.sx = 10;
        conf.OI.sy = 10;
        conf.OI.tempthresh = 1/24/2-eps;
        conf.OI.cleanTotalsVarargin = {{'OIuncert','Uerr',0.6}, {'OIuncert','Verr',0.6}};
end
conf.Radials.RangeBearSlop = repmat( 1e-10, [ numel(conf.Radials.Sites), 2 ] );


