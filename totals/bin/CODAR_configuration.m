function conf = CODAR_configuration(systemType, sites, pattType, baseDir, saveDir, region)
%-------------------------------------------------
%  HF-Radar Processing Toolbox
%  Create Codar configuration file
%-------------------------------------------------

% Setup the configuration matrix directories, Sites, and Pattern Types

if systemType == 2
    domain = 'lr';
elseif systemType == 3
    domain = 'mr';
elseif systemType == 4
    domain = 'sr';
end

if num2str(region) == '3'
    regName = 'CARA';
elseif num2str(region) == '7'
    regName = 'MARA';
elseif num2str(region) == '13'
    regName = 'PLDP';
end

% Fill in Sites if the input was empty.. For off-line processing.
if isempty(sites)
    switch regName
        case 'MARA'
            switch domain
                case 'sr'
                    conf.Radials.Sites = {'PORT', 'SILD', 'MISQ', 'SUNS', 'SLTR', 'PORT', 'CBBT', 'GCAP', 'STLI', 'VIEW', 'BISL', 'CPHN'};
                    conf.Radials.Types = {'RDLm', 'RDLi', 'RDLi', 'RDLm', 'RDLi', 'RDLi', 'RDLm', 'RDLi', 'RDLi', 'RDLm', 'RDLi', 'RDLm'};
                case 'mr'
                    conf.Radials.Sites = {'BRAD', 'SPRK', 'BRNT','BRMR','RATH','WOOD'};
                    conf.Radials.Types = {'RDLm', 'RDLm', 'RDLm','RDLm','RDLm','RDLm'}; %remove belm august 31, 2012
                case 'lr'
                    conf.Radials.Types = {'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi', 'RDLi'}; 
                    conf.Radials.Sites = {'ASSA', 'BLCK', 'BRIG', 'CEDR', 'DUCK', 'HATY', 'HEMP', 'HOOK', 'LISL', 'LOVE', 'MRCH', 'MVCO', 'NANT', 'ERRA', 'NAUS', 'WILD'}; %Current
            end
        case 'CARA'
            conf.Radials.Sites = {'CDDO', 'FURA'};
            conf.Radials.Types = {'RDLi', 'RDLi'}; %remove belm august 31, 2012
        case 'PLDP'
            conf.Radials.Sites = {'JOUB', 'PALM', 'WAUW'};
            conf.Radials.Types = {'RDLi', 'RDLi', 'RDLi'};
    end
    suffix = '/';
else
    conf.Radials.Sites = sites;
    conf.Radials.Types = cell(size(conf.Radials.Sites));

    
    if isempty(pattType)
        emptyIndex = cellfun(@isempty, conf.Radials.Types);
        conf.Radials.Types(emptyIndex) = {'RDLi'}; %make ideal patterns default when pattType is left empty.
        suffix = '/ideal/';
    else
        if length(pattType) == 1
            emptyIndex = cellfun(@isempty, conf.Radials.Types);
            conf.Radials.Types(emptyIndex) = pattType;
            if strcmpi(pattType, 'RDLi')
                suffix = '/ideal/';
            else
                suffix = '/measured/';
            end
        else
           conf.Radials.Types = pattType; 
           suffix = '/';
        end
    end
end

% % % % % Radial Configuration
% % % % if isempty(pattType)
% % % %     emptyIndex = cellfun(@isempty, conf.Radials.Types);
% % % %     % Default processing - Ideal pattern
% % % %     conf.Radials.Types(emptyIndex) = 'RDLi';
% % % % else
% % % %     conf.Radials.Types = pattType;
% % % % end

% suffix = '/ideal/';
% suffix = '/measured/';


conf.Radials.BaseDir = [baseDir 'radials/'];
conf.Radials.RangeLims = [];
conf.Radials.BearLims = [];
conf.Radials.RangeBearSlop = repmat( 1e-10, [ numel(conf.Radials.Sites), 2 ] );
switch regName
    case 'MARA'
        switch domain
            case 'sr'
                conf.Radials.MaskDir = ['/Users/codar/Documents/MATLAB/totals_toolbox/mask_files/25MHz_1kmMask.txt'];% MASK FILE
                conf.Radials.MaxRadSpeed = 100;
                conf.Radials.MonthFlag = false;
                conf.Radials.TypeFlag=false;

                % Total Configuration
                conf.Totals.DomainName = 'MARASR';
                conf.Totals.CreationInfo = '25MHz/MARACOOS Domain';
                conf.Totals.BaseDir = [saveDir 'totals/maracoos/lsq/25MHz/'];
                conf.Totals.FilePrefix = strcat('tuv_',conf.Totals.DomainName,'_');
                conf.Totals.GridFile = ['/Users/codar/Documents/MATLAB/totals_toolbox/grid_files/old_grids/NY_1kmgrid.txt']; % GRID FILE
                conf.Totals.spatthresh = 1.6; %km 1.6
                conf.Totals.tempthresh = 1/24/2-eps; %
                conf.Totals.MaxTotSpeed = 150;
                conf.Totals.cleanTotalsVarargin = {{'GDOP','TotalErrors',1.05}};%1.25
                
                % OI Configuration
                conf.OI.BaseDir = [saveDir 'totals/maracoos/25MHz/1km/oi/mat' suffix];
                conf.OI.AsciiDir = [saveDir 'totals/maracoos/25MHz/1km/oi/ascii' suffix];
                conf.OI.FilePrefix = strcat('tuv_oi_',conf.Totals.DomainName,'_');
                conf.OI.FileSuffix = '.mat';
                conf.OI.mdlvar = 420;
                conf.OI.errvar = 66;
                conf.OI.sx = 1.6;% km 1.6
                conf.OI.sy = 1.6; % km 1.6
                conf.OI.tempthresh = 1/24/2-eps;
                conf.OI.cleanTotalsVarargin = {{'OIuncert','Uerr',0.7}, {'OIuncert','Verr',0.7}};%0.6
            case 'mr'
                conf.Radials.MaskDir = ['/home/codaradm/data/mask_files/BPU_2km_extendedMask.txt'];% MASK FILE
                conf.Radials.MaxRadSpeed = 150;
                conf.Radials.MonthFlag = false;

                % Total Configuration
                conf.Totals.DomainName = 'BPU';
                conf.Totals.CreationInfo = 'BPU/Domain';
                conf.Totals.BaseDir = [saveDir 'totals/maracoos/lsq/13MHz' suffix];
                conf.Totals.FilePrefix = strcat('tuv_',conf.Totals.DomainName,'_');
                conf.Totals.GridFile = ['/home/codaradm/data/grid_files/OI_2km_Grid.txt']; % GRID FILE expanded on September 7 2011
                conf.Totals.MaskFile = ['/home/codaradm/data/mask_files/BPU_2km_extendedMask.txt'];
                conf.Totals.spatthresh = 5; %km 1.6
                conf.Totals.tempthresh = 1/24/2-eps; %
                conf.Totals.MaxTotSpeed = 150;
                conf.Totals.cleanTotalsVarargin = {{'GDOP','TotalErrors',1.25}};
                conf.Radials.TypeFlag=false;

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
            case 'lr'
                conf.Radials.MaskDir = ['./mask_files/MARACOOS_6kmMask.txt'];
                conf.Radials.MaxRadSpeed = 300; %Please See notes
                conf.Radials.MonthFlag = true;
                conf.Radials.MonthSeperatorFlag = false;
                conf.Radials.TypeFlag=false;

                % Total Configuration
                conf.Totals.DomainName = 'MARA';
                conf.Totals.CreationInfo = 'Rutgers/MARACOOS Domain';
                conf.Totals.BaseDir = [saveDir 'totals/maracoos/lsq/5MHz/'];
                conf.Totals.FilePrefix = strcat('tuv_',conf.Totals.DomainName,'_');
                conf.Totals.GridFile = ['./grid_files/OI_6km_Grid.txt'];
                conf.Totals.MaskFile = ['./mask_files/MARACOOS_6kmMask.txt'];
                conf.Totals.spatthresh = 10; %km
                conf.Totals.tempthresh = 1/24/2-eps;
                conf.Totals.MaxTotSpeed = 300;% Please See Notes
                conf.Totals.cleanTotalsVarargin = {{'GDOP','TotalErrors',1.25}};

                % OI Configuration
                conf.OI.BaseDir = [saveDir 'totals/maracoos/oi/mat/5MHz' suffix];
                conf.OI.AsciiDir = [saveDir 'totals/maracoos/oi/ascii/5MHz' suffix];
                conf.OI.FilePrefix = strcat('tuv_oi_',conf.Totals.DomainName,'_');
                conf.OI.FileSuffix = '.mat';
                conf.OI.mdlvar = 420;
                conf.OI.errvar = 66;
                conf.OI.sx = 10;
                conf.OI.sy = 25;
                conf.OI.tempthresh = 1/24/2-eps;
                conf.OI.cleanTotalsVarargin = {{'OIuncert','Uerr',0.6}, {'OIuncert','Verr',0.6}};
        end
    case 'CARA'
        conf.Radials.MaskDir = ['/home/codaradm/data/mask_files/PR_coast.mask'];
        conf.Radials.MaxRadSpeed = 150;
        conf.Radials.MonthFlag = false;
        conf.Radials.MonthSeperatorFlag = false;
        conf.Radials.TypeFlag=false;

        % Total Configuration
        conf.Totals.DomainName = 'CARA';
        conf.Totals.CreationInfo = 'UPRM/CARACOOS Domain';
        conf.Totals.BaseDir = [saveDir 'totals/caracoos/lsq/5MHz/'];
        conf.Totals.FilePrefix = strcat('tuv_',conf.Totals.DomainName,'_');
        conf.Totals.GridFile = ['/home/codaradm/data/grid_files/UPRM_2km_grid.txt'];
        conf.Totals.MaskFile = ['/home/codaradm/data/mask_files/PR_coast.mask'];
        conf.Totals.spatthresh = 5; %km
        conf.Totals.tempthresh = 1/24/2-eps;
        conf.Totals.MaxTotSpeed = 150;
        conf.Totals.cleanTotalsVarargin = {{'GDOP','TotalErrors',1.25}};

        % OI Configuration
        conf.OI.BaseDir = [saveDir 'totals/caracoos/oi/mat/13MHz' suffix];
        conf.OI.AsciiDir = [saveDir 'totals/caracoos/oi/ascii/13MHz' suffix];
        conf.OI.FilePrefix = strcat('tuv_OI_',conf.Totals.DomainName,'_');
        conf.OI.FileSuffix = '.mat';
        conf.OI.mdlvar = 420;
        conf.OI.errvar = 66;
        conf.OI.sx = 5;
        conf.OI.sy = 8;
        conf.OI.tempthresh = 1/24/2-eps;
        conf.OI.cleanTotalsVarargin = {{'OIuncert','Uerr',0.6}, {'OIuncert','Verr',0.6}};
    case 'PLDP'
        conf.Radials.MaskDir = ['/Users/codar/Documents/MATLAB/totals_toolbox/mask_files/25MHz_1kmMask.txt'];% MASK FILE
        conf.Radials.MaxRadSpeed = 100;
        conf.Radials.MonthFlag = false;
        conf.Radials.TypeFlag=false;

        % Total Configuration
        conf.Totals.DomainName = 'PLDP';
        conf.Totals.CreationInfo = '25MHz/Palmer Deep Domain';
        conf.Totals.BaseDir = [saveDir 'totals/pldp/25MHz/lsq/'];
        conf.Totals.FilePrefix = strcat('tuv_',conf.Totals.DomainName,'_');
        conf.Totals.GridFile = {'/Users/codar/Documents/MATLAB/totals_toolbox/grid_files/PLDP_0.5km_grid.txt',...
            '/Users/codar/Documents/MATLAB/totals_toolbox/grid_files/PLDP_1km_grid.txt' };
        conf.Totals.spatthresh = 1.6; %km 1.6
        conf.Totals.tempthresh = 1/24/2-eps; %
        conf.Totals.MaxTotSpeed = 150;
        conf.Totals.cleanTotalsVarargin = {{'GDOP','TotalErrors',1.05}};%1.25

        % OI Configuration
        conf.OI.BaseDir = {[saveDir 'totals/pldp/25MHz/0.5km/oi/mat' suffix],...
            [saveDir 'totals/pldp/25MHz/1km/oi/mat' suffix]};
        conf.OI.AsciiDir = {[saveDir 'totals/pldp/25MHz/0.5km/oi/ascii' suffix],...
            [saveDir 'totals/pldp/25MHz/1km/oi/ascii' suffix]};
        conf.OI.FilePrefix = strcat('tuv_oi_',conf.Totals.DomainName,'_');
        conf.OI.FileSuffix = '.mat';
        conf.OI.mdlvar = 420;
        conf.OI.errvar = 66;
        conf.OI.sx = 1.6;% km 1.6
        conf.OI.sy = 1.6; % km 1.6
        conf.OI.tempthresh = 1/24/2-eps;
        conf.OI.cleanTotalsVarargin = {{'OIuncert','Uerr',0.7}, {'OIuncert','Verr',0.7}};%0.6
end
