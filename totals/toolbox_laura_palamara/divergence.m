function [lon,lat,u,v,div,vor]=divergence(inputtype,varargin)

% Get divergence and vorticity from u and v velocity on a rectangular grid
% at a single timestamp.
%
% Outputs:
%     div: instantaneous divergence
%     vor: instantaneous vorticity
%     lon: longitude used in calculation
%     lat: latitude used in calculation
%     u: eastward velocity used in calculation
%     v: northward velocity used in calculation
%
% Inputs:
%     inputtype: format of input variable
%         'manual': lon, lat, u, and v are put directly into the function
%             varargin options:
%             'lon': longitude 1-D array (required)
%             'lat': latitude 1-D array (required)
%             'u': eastward velocity 2-D array (required)
%             'v': northward velocity 2-D array (required)
%         'netcdf': lon, lat, u, and v are read from a netcdf file for a single timestamp
%             varargin options:
%             'ncfile': string containing netcdf file
%             'currenttype': variable name to read in and filter, formatted
%               with 'X' replacing 'u' or 'v' (i.e. to read in 'u_detided'
%               and 'v_detided', option would be 'X_detided'; default: 'X')
%             'adddivergence': logical value indicating whether to add
%               instantaneous divergence to the original netcdf file
%               (default: false)
%             'addvorticity': logical value indicating whether to add
%               instantaneous vorticity to the original netcdf file
%               (default: false)
%             'overwrite': logical value indicating whether to overwrite
%               divergence and vorticity if they already exist in the
%               original netcdf file (default: false)
%             'lonbounds': min and max longitude to read in (formatted
%               [minlon maxlon]) (default: read in all)
%             'latbounds': min and max latitude to read in (formatted
%               [minlat maxlat]) (default: read in all)
%         'tuv': lon, lat, u, and v are obtained from a TUV input
%             varargin options:
%             'tuv': TUV struct
%

lon=[];
lat=[];
u=[];
v=[];
div=[];
vor=[];

caller = [mfilename '.m'];

validinputtypes={'manual';'netcdf';'tuv'};

if nargin == 0
    fprintf(2,...
        '%s E: not enough arguments.\n',...
        caller);
    return;
elseif ~ischar(inputtype) | ~ismember(lower(inputtype),validinputtypes)
    fprintf(2,...
        '%s E: invalid input type.\n',...
        caller);
    return;
end

if(strcmpi(inputtype,'manual'))
    for x=1:2:length(varargin)
        name=varargin{x};
        value=varargin{x+1};
        switch lower(name)
            case 'lon'
                if(~isnumeric(value))
                    fprintf(2,'%s E: Value for %s must be a numeric array\n',caller,name);
                    return;
                end
                if(length(size(value))>2|~ismember(1,size(value)))
                    fprintf(2,'%s E: Value for %s must be a 1-D array\n',caller,name);
                    return;
                end
                lon=value;
            case 'lat'
                if(~isnumeric(value))
                    fprintf(2,'%s E: Value for %s must be a numeric array\n',caller,name);
                    return;
                end
                if(length(size(value))>2|~ismember(1,size(value)))
                    fprintf(2,'%s E: Value for %s must be a 1-D array\n',caller,name);
                    return;
                end
                lat=value;
            case 'u'
                if(~isnumeric(value))
                    fprintf(2,'%s E: Value for %s must be a numeric array\n',caller,name);
                    return;
                end
                if(length(size(value))~=2|~ismember(length(lon),size(value))|~ismember(length(lat),size(value)))
                    fprintf(2,'%s E: Dimensions for %s do not match lat and lon arrays\n',caller,name);
                    return;
                end
                u=value;
            case 'v'
                if(~isnumeric(value))
                    fprintf(2,'%s E: Value for %s must be a numeric array\n',caller,name);
                    return;
                end
                if(size(v)~=size(u))
                    fprintf(2,'%s E: Dimensions for %s do not match lat and lon arrays\n',caller,name);
                    return;
                end
                v=value;
                
            otherwise
                fprintf(2,...
                    '%s E: Unknown option specified for %s input: %s\n',...
                    caller,inputtype,...
                    name);
                return;
        end
    end
end


adddiv=false;
addvor=false;

if(strcmpi(inputtype,'netcdf'))
    inputvariable='X';
    replacevals=false;
    lonbounds=[-360 360];
    latbounds=[-90 90];
    
    for x=1:2:length(varargin)
        name=varargin{x};
        value=varargin{x+1};
        switch lower(name)
            case 'currenttype'
                if(~ischar(value))
                    fprintf(2,'%s E: Value for %s must be a string\n',caller,name);
                    return;
                end
                if(value(1)~='X')
                    fprintf(2,'%s E: Value for %s must begin with ''X'' (will be replaced with ''u'' and ''v'' in variable name)\n',caller,name);
                    return;
                end
                inputvariable=value;
            case 'lonbounds'
                if(~isnumeric(value))
                    fprintf(2,'%s E: Value for %s must be numeric\n',caller,name);
                    return;
                end
                if(length(value)~=2)
                    fprintf(2,'%s E: Value for %s must be 2 elements (min and max)\n',caller,name);
                    return;
                end
                lonbounds=value;
            case 'latbounds'
                if(~isnumeric(value))
                    fprintf(2,'%s E: Value for %s must be numeric\n',caller,name);
                    return;
                end
                if(length(value)~=2)
                    fprintf(2,'%s E: Value for %s must be 2 elements (min and max)\n',caller,name);
                    return;
                end
                latbounds=value;
        end
    end
    for x=1:2:length(varargin)
        name=varargin{x};
        value=varargin{x+1};
        switch lower(name)
            case 'currenttype'
            case 'lonbounds'
            case 'latbounds'
            case 'ncfile'
                if(~ischar(value))
                    fprintf(2,'%s E: Value for %s must be a string\n',caller,name);
                    return;
                end
                if(exist(value,'file')~=2|~strcmp(value(end-2:end),'.nc'))
                    fprintf(2,'%s E: Invalid netcdf file\n',caller);
                    return;
                end
                inputfile=value;
                
                varexist=false;
                divexist=false;
                vorexist=false;
                inputinfo=ncinfo(value);
                for k=1:length(inputinfo.Variables)
                    if(strcmp(inputinfo.Variables(k).Name,['u' inputvariable(2:end)]))
                        varexist=true;
                    elseif(strcmp(inputinfo.Variables(k).Name,['div' inputvariable(2:end)]))
                        divexist=true;
                    elseif(strcmp(inputinfo.Variables(k).Name,['vor' inputvariable(2:end)]))
                        vorexist=true;
                    end
                end
                
                if(~varexist)
                    fprintf(2,'%s E: Variables ''u %s'' and ''v %s'' not found in %s\n',caller,inputvariable(2:end),inputvariable(2:end),value);
                    return;
                end
                
                ncid=netcdf.open(value,'NC_NOWRITE');
                
                id=netcdf.inqVarID(ncid,'time');
                t=netcdf.getVar(ncid,id);
                if(isempty(t))
                    fprintf(2,'%s E: %s has no data\n',caller,inputfile);
                    netcdf.close(ncid)
                    return;
                end
                
                id=netcdf.inqVarID(ncid,'lon');
                lon=double(netcdf.getVar(ncid,id));
                if(length(size(lon))>2|~ismember(1,size(lon)))
                    fprintf(2,'%s E: Value for lon must be a 1-D array\n',caller);
                    return;
                end
                
                id=netcdf.inqVarID(ncid,'lat');
                lat=double(netcdf.getVar(ncid,id));
                if(length(size(lat))>2|~ismember(1,size(lat)))
                    fprintf(2,'%s E: Value for lat must be a 1-D array\n',caller);
                    return;
                end
                
                indlon=find(lon>=min(lonbounds)&lon<=max(lonbounds));
                if(isempty(indlon))
                    indlon=find(lon>=min(lonbounds+360)&lon<=max(lonbounds+360));
                    if(isempty(indlon))
                        indlon=find(lon>=min(lonbounds-360)&lon<=max(lonbounds-360));
                    end
                end
                indlat=find(lat>=min(latbounds)&lat<=max(latbounds));
                if(isempty(indlon)|isempty(indlat))
                    netcdf.close(ncid)
                    fprintf(2,'%s E: no gridpoints found between latitudes %f and %f and longitudes %f and %f\n',caller,min(latbounds),max(latbounds),min(lonbounds),max(lonbounds));
                    return;
                end
                
                sizelon=length(lon);
                sizelat=length(lat);
                
                lon=lon(indlon);
                lat=lat(indlat);
                
                
                id=netcdf.inqVarID(ncid,['u' inputvariable(2:end)]);
                u=squeeze(double(netcdf.getVar(ncid,id,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1])));
                u(u==netcdf.getAtt(ncid,id,'_FillValue'))=nan;
                udescription=netcdf.getAtt(ncid,id,'longname');
                if(length(size(u))~=2|~ismember(length(lon),size(u))|~ismember(length(lat),size(u)))
                    fprintf(2,'%s E: Dimensions for %s do not match lat and lon arrays\n',caller,name);
                    return;
                end
                
                id=netcdf.inqVarID(ncid,['v' inputvariable(2:end)]);
                v=squeeze(double(netcdf.getVar(ncid,id,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1])));
                v(v==netcdf.getAtt(ncid,id,'_FillValue'))=nan;
                vdescription=netcdf.getAtt(ncid,id,'longname');
                if(size(v)~=size(u))
                    fprintf(2,'%s E: Dimensions for v do not match lat and lon arrays\n',caller);
                    return;
                end
                
                netcdf.close(ncid)
                
            case 'adddivergence'
                if(~islogical(value))
                    fprintf(2,'%s E: Value for %s must be logical\n',caller,name);
                    return;
                end
                adddiv=value;
            case 'addvorticity'
                if(~islogical(value))
                    fprintf(2,'%s E: Value for %s must be logical\n',caller,name);
                    return;
                end
                addvor=value;
            case 'overwrite'
                if(~islogical(value))
                    fprintf(2,'%s E: Value for %s must be logical\n',caller,name);
                    return;
                end
                if(value&~adddiv&~addvor)
                    fprintf(1,'%s Warning: Either ''adddivergence'' or ''addvorticity'' must be set to true for %s to be used; no data will be added to file\n',caller,name);
                end
                replacevals=value;
                if(adddiv&~replacevals)
                    fprintf(1,'%s Warning: Divergence already exists in %s and ''replace'' is set to false; tidal component will not be added\n',caller,inputfile);
                end
                if(addvor&~replacevals)
                    fprintf(1,'%s Warning: Vorticity already exists in %s and ''replace'' is set to false; detided component will not be added\n',caller,inputfile);
                end
            otherwise
                fprintf(2,...
                    '%s E: Unknown option specified for %s input: %s\n',...
                    caller,inputtype,...
                    name);
                return;
                
        end
    end
end


if(strcmpi(inputtype,'tuv'))
    for x=1:2:length(varargin)
        name=varargin{x};
        value=varargin{x+1};
        switch lower(name)
            case 'tuv'
                if(~isstruct(value))
                    fprintf(2,'%s E: TUV input must be a structure\n',caller);
                    return;
                end
                lon=unique(value.LonLat(:,1));
                lat=unique(value.LonLat(:,2));
                
                if(abs((max(diff(lon))-min(diff(lon)))/min(diff(lon)))>1|abs((max(diff(lat))-min(diff(lat)))/min(diff(lat)))>1)
                    fprintf(2,'%s E: Grid must be rectangular\n',caller);
                    return;
                end
                
                u=nan(length(lon),length(lat));
                v=nan(length(lon),length(lat));
                
                for a=1:length(TUV.LonLat)
                    i=find(lon==TUV.LonLat(a,1));
                    j=find(lat==TUV.LonLat(a,2));
                    u(i,j)=TUV.U(a,1);
                    v(i,j)=TUV.V(a,1);
                end
                
            otherwise
                fprintf(2,...
                    '%s E: Unknown option specified for %s input: %s\n',...
                    caller,inputtype,...
                    name);
                return;
                
        end
    end
end



if(length(lon)==length(lat))
    fprintf(1,'%s Warning: Lat and Lon have same dimensions. Assuming Lon is 1st dimension and Lat is 2nd.\n',caller);
    londim=[1 0];
    latdim=[0 1];
else
    londim=double(size(u)==length(lon));
    latdim=double(size(u)==length(lat));
end

div=nan(size(u));
vor=nan(size(u));


% loop through points and get instantaneous divergence and vorticity
for i=2:size(u,1)-1
    for j=2:size(u,2)-1
        
        uy1=u(i,j+1); % u in y direction from center point
        uy2=u(i,j-1); % u in opposite y direction from center point
        vy1=v(i,j+1); % v in y direction from center point
        vy2=v(i,j-1); % v in opposite y direction from center point
        ux1=u(i+1,j); % u in x direction from center point
        ux2=u(i-1,j); % u in opposite x direction from center point
        vx1=v(i+1,j); % v in x direction from center point
        vx2=v(i-1,j); % v in opposite x direction from center point
        
        if(~any([uy1 uy2 vy1 vy2 ux1 ux2 vx1 vx2]==-999)) % make sure no surrounding points have a NaN
            
            dy1minus2=(lat(j+1)-lat(j-1))*111.12*100000; % calculate dy (from point -y1 to -y2)
            dx1minus2=(lon(i+1)-lon(i-1))*111.12*100000*cosd(lat(j)); % calculate dx (from point -x1 to -x2)
            
            % calculate dudx, dvdx, dudy, dvdy at one gridpoint distance
            % from center point
            dudx=(ux1-ux2)./dx1minus2;
            dvdx=(vx1-vx2)./dx1minus2;
            dudy=(uy1-uy2)./dy1minus2;
            dvdy=(vy1-vy2)./dy1minus2;
            
            f=4*pi*sind(lat(j))/86400;
            
            div(i,j)=(dudx+dvdy)*86400; % get surface divergence in m/d
            vor(i,j)=(dvdx-dudy)/f; % get surface vorticity normalized by f
            
        end
    end
end


if(adddiv|addvor)
    outnc=netcdf.open(inputfile,'NC_WRITE');
    outlondimid=netcdf.inqDimID(outnc,'lon');
    outlatdimid=netcdf.inqDimID(outnc,'lat');
    outtdimid=netcdf.inqDimID(outnc,'time');
    
    if(adddiv)
        if(~divexist)
            netcdf.reDef(outnc)
            % define divergence variable
            outdivid=netcdf.defVar(outnc,['div' inputvariable(2:end)],'float',[outlondimid outlatdimid outtdimid]);
            netcdf.putAtt(outnc,outdivid,'longname', 'Surface Divergence');
            netcdf.putAtt(outnc,outdivid,'units','m/day');
            netcdf.putAtt(outnc,outdivid,'description',['instantaneous surface divergence calculated using ' udescription ' and ' vdescription ' at neighboring gridpoints']);
            netcdf.putAtt(outnc,outdivid,'_FillValue',single(-999.0));
            netcdf.defVarFill(outnc,outdivid,false,single(-999));
            netcdf.endDef(outnc)
        end
        
        if(divexist&replacevals)
            outdivid=netcdf.inqVarID(outnc,['div' inputvariable(2:end)]);
            netcdf.reDef(outnc)
            netcdf.delAtt(outnc,outdivid,'description');
            netcdf.putAtt(outnc,outdivid,'description',['instantaneous surface divergence calculated using ' udescription ' and ' vdescription ' at neighboring gridpoints']);
            netcdf.endDef(outnc)
        end
        
        if(~divexist|replacevals)
            div(isnan(div))=-999;
            netcdf.putVar(outnc,outdivid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1],single(div));
            div(div==-999)=nan;
        end
    end
    
    if(addvor)
        if(~vorexist)
            netcdf.reDef(outnc)
            % define vorticity variable
            outvorid=netcdf.defVar(outnc,['vor' inputvariable(2:end)],'float',[outlondimid outlatdimid outtdimid]);
            netcdf.putAtt(outnc,outvorid,'longname', 'Surface Vorticity');
            netcdf.putAtt(outnc,outvorid,'units','unitless');
            netcdf.putAtt(outnc,outvorid,'description',['instantaneous surface vorticity (normalized by local f) calculated using ' udescription ' and ' vdescription ' at neighboring gridpoints']);
            netcdf.putAtt(outnc,outvorid,'_FillValue',single(-999.0));
            netcdf.defVarFill(outnc,outvorid,false,single(-999));
            netcdf.endDef(outnc)
        end
        
        if(vorexist&replacevals)
            outvorid=netcdf.inqVarID(outnc,['vor' inputvariable(2:end)]);
            netcdf.reDef(outnc)
            netcdf.delAtt(outnc,outvorid,'description');
            netcdf.putAtt(outnc,outvorid,'description',['instantaneous surface vorticity (normalized by local f) calculated using ' udescription ' and ' vdescription ' at neighboring gridpoints']);
            netcdf.endDef(outnc)
        end
        
        if(~vorexist|replacevals)
            vor(isnan(vor))=-999;
            netcdf.putVar(outnc,outvorid,[min(indlon)-1 min(indlat)-1 0],[length(indlon) length(indlat) 1],single(vor));
            vor(vor==-999)=nan;
        end
    end
    
    netcdf.close(outnc)
end

