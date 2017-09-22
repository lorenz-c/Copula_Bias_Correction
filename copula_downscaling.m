%function [] = copula_downscaling(fnme_x_cal, fnme_y_cal, fnme_x_val, varnme, maskfile, maskval, outdir, remssnl, node, ncores)
function [] = copula_downscaling(fnme_settings, nnodes, node, newrun)

% /app/matlab/2015b/bin/./matlab -nodisplay -nodesktop -r "cd('/home/lorenz-c/SMOS_Downscaling'); copula_downscaling('/home/lorenz-c/SMOS_Downscaling/SMOS_val_downscaling_settings.mat', 4, 1)"
% Load the settings from the submitted filename
load(fnme_settings)

fnme_x_cal       = settings.fnme_x_cal;
fnme_y_cal       = settings.fnme_y_cal;
fnme_x_val       = settings.fnme_x_val;
varnme           = settings.varnme;
maskfile         = settings.maskfile;
maskval          = settings.maskval;
fnme_out         = settings.fnme_out;
varnme_out       = settings.varnme_out;
remssnl          = settings.remssnl;
merge_settings   = settings.merge_settings;
warnings         = settings.warnings;
combine_calval   = 0;



% Set number of cores etc.
ncores                = str2num(getenv('SLURM_NTASKS_PER_NODE'));
user_clust            = parcluster('local');
user_clust.NumWorkers = ncores;
parpool(user_clust, ncores)

% Turn warnings on or off
warning('off', 'all');
pctRunOnAll warning off


%% Add the path with copula-tools
addpath('/home/lorenz-c/bin/Copula_BCS/');
addpath('/home/lorenz-c/bin/Matlab_GeoTS_Tools/TS_Tools/')
addpath('/home/lorenz-c/bin/Matlab_GeoTS_Tools/TS_Base/')
addpath('/home/lorenz-c/bin/Matlab_GeoTS_Tools/Read_Write_Tools/')

%% Load the calibration- and validation-periods
tme_val     = ncread(fnme_x_val, 'time');
tme_cal     = ncread(fnme_x_cal, 'time');
tme_units   = ncreadatt(fnme_x_val, 'time', 'units');

% Combine (if desired) the calibration and validation periods
if combine_calval == 1
    tme_out = [tme_val; tme_cal];
elseif combine_calval == 2
    tme_out = [tme_cal; tme_val];
else
    tme_out = tme_val;
end


%% Create the output files
fnme_std    = [fnme_out, '_std.nc'];
fnme_cpla   = [fnme_out, '_copula.nc'];
fnme_theta  = [fnme_out, '_theta.nc'];

backup_file = [fnme_out, '_Backup.nc'];     


%% Load a land/sea mask
mask = ncread(maskfile, 'mask');
lat  = ncread(maskfile, 'lat');
long = ncread(maskfile, 'lon');
mask = mask';

nlat = length(lat);
nlon = length(long);

% Get the land ids
ids = find(mask == maskval);

% Delete old files 
if node == 1 
    if newrun == true
        if exist(fnme_out, 'file')
            delete(fnme_out);
        end
    
        if exist(fnme_std, 'file')
            delete(fnme_std);
        end

        if exist(fnme_cpla, 'file')
            delete(fnme_cpla);
        end
    
        if exist(fnme_theta, 'file')
            delete(fnme_theta);
        end
        
        if exist(backup_file, 'file')
            delete(backup_file);
        end
    
        % Create new output files
        create_3d_netcdf(fnme_out, varnme_out, tme_out, lat, long, ...
                                                                tme_units);
        create_3d_netcdf(fnme_std, [varnme_out, '_std'], ...
                                            tme_out, lat, long, tme_units);
        create_2d_netcdf(fnme_cpla, 'Copula', lat, long, ...
                                           settings.merge_settings.copnms);
        create_2d_netcdf(fnme_theta, 'Theta', lat, long);
        
    elseif newrun == false      
        if exist(backup_file, 'file')
            delete(backup_file);
        end
        copyfile(fnme_cpla, backup_file);
    end     
end

if newrun == false
    processed_pixels = ncread(backup_file, 'Copula');
    % Get the pixels corresponding to the masked area
    processed_pixels = processed_pixels(ids);
    % Set the pixels where no Copula could be fitted to NaN
    processed_pixels(processed_pixels == 0) = NaN;
    % Remove the ids where ~isnan(processed_pixels)
    ids(~isnan(processed_pixels)) = [];
end

if nnodes > 1
    ids_per_node = nearest(length(ids)/nnodes);
    start_id     = (node - 1)*ids_per_node + 1;

    if node == nnodes
        end_id = length(ids);
    else 
        end_id = node*ids_per_node;
    end
else
    start_id = 1;
    end_id   = length(ids);
end

ids = ids(start_id:end_id);

%% Loop over all  pixels
%parfor i = start_id:end_id
parfor i = 1:length(ids)
    % Transform the ID to row- and column-indices
    [~,out]=system('vmstat -s -S M | grep "free memory"');
    mem=sscanf(out,'%f  free memory');
    fprintf('Ind2sub...')
    [rw, clm]  = ind2sub([nlat nlon], ids(i));
    fprintf('Ncread...')
    % Load data for cal/val
    x_cal      = squeeze(ncread(fnme_x_cal, varnme, [clm, rw, 1], ...
                                                               [1 1 Inf]));
    y_cal      = squeeze(ncread(fnme_y_cal, varnme, [clm, rw, 1], ...
                                                               [1 1 Inf]));
    x_val      = squeeze(ncread(fnme_x_val, varnme, [clm, rw, 1], ...
                                                               [1 1 Inf]));
    
    fprintf('Combine...')                                    
    % Combine the calibration and validation period
    if combine_calval == 1
        x_val      = [x_val; x_cal];
    elseif combine_calval == 2
        x_val      = [x_cal; x_val];
    end
    
%     % Search for valid data in order to keep the output files small
%     val_ids = find(~isnan(x_val));

    if ~all(isnan(x_val))
        try 
            fprintf('Copula_merge...')
            [x_val_out, x_val_std, theta, cpla] = ...
                         copula_merge(x_cal, y_cal, x_val, merge_settings);
        catch
            x_val_out  = 1e+20*ones(length(tme_out), 1);
            x_val_std  = 1e+20*ones(length(tme_out), 1);
            val_ids    = 1e+20*ones(length(tme_out), 1);
            theta      = 1e+20;
            cpla       = 0;
        end
    else
        x_val_out  = 1e+20*ones(length(tme_out), 1);
        x_val_std  = 1e+20*ones(length(tme_out), 1);
        val_ids    = 1e+20*ones(length(tme_out), 1);
        theta      = 1e+20;
        cpla       = 1e+20;  
    end
    fprintf(['Memory: ' num2str(mem), 'MB...']) 
    fprintf('Writing.!\n')
    fprintf(['Writing ID ', num2str(ids(i)), '...'])
    
    ncwrite(fnme_out, varnme, x_val_out, [1, rw, clm]);
    ncwrite(fnme_std, [varnme, '_std'], x_val_out, [1, rw, clm]);

    ncwrite(fnme_cpla, 'Copula', cpla, [rw, clm]);
    ncwrite(fnme_theta, 'Theta', theta, [rw, clm]);
    
    fprintf(['Done! \n'])

end

end


% function res = remssnl(tme, dta)
% 
% res = NaN(length(dta), 1);
% 
% for i = 1:12
%     mnth_indx = find(tme(:, 2) == i);
%     mnth_mn   = nanmean(dta(mnth_indx));
%     
%     res(mnth_indx) = dta(mnth_indx) - mnth_mn;
% end
% 
% end

function [] = create_3d_netcdf(fnme, varnme, times, lat, lon, time_units)   
    
    ntimes = length(times);
    nlat   = length(lat);
    nlon   = length(lon);
    
    ncid        = netcdf.create(fnme, 'NETCDF4');
    
    time_dim_id = netcdf.defDim(ncid, 'time', ntimes);
    lat_dim_id  = netcdf.defDim(ncid, 'lat', nlat);
    lon_dim_id  = netcdf.defDim(ncid, 'lon', nlon);
    
    time_id     = netcdf.defVar(ncid, 'time', 'double', time_dim_id);
    lat_id      = netcdf.defVar(ncid, 'lat', 'double', lat_dim_id);
    lon_id      = netcdf.defVar(ncid, 'lon', 'double', lon_dim_id);
    
    netcdf.putAtt(ncid, time_id, 'units', time_units);
    netcdf.putAtt(ncid, lat_id, 'units', 'degrees_north');
    netcdf.putAtt(ncid, lon_id, 'units', 'degrees_east');
    
    var_dims    = [time_dim_id, lat_dim_id, lon_dim_id];
    var_id      = netcdf.defVar(ncid, varnme, 'double', var_dims);
    
    netcdf.defVarFill(ncid, var_id, false, 1e+20);
    
%     for i = 0:4
%          netcdf.defVarChunking(ncid, i, 'CONTIGUOUS');
%     end
    
    netcdf.endDef(ncid)

    ncwrite(fnme, 'time', times);
    ncwrite(fnme, 'lat', lat);
    ncwrite(fnme, 'lon', lon);      
end

function [] = create_2d_netcdf(fnme, varnme, lat, lon, varargin) 

    nlat   = length(lat);
    nlon   = length(lon);
    
    ncid = netcdf.create(fnme, 'NETCDF4');
    
    lat_dim_id  = netcdf.defDim(ncid, 'lat', nlat);
    lon_dim_id  = netcdf.defDim(ncid, 'lon', nlon);
    
    lat_id  = netcdf.defVar(ncid, 'lat', 'double', lat_dim_id);
    lon_id  = netcdf.defVar(ncid, 'lon', 'double', lon_dim_id);
    
    netcdf.putAtt(ncid, lat_id, 'units', 'degrees_north');
    netcdf.putAtt(ncid, lon_id, 'units', 'degrees_east');
    
    var_dims    = [lat_dim_id, lon_dim_id];
    var_id      = netcdf.defVar(ncid, varnme, 'double', var_dims);
    netcdf.defVarFill(ncid, var_id, false, 1e+20);
    
    if strcmp(varnme, 'Copula')
        if size(varargin{1}, 1) > 1
            flag_meanings = char(strjoin(['no_fit'; varargin{1}]));
        else
            flag_meanings = char(strjoin(['no_fit', varargin{1}]));
        end
        netcdf.putAtt(ncid, var_id, 'flag_values', 0:length(varargin{1}));
        netcdf.putAtt(ncid, var_id, 'flag_meanings', flag_meanings);
     end
        
    
%     for i = 0:4
%          netcdf.defVarChunking(ncid, i, 'CONTIGUOUS');
%     end
    
    netcdf.endDef(ncid)
    
    ncwrite(fnme, 'lat', lat);
    ncwrite(fnme, 'lon', lon);      
end




