%----------------------------------------------------------------------------------
% Compute quantities of interest from E3SM sensitivity runs.
%
%  This script processes raw data from a E3SM sensitivity run and
%  computes the following quantities:
%
%    SIE - Arctic sea ice extent (total area of sea ice concentration > 15 percent)
%    SIA - Arctic sea ice area 
%    SIV - Arctic sea ice volume
%    SST - Sea surface temperature (average over 60-90 N)
%    TS - Surface air temperature (average over 60-90 N)
%    QS - Surface specific humidity (average over 60-90 N)
%    PRS - Snow precipitation (average over 60-90 N)
%    CLD - Low cloud coverage
%    FLNS - Net longwave flux at surface
%    SH - Sea level pressure over the Siberian High (40-65N and 80-120E)
%    AL - Sea level pressure over the Aleutian Low (30-65N and 160-220E)
%    BH - Sea level pressure over the Beaufort Sea (72.5-80N and 180-225E)
%
%  Sensitivity runs are 100 years long and output from this script includes
%  time series data in the form of (year,month) for each of the variables
%  listed above. In cases where the run did not complete a full 100 years
%  this script will compute values for the completed years are fill the
%  yearly output file with -999 to indicate missing data. 
%  
%  The time series data are located in ascii files: (QOI_short_name)_yearly.txt
%
%  This script will also compute scalar values of each QOI averaged over the last
%  30 years of the run and output to the file: QOIs.txt.  In the case where
%  less than 30 years total are available, all years will be used to compute the 
%  averaged QOI.
%
%----------------------------------------------------------------------------------

   clear all

%----------------------------------------------------------------------------------
% ikalash INPUT: Set path to E3SM run directory where output files are located and case name root
%----------------------------------------------------------------------------------
  filepath='/nscratch/ikalash/acme_scratch/sandiatoss3/master.A_WCYCL1850.ne4_oQU240.sensitivity_102/run/'
  fileroot='master.A_WCYCL1850.ne4_oQU240.sensitivity_102.cam.h0.'


%----------------------------------------------------------------------------------
% Starting year for sensitivity run and total number of years
% If the run did not complete a full 100 years this script will use
% the total number of years available
%----------------------------------------------------------------------------------
  nyears = 100;
  ystart = 676;

%----------------------------------------------------------------------------------
% Open a sample atmosphere file and get dimensions and coordinates
%----------------------------------------------------------------------------------

  % get filename
  z='0';
  if (ystart < 10)
    z='000';
  elseif (ystart < 100)
    z='00';
  end
  filename = strcat(filepath,fileroot,z,num2str(ystart),'-01.nc');

  % open file
  ncid=netcdf.open(char(filename),'NC_NOWRITE');

  % get dimensions
  dimid = netcdf.inqDimID(ncid,'ncol');   
  [name,ncol] = netcdf.inqDim(ncid,dimid)
  dimid = netcdf.inqDimID(ncid,'ilev');   
  [name,nlev] = netcdf.inqDim(ncid,dimid)

  % get lat/lon and grid cell area
  varid = netcdf.inqVarID(ncid,'lat');   
  lat = netcdf.getVar(ncid,varid);
  varid = netcdf.inqVarID(ncid,'lon');   
  lon = netcdf.getVar(ncid,varid);
  varid = netcdf.inqVarID(ncid,'area');   
  area = netcdf.getVar(ncid,varid);

  %[lat lon]

  % close file
  netcdf.close(ncid);

%----------------------------------------------------------------------------------
% Get indices for regions to sum over
%----------------------------------------------------------------------------------
  
  % indN - corresponds to 60-90N latitude
  indN = find(lat >= 60);

  % indBH - corresponds to 72.5-80N latitude and 180-225 E longitude
  ind0 = find(lat >= 72.5 & lat <=80);
  ind1 = find(lon >= 180 & lon <=225);
  indBH = intersect(ind0,ind1);

  % indAL - corresponds to 30.0-65N latitude and 160-220 E longitude
  ind0 = find(lat >= 30 & lat <=65);
  ind1 = find(lon >= 160 & lon <=220);
  indAL = intersect(ind0,ind1);

  % indSH - corresponds to 40.0-65N latitude and 80-120 E longitude
  ind0 = find(lat >= 40 & lat <=65);
  ind1 = find(lon >= 80 & lon <=120);
  indSH = intersect(ind0,ind1);

  % get cell areas in those regions for cell average
  areaN = area(indN);
  areaBH = area(indBH);
  areaAL = area(indAL);
  areaSH = area(indSH);

%----------------------------------------------------------------------------------
%  Open sample sea ice file and Get cell area, lat, and lon
%----------------------------------------------------------------------------------

 % set filename and open file
   if (ystart < 10)
      filename = strcat(filepath,'mpascice.hist.000',num2str(ystart),'-01-01_00000.nc');
   elseif (ystart < 100)
      filename = strcat(filepath,'mpascice.hist.00',num2str(ystart),'-01-01_00000.nc');
   else
      filename = strcat(filepath,'mpascice.hist.0',num2str(ystart),'-01-01_00000.nc');
   end
   ncid=netcdf.open(filename,'NC_NOWRITE');

 % get dimensions
   dimid = netcdf.inqDimID(ncid,'nVertices');
   [dimname,num_verts] = netcdf.inqDim(ncid,dimid);
   dimid = netcdf.inqDimID(ncid,'nCells');
   [dimname,num_cells] = netcdf.inqDim(ncid,dimid);

 % get cell area, latitude and longitude
   varid = netcdf.inqVarID(ncid,'latCell');        % cell latitude
   lat_cell = netcdf.getVar(ncid,varid);
   varid = netcdf.inqVarID(ncid,'lonCell');        % cell longitude
   lon_cell = netcdf.getVar(ncid,varid);
   varid = netcdf.inqVarID(ncid,'areaCell');       % cell area
   area_cell = netcdf.getVar(ncid,varid);

 % close file
   netcdf.close(ncid);

  % Get northern hemisphere cells for computing ice coverage (use lat > 30)
  % Use lower latitute cutoff to make sure all ice is counted
  m2tokm2 = 1000*1000;
  indN_ice = find(lat_cell>pi/6);
  areaN_ice = area_cell(indN_ice)/m2tokm2;

  % Get northern hemisphere cells for computing sea surface temp (lat > 60)
  indN_sst = find(lat_cell>pi/3);
  areaN_sst = area_cell(indN_sst);

%----------------------------------------------------------------
% Read output files and compute the QOIs
%----------------------------------------------------------------

  % Arrays for time series by year and month
  sie = zeros(nyears,12); % sea ice extent
  sia = zeros(nyears,12); % sea ice area
  siv = zeros(nyears,12); % sea ice volume
  ave_sst  = zeros(nyears,12); % sea surface temperature
  ave_ts  = zeros(nyears,12); % surface temperature
  ave_qs  = zeros(nyears,12); % surface specific humidity
  ave_cld = zeros(nyears,12); % low cloud fraction
  ave_flns = zeros(nyears,12); % net longwave flux at surface
  ave_prs = zeros(nyears,12); % snow precipitation
  ave_bh  = zeros(nyears,12); % sea level pressure over Beaufort high
  ave_al  = zeros(nyears,12); % sea level pressure over Aleutian low
  ave_sh  = zeros(nyears,12); % sea level pressure over Siberian high
 
  % This is a string array that only works in version R2016b or highter
  %month = ["01","02","03","04","05","06","07",...
  %         "08","09","10","11","12"];

  % This is a char array that will work with older versions of Matlab
  month = {'01','02','03','04','05','06','07',...
           '08','09','10','11','12'};

  for i=1:nyears

     year = i + ystart -1

     for j=1:12

       % get file name for sea ice quantities
        ind = i + ystart-1;
        if (ind < 10)
           filename = strcat(filepath,'mpascice.hist.000',num2str(ind),'-',month(j),'-01_00000.nc');
        elseif (ind < 100)
           filename = strcat(filepath,'mpascice.hist.00',num2str(ind),'-',month(j),'-01_00000.nc');
        else
           filename = strcat(filepath,'mpascice.hist.0',num2str(ind),'-',month(j),'-01_00000.nc');
        end

       % if (exist(filename)==2)
        if (isfile(filename))

           iyear1=i;
           imonth1=j;

           ncid=netcdf.open(char(filename),'NC_NOWRITE');
           varid = netcdf.inqVarID(ncid,'iceAreaCell');
           conc = netcdf.getVar(ncid,varid);
           varid = netcdf.inqVarID(ncid,'iceVolumeCell');
           vol = netcdf.getVar(ncid,varid);
           netcdf.close(ncid);

           % for extent, sum areas of cells with ice concentration greater than 0.15
           ice_concN = conc(indN_ice);
           ind15 = find(ice_concN > 0.15);
           sie(i,j) = sum(areaN_ice(ind15));

           % to get ice area, sum ice concentration times area
           sia(i,j) = sum(areaN_ice.*ice_concN);

           % ice volume
           ice_volN = vol(indN_ice);
           siv(i,j) = sum(ice_volN.*areaN_ice./1000);

        else
           sie(i,j) = -999;
           sia(i,j) = -999;
           siv(i,j) = -999;
        end

       % get file name for sea surface temp
        ind = i + ystart - 1;
        if (ind < 10)
           filename = strcat(filepath,'mpascice.hist.am.timeSeriesStatsMonthly.000',num2str(ind),'-',month(j),'-01.nc');
        elseif (ind < 100)
           filename = strcat(filepath,'mpascice.hist.am.timeSeriesStatsMonthly.00',num2str(ind),'-',month(j),'-01.nc');
        else
           filename = strcat(filepath,'mpascice.hist.am.timeSeriesStatsMonthly.0',num2str(ind),'-',month(j),'-01.nc');
        end

        %if (exist(filename)==2)
        if (isfile(filename))

           iyear2=i;
           imonth2=j;

           ncid=netcdf.open(char(filename),'NC_NOWRITE');
           varid = netcdf.inqVarID(ncid,'timeMonthly_avg_seaSurfaceTemperature');
           sst_cell = netcdf.getVar(ncid,varid);
           netcdf.close(ncid);
 
           % Non-ocean cells are not included in output, so sum over all NH cells
           ave_sst(i,j) = sum(sst_cell(indN_sst).*areaN_sst)/sum(areaN_sst);
 
        else
           ave_sst(i,j) = -999;
        end

       % get file name for atmosphere quantities
        if (year < 10)
           filename = strcat(filepath, fileroot,'000',num2str(year),'-',month(j),'.nc');
        elseif (year < 100)
           filename = strcat(filepath, fileroot,'00',num2str(year),'-',month(j),'.nc');
        else
           filename = strcat(filepath, fileroot,'0',num2str(year),'-',month(j),'.nc');
        end

        %if (exist(filename)==2)
        if (isfile(filename))

           iyear3=i;
           imonth3=j;

           %open file
           ncid=netcdf.open(char(filename),'NC_NOWRITE');

           % get variables
           varid = netcdf.inqVarID(ncid,'TS'); % surface temperatuere
           ts = netcdf.getVar(ncid,varid); 
           varid = netcdf.inqVarID(ncid,'Q'); % specific humidity
           q = netcdf.getVar(ncid,varid); 
           varid = netcdf.inqVarID(ncid,'PSL'); % sea level pressure
           psl = netcdf.getVar(ncid,varid); 
           varid = netcdf.inqVarID(ncid,'CLDLOW'); % cloud low
           cld = netcdf.getVar(ncid,varid); 
           varid = netcdf.inqVarID(ncid,'PRECSL'); % snow precipitation
           prs = netcdf.getVar(ncid,varid); 
           varid = netcdf.inqVarID(ncid,'FLNS'); % net longwave flux at surface
           flns = netcdf.getVar(ncid,varid); 

           % close file
           netcdf.close(ncid);

           % get values in regions of interest
           tsN = ts(indN,1);
           q0 = q(:,nlev-1,1);
           qsN = q0(indN);
           cldN = cld(indN,1);
           prsN = prs(indN,1);
           flnsN = flns(indN,1);
           pslBH = psl(indBH,1);
           pslAL = psl(indAL,1);
           pslSH = psl(indSH,1);
       
           % get area averages of monthly data
           ave_ts(i,j) = sum(tsN.*areaN)/sum(areaN);
           ave_qs(i,j) = sum(qsN.*areaN)/sum(areaN);
           ave_cld(i,j) = sum(cldN.*areaN)/sum(areaN);
           ave_prs(i,j) = sum(prsN.*areaN)/sum(areaN);
           ave_flns(i,j) = sum(flnsN.*areaN)/sum(areaN);
           ave_bh(i,j) = sum(pslBH.*areaBH)/sum(areaBH);
           ave_al(i,j) = sum(pslAL.*areaAL)/sum(areaAL);
           ave_sh(i,j) = sum(pslSH.*areaSH)/sum(areaSH);

        else
           ave_ts(i,j) = -999;
           ave_qs(i,j) = -999;
           ave_cld(i,j) = -999;
           ave_prs(i,j) = -999;
           ave_flns(i,j) = -999;
           ave_bh(i,j) = -999;
           ave_al(i,j) = -999;
           ave_sh(i,j) = -999;
 
        end
     end
  end

  iyear1
  iyear2
  iyear3

  if (imonth1 < 12)
    iyear1 = iyear1-1;
  end
  if (imonth2 < 12)
    iyear2 = iyear2-1;
  end
  if (imonth3 < 12)
    iyear3 = iyear3-1;
  end


%---------------------------------------------------
% Write to text file
%---------------------------------------------------

  csvwrite('SIE_yearly.txt',sie(1:iyear1,:))
  csvwrite('SIA_yearly.txt',sia(1:iyear1,:))
  csvwrite('SIV_yearly.txt',siv(1:iyear1,:))
  csvwrite('SST_yearly.txt',ave_sst(1:iyear2,:))
  csvwrite('TS_yearly.txt',ave_ts(1:iyear3,:))
  csvwrite('QS_yearly.txt',ave_qs(1:iyear3,:))
  csvwrite('CLD_yearly.txt',ave_cld(1:iyear3,:))
  csvwrite('PRS_yearly.txt',ave_prs(1:iyear3,:))
  csvwrite('FLNS_yearly.txt',ave_flns(1:iyear3,:))
  csvwrite('BH_yearly.txt',ave_bh(1:iyear3,:))
  csvwrite('AL_yearly.txt',ave_al(1:iyear3,:))
  csvwrite('SH_yearly.txt',ave_sh(1:iyear3,:))

  % Create list of QOI values averaged over all years
  if (iyear1 > 30)
    siemean0 = mean(sie(iyear1-30:iyear1,:),2);
    siamean0 = mean(sia(iyear1-30:iyear1,:),2);
    sivmean0 = mean(siv(iyear1-30:iyear1,:),2);
  else
    siemean0 = mean(sie(1:iyear1,:),2);
    siamean0 = mean(sia(1:iyear1,:),2);
    sivmean0 = mean(siv(1:iyear1,:),2);
  end

  if (iyear2 > 30)
    sstmean0 = mean(ave_sst(iyear2-30:iyear2,:),2);
  else
    sstmean0 = mean(ave_sst(1:iyear2,:),2);
  end
 
  if (iyear3 > 30)
    tsmean0 = mean(ave_ts(iyear3-30:iyear3,:),2);
    qsmean0 = mean(ave_qs(iyear3-30:iyear3,:),2);
    cldmean0 = mean(ave_cld(iyear3-30:iyear3,:),2);
    prsmean0 = mean(ave_prs(iyear3-30:iyear3,:),2);
    flnsmean0 = mean(ave_flns(iyear3-30:iyear3,:),2);
    bhmean0 = mean(ave_bh(iyear3-30:iyear3,:),2);
    almean0 = mean(ave_al(iyear3-30:iyear3,:),2);
    shmean0 = mean(ave_sh(iyear3-30:iyear3,:),2);
  else
    tsmean0 = mean(ave_ts(1:iyear3,:),2);
    qsmean0 = mean(ave_qs(1:iyear3,:),2);
    cldmean0 = mean(ave_cld(1:iyear3,:),2);
    prsmean0 = mean(ave_prs(1:iyear3,:),2);
    flnsmean0 = mean(ave_flns(1:iyear3,:),2);
    bhmean0 = mean(ave_bh(1:iyear3,:),2);
    almean0 = mean(ave_al(1:iyear3,:),2);
    shmean0 = mean(ave_sh(1:iyear3,:),2);
  end

  qoi = zeros(12,1);
  
  format long;
  qoi(1) = mean(siemean0);
  qoi(2) = mean(siamean0);
  qoi(3) = mean(sivmean0);
  qoi(4) = mean(sstmean0);
  qoi(5) = mean(tsmean0);
  qoi(6) = mean(qsmean0);
  qoi(7) = mean(cldmean0);
  qoi(8) = mean(prsmean0);
  qoi(9) = mean(flnsmean0);
  qoi(10) = mean(bhmean0);
  qoi(11) = mean(almean0);
  qoi(12) = mean(shmean0);

  csvwrite('QOIs.txt',qoi);

