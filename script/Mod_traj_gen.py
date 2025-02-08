from __future__ import division
import os
import shutil
from subprocess import call
import itertools
import fnmatch
from calendar import monthrange
from Mod_traj_util import (_populate_control, _try_to_remove,
                                   _day2filenum, _cliptraj)
import datetime as dt

def generate_bulktraj(basename, hysplit_working, output_dir, meteo_dir, coordinates, run, hours, altitudes, startdate, enddate, meteoyr_2digits=True, outputyr_2digits=False, get_reverse=False, get_clipped=False, hysplit="../hyts_std"):
    """
    Generate a sequence of trajectories within the period defined by startdate and enddate.
    
    Parameters
    ----------
    basename : string
        Base name for all output files.
    hysplit_working : string
        Path to the HYSPLIT working directory.
    output_dir : string
        Output directory where trajectories are saved.
    meteo_dir : string
        Directory with meteorological files.
    coordinates : tuple of floats
        (latitude, longitude) for parcel launch.
    run : int
        Simulation duration (negative for back trajectories).
    hours : list
        List of hours (as strings or integers) at which to launch trajectories.
    altitudes : list
        List of altitudes (in meters) from which to launch parcels.
    startdate : string
        Start date in YYYYMMDD format.
    enddate : string
        End date in YYYYMMDD format.
    meteoyr_2digits : bool, optional
        Use 2-digit year for meteorology file naming.
    outputyr_2digits : bool, optional
        Use 2-digit year for trajectory file naming.
    get_reverse : bool, optional
        Whether to generate reverse trajectories.
    get_clipped : bool, optional
        Whether to generate clipped trajectories.
    hysplit : string, optional
        Path to the HYSPLIT executable.
    
    """
    # Convert start/end dates to datetime objects
    #start_dt = dt.datetime.strptime(startdate, "%Y%m%d")
    #end_dt = dt.datetime.strptime(enddate, "%Y%m%d")
    start_dt = dt.datetime.strptime(str(startdate), "%Y%m%d")
    end_dt = dt.datetime.strptime(str(enddate), "%Y%m%d")
    # Get directory information and create output subdirectories
    cwd = os.getcwd()
    traj_output_dir = os.path.join(output_dir, basename, "traj")
    os.makedirs(traj_output_dir, exist_ok=True)
    output_rdir = os.path.join(traj_output_dir, 'reversetraj')
    output_cdir = os.path.join(traj_output_dir, 'clippedtraj')
    meteo_dir = meteo_dir.replace('\\', '/')
    if get_reverse and not os.path.isdir(output_rdir):
        os.mkdir(output_rdir)
    if get_clipped and not os.path.isdir(output_cdir):
        os.mkdir(output_cdir)
    
    # Determine hemisphere and month dictionary for meteorology file selection.
    n_hemisphere = True
    if coordinates[0] < 0:
        n_hemisphere = False
    mon_dict = _mondict(n_hem=n_hemisphere)
    
    # Year formatting functions.
    yr_is2digits = {True: _year2string, False: str}
    controlyearfunc = yr_is2digits[True]       # always 2-digit for CONTROL file
    fnameyearfunc = yr_is2digits[outputyr_2digits]
    meteoyearfunc = yr_is2digits[meteoyr_2digits]
    
    controlfname = 'CONTROL'
    
    try:
        os.chdir(hysplit_working)
        current_dt = start_dt
        # Iterate through all dates from start_dt to end_dt (inclusive)
        while current_dt <= end_dt:
            # Loop over the specified hours
            for h in hours:
                # Create simulation start datetime for the given day and hour.
                sim_dt = current_dt.replace(hour=int(h))
                if sim_dt < start_dt or sim_dt > end_dt:
                    continue
                for a in altitudes:
                    y = sim_dt.year
                    m = sim_dt.month
                    d = sim_dt.day
                    hr = sim_dt.hour
                    
                    controlyr = controlyearfunc(y)
                    fnameyr = fnameyearfunc(y)
                    # Create unique trajectory filename.
                    trajname = f"{basename}_ALT{'{:04}'.format(a)}m.{fnameyr}{m:02}{d:02}{hr:02}"
                    final_trajpath = os.path.join(traj_output_dir, trajname)
                    
                    # Remove any existing CONTROL and temporary trajectory files.
                    _try_to_remove(controlfname)
                    _try_to_remove(trajname)
                    _try_to_remove(final_trajpath)
                    
                    # Assemble the meteorology file list for the month of sim_dt.
                    meteofiles = _meteofinder(meteo_dir, ([4, 5], [1]), m, y, mon_dict, meteoyearfunc)
                    
                    # Populate CONTROL file with trajectory initialization data.
                    _populate_control(coordinates, controlyr, m, d, hr, a, meteo_dir, meteofiles, run, controlfname, trajname)
                    
                    # Call the HYSPLIT executable to generate the trajectory.
                    call(hysplit)
                    
                    # Generate reverse and/or clipped trajectories if requested.
                    if get_reverse:
                        _reversetraj_whilegen(trajname, run, hysplit, output_rdir, meteo_dir, meteofiles, controlfname)
                    if get_clipped:
                        _cliptraj(output_cdir, trajname)
                    
                    # Move the computed trajectory file to the output directory.
                    shutil.move(trajname, final_trajpath)
            current_dt += dt.timedelta(days=1)
    finally:
        os.chdir(cwd)

def _reversetraj_whilegen(trajname, run, hysplit, output_rdir, meteo_dir,
                          meteofiles, controlfname):
    """
    Calculate reverse trajectory during main trajectory generation.

    Calculates a new trajectory ('reverse trajectory') from the endpoint of the
    trajectory just calculated in the main generation sequence and running
    in the opposite direction.

    Parameters
    ----------
    trajname : string
        The file name of the just-calculated trajectory.  New backwards
        trajectory will be named ``trajname`` + 'REVERSE'
    run : int
        The length in hours of the trajectory simulation
    hysplit : string
        The location of the executable that calculates trajectories
    output_rdir : string
        The subdirectory for reverse trajectories.
    meteo_dir : string
        The location of the meteorology files
    meteofiles : string
        The list of meteorology files required to calculate the
        reverse trajectory.
    controlfname : string
        The name of the control file, which should be 'CONTROL'

    """
    # Initialize name, path
    reversetrajname = trajname + 'REVERSE'
    final_rtrajpath = os.path.join(output_rdir, reversetrajname)

    with open(trajname) as traj:
        contents = traj.readlines()

        last_timepoint = -1

        # for multiline files, critical information is in contents[-2]
        if len(contents[-1]) < len(contents[-2]):
            last_timepoint = -2

        data = contents[last_timepoint].split()

    # Get reverse trajectory start information
    year = int(data[2])
    mon  = int(data[3])
    day  = int(data[4])
    hour = int(data[5])
    lat  = float(data[9])
    lon  = float(data[10])
    alt  = float(data[11])
    run  = run * -1

    # Sometimes start height is greater than 10000 m model top
    if alt >= 10000:
        alt = 9999

    # Always 2 digit year for CONTROL
    yr = '{:02}'.format(year)

    # Remove (if present) any existing CONTROL or temp files
    _try_to_remove(controlfname)
    _try_to_remove(reversetrajname)
    _try_to_remove(final_rtrajpath)

    # Populate control text
    _populate_control((lat, lon), yr, mon, day, hour, alt, meteo_dir,
                      meteofiles, run, controlfname, reversetrajname)

    # Call executable
    call(hysplit)

    # Move the trajectory file to the desired output directory
    shutil.move(reversetrajname, final_rtrajpath)

def _meteofinder(meteo_dir, meteo_bookends, mon, year, mon_dict,
                 meteoyearfunc):
    """
    Get list of meteorology files.

    Creates list of files in storage location ``meteo_dir`` that belong
    to the given month and year, plus the necessary files from previous
    and the next months (``meteo_bookends``).

    For successful meteofinding, separate different meteorology types into
    different folders and name weekly or semi-monthly files according to the
    following convention:
        *mon*YY*#
    where the * represents a Bash wildcard.

    Parameters
    ----------
    meteo_dir : string
        Full or relative path to the location of the meteorology files.
    meteo_bookends : tuple of lists of ints
        To calculate a month of trajectories, files from the previous and next
        month must be included.  This indicates which file numbers from the
        previous month and which from the next month are necessary.
        The user is responsible for making sure the correct bookends for their
        trajectory length and meteorology file periods are provided.
    mon : int
        The integer representation of the current month.  Converted to a
        3-letter string to find meteorology files.
    year : int
        The integer representation of the current year.  Converted to a length
        2 string to find meteorology files.
    mon_dict : dictionary
        Dictionary keyed by month integer, with lists of [season, mon]
    meteoyearfunc : function
        Function that formats the year string to length 2 or 4 to identify
        appropriate meteorology files

    Returns
    -------
    meteofiles : list of strings
        List of strings representing the names of the required
        meteorology files

    """
    # Current working directory set in generate_bulktraj() environment
    orig_dir = os.getcwd()

    # Initialize lists, count
    meteofiles = []
    file_number = -1

    # Get the strings that will match files for the previous, next,
    # and current months
    prv, nxt, now = _monyearstrings(mon, year, mon_dict, meteoyearfunc)

    # Change directory and walk through files
    try:
        os.chdir(meteo_dir)

        _, _, files = next(os.walk('.'))

        # Order of files to CONTROL doesn't matter
        for each_file in files:
            if fnmatch.fnmatch(each_file, now):
                meteofiles.append(each_file)
            elif fnmatch.fnmatch(each_file, prv):
                if int(each_file[file_number]) in meteo_bookends[0]:
                    meteofiles.append(each_file)
            elif fnmatch.fnmatch(each_file, nxt):
                if int(each_file[file_number]) in meteo_bookends[1]:
                    meteofiles.append(each_file)

    finally:
        os.chdir(orig_dir)

    num_files = len(meteofiles)

    if num_files == 0:
        raise OSError('0 files found for month/year %(mon)d / %(year)d'
                      %{'mon': mon, 'year': year})
        
    if num_files > 12:
        print(meteofiles)
        raise OSError('%(f)d files found for month/year %(mon)d / %(year)d.'\
                      '  Maximum 12 allowed.  If wrong years are included, '\
                      'identify files by 4 digit years (meteoyr_2digits=True).'\
                      '  May require renaming meteorology files.'
                      %{'f': num_files, 'mon': mon, 'year': year})

    return meteofiles

def _year2string(year):
    """
    Helper function, takes a four digit integer year, makes a length-2 string.

    Parameters
    ----------
    year : int
        The year.

    Returns
    -------
    Length-2 string representation of ``year``

    """
    return '{0:02}'.format(year % 100)

def _monyearstrings(mon, year, mon_dict, meteoyearfunc):
    """
    Increment the months and potentially the years.

    Assemble the strings that will allow ``_meteofinder`` to get correct files.

    Parameters
    ----------
    mon : int
        Integer representation of the month
    year : int
        Integer representation of the year
    mon_dict : dictionary
        Dictionary keyed by month integer, with lists of [season, mon]
    meteoyearfunc : function
        Function that formats the year string to length 2 or 4 to identify
        appropriate meteorology files

    Returns
    -------
    prv : string
        Signature for gathering the meteorology files for the previous month
    nxt : string
        Signature for gathering the meteorology files for the next month
    now : string
        Signature for gathering the meteorology files for the current month

    """
    next_year = year
    prev_year = year

    next_mon = mon + 1
    prev_mon = mon - 1

    if prev_mon == 0:
        prev_mon = 12
        prev_year = year - 1
    if next_mon == 13:
        next_mon = 1
        next_year = year + 1

    w = '*'

    prv = w + mon_dict[prev_mon][1] + w + meteoyearfunc(prev_year) + w
    nxt = w + mon_dict[next_mon][1] + w + meteoyearfunc(next_year) + w
    now = w + mon_dict[mon][1] + w + meteoyearfunc(year) + w

    return prv, nxt, now

def _mondict(n_hem=True):
    """
    Get a dictionary of season and month string.

    Parameters
    ----------
    n_hem : Boolean
        Default True.  Indicates hemisphere of parcel launch and thus
        actual season.

    Returns
    -------
    season_month_dict : dictionary
        Dictionary keyed by month integer, with lists of [season, mon]

    """
    if n_hem:
        season_month_dict = {12: ['winter', 'dec'],
                             1 : ['winter', 'jan'],
                             2 : ['winter', 'feb'],
                             3 : ['spring', 'mar'],
                             4 : ['spring', 'apr'],
                             5 : ['spring', 'may'],
                             6 : ['summer', 'jun'],
                             7 : ['summer', 'jul'],
                             8 : ['summer', 'aug'],
                             9 : ['autumn', 'sep'],
                             10: ['autumn', 'oct'],
                             11: ['autumn', 'nov']}
    else:
        season_month_dict = {12: ['summer', 'dec'],
                             1 : ['summer', 'jan'],
                             2 : ['summer', 'feb'],
                             3 : ['autumn', 'mar'],
                             4 : ['autumn', 'apr'],
                             5 : ['autumn', 'may'],
                             6 : ['winter', 'jun'],
                             7 : ['winter', 'jul'],
                             8 : ['winter', 'aug'],
                             9 : ['spring', 'sep'],
                             10: ['spring', 'oct'],
                             11: ['spring', 'nov']}

    return season_month_dict

