import numpy as np
import os
import glob
import pandas as pd
import importlib
import re
import sys
from typing import Dict, Any, Tuple, List, Union

class NamelistReader:
    """
    A class for reading and processing namelist files.
    """

    def __init__(self, project_root: str = None):
        """
        Initialize the NamelistReader with metadata and error settings.
        """
        self.name = 'namelist_read'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        self.project_root = project_root if project_root else os.getcwd()
        np.seterr(all='ignore')

    def get_namelist_attr(self, namelist_file: str) -> Dict[str, Any]:
        """
        Reads, processes, and returns the fully structured namelist.
        """
        raw_nml = self._read_namelist_from_file(namelist_file)
        processed_nml = self._process_namelist(raw_nml)
        self.write_setup_cfg(processed_nml)
        print("Successfully created SETUP.CFG")
        return processed_nml

    def _read_namelist_from_file(self, file_path: str) -> Dict[str, Dict[str, Any]]:
        """
        Reads a namelist from a text file.
        """
        namelist = {}
        current_dict = None

        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if line.startswith('&'):
                    dict_name = line[1:]
                    current_dict = {}
                    namelist[dict_name] = current_dict
                elif line.startswith('/'):
                    current_dict = None
                elif current_dict is not None:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.split('#')[0].strip()
                    current_dict[key] = self._parse_value(key, value)
        return namelist

    def _process_namelist(self, nml: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """
        Processes the raw namelist dictionary to clean and structure the data.
        """
        if 'general' in nml:
            nml['general'] = self._process_general_section(nml['general'])
        if 'plot' in nml:
            nml['plot'] = self._process_plot_section(nml['plot'])
        
        return nml

    def _process_general_section(self, general_nml: Dict[str, Any]) -> Dict[str, Any]:
        """
        Processes the [general] section of the namelist.
        """
        # Ensure paths are single strings and absolute relative to project_root
        for key in ['working_dir', 'storage_dir', 'meteo_dir', 'isotope_dir', 'hysplit']:
            if isinstance(general_nml.get(key), list):
                general_nml[key] = general_nml[key][0]
            if general_nml.get(key):
                general_nml[key] = os.path.join(self.project_root, general_nml[key])

        # Process location
        location = general_nml.get('location', '0,0')
        if isinstance(location, list):
            general_nml['location'] = tuple(float(c) for c in location)
        else:
            general_nml['location'] = tuple(float(c) for c in location.split(','))

        # Process startdate and enddate
        startdate = general_nml.get('startdate', '20000101')
        general_nml['startdate'] = str(startdate[0] if isinstance(startdate, list) else startdate)
        enddate = general_nml.get('enddate', '20000101')
        general_nml['enddate'] = str(enddate[0] if isinstance(enddate, list) else enddate)

        # Process runtime
        runtime = general_nml.get('runtime', -240)
        general_nml['runtime'] = int(runtime[0] if isinstance(runtime, list) else runtime)

        # Process altitudes
        altitudes = general_nml.get('altitudes', '1000')
        if isinstance(altitudes, list):
            general_nml['altitudes'] = [int(alt) for alt in altitudes]
        else:
            general_nml['altitudes'] = [int(alt.strip()) for alt in str(altitudes).split(',')]
        
        os.makedirs(general_nml['storage_dir'], exist_ok=True)
        return general_nml

    def _process_plot_section(self, plot_nml: Dict[str, Any]) -> Dict[str, Any]:
        """
        Processes the [plot] section of the namelist.
        """
        def get_coord(value):
            if isinstance(value, list):
                return float(value[0])
            return float(value)

        maxlon = get_coord(plot_nml.get('maxlon', 180))
        minlon = get_coord(plot_nml.get('minlon', -180))
        maxlat = get_coord(plot_nml.get('maxlat', 90))
        minlat = get_coord(plot_nml.get('minlat', -90))
        plot_nml['mapcorners'] = [minlon, minlat, maxlon, maxlat]
        return plot_nml

    def _parse_value(self, key: str, value: str) -> Union[bool, int, float, list, str]:
        """
        Parse a string value into its appropriate type.
        """
        value = value.strip()
        if key in ['suffix', 'prefix']:
            return value
        if value.lower() in ['true', 'false']:
            return self._strtobool(value)
        elif value.replace('-', '', 1).isdigit():
            return int(value)
        elif value.replace('.', '', 1).replace('-', '', 1).isdigit():
            return float(value)
        elif ',' in value:
            return [v.strip() for v in value.split(',')]
        else:
            return value

    @staticmethod
    def _strtobool(val: str) -> bool:
        """
        Convert a string representation of truth to a boolean.
        """
        return val.lower() in ('y', 'yes', 't', 'true', 'on', '1')

    def write_setup_cfg(self, namelist: Dict[str, Dict[str, Any]]) -> None:
        """
        Write SETUP.CFG in namelist format with all required parameters.
        """
        default_setup = {
            'tratio': 0.75, 'delt': 0.0, 'mgmin': 10, 'khmax': 9999,
            'kmixd': 0, 'kmsl': 0, 'kagl': 1, 'k10m': 1, 'nstr': 0,
            'mhrs': 9999, 'nver': 0, 'tout': 60, 'tm_pres': 1,
            'tm_tpot': 1, 'tm_tamb': 1, 'tm_rain': 1, 'tm_mixd': 1,
            'tm_relh': 1, 'tm_sphu': 1, 'tm_mixr': 1, 'tm_dswf': 1,
            'tm_terr': 1, 'dxf': 1.00, 'dyf': 1.00, 'dzf': 0.01,
            'messg': 'MESSAGE'
        }
        
        setup_values = namelist.get('setup', {})
        for key in default_setup:
            if key in setup_values:
                value = setup_values[key]
                default_setup[key] = value[0] if isinstance(value, list) else value
        
        working_dir = namelist.get('general', {}).get('working_dir', '../working')
        os.makedirs(working_dir, exist_ok=True)
        with open(os.path.join(working_dir, 'SETUP.CFG'), 'w') as f:
            f.write('&SETUP\n')
            for key, value in default_setup.items():
                if isinstance(value, str):
                    f.write(f" {key} = '{value}',\n")
                elif isinstance(value, float):
                    f.write(f" {key} = {value:.2f},\n")
                else:
                    f.write(f" {key} = {value},\n")
            f.write(' /\n')