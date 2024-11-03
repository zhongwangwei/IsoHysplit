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

    def __init__(self):
        """
        Initialize the NamelistReader with metadata and error settings.
        """
        self.name = 'namelist_read'
        self.version = '0.1'
        self.release = '0.1'
        self.date = 'Mar 2023'
        self.author = "Zhongwang Wei / zhongwang007@gmail.com"
        
        # Ignore all numpy warnings
        np.seterr(all='ignore')

    @staticmethod
    def strtobool(val: str) -> int:
        """
        Convert a string representation of truth to 1 (true) or 0 (false).

        Args:
            val (str): The string to convert.

        Returns:
            int: 1 for true values, 0 for false values.

        Raises:
            ValueError: If the input string is not a valid truth value.
        """
        val = val.lower()
        if val in ('y', 'yes', 't', 'true', 'on', '1'):
            return 1
        elif val in ('n', 'no', 'f', 'false', 'off', '0'):
            return 0
        else:
            raise ValueError(f"Invalid truth value: {val}")

    @staticmethod
    def select_variables(namelist: Dict[str, Any]) -> Dict[str, Any]:
        """
        Select variables from namelist if the value is truthy.

        Args:
            namelist (Dict[str, Any]): The namelist dictionary.

        Returns:
            Dict[str, Any]: A dictionary containing only the truthy values.
        """
        return {k: v for k, v in namelist.items() if v}

    def read_namelist(self, file_path: str) -> Dict[str, Dict[str, Any]]:
        """
        Read a namelist from a text file.

        Args:
            file_path (str): The path to the namelist file.

        Returns:
            Dict[str, Dict[str, Any]]: A nested dictionary representing the namelist structure.
        """
        namelist = {}
        current_dict = None

        def parse_value(key: str, value: str) -> Union[bool, int, float, list, str]:
            """
            Parse a string value into its appropriate type.

            Args:
                key (str): The key of the value being parsed.
                value (str): The string value to parse.

            Returns:
                Union[bool, int, float, list, str]: The parsed value.
            """
            value = value.strip()
            if key in ['suffix', 'prefix']:
                return value  # Return as string for suffix and prefix
            if value.lower() in ['true', 'false']:
                return bool(self.strtobool(value))
            elif value.replace('-', '', 1).isdigit():
                return int(value)
            elif value.replace('.', '', 1).replace('-', '', 1).isdigit():
                return float(value)
            elif ',' in value:
                return [v.strip() for v in value.split(',')]
            else:
                return value

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
                    value = value.split('#')[0].strip()  # Remove inline comments
                    current_dict[key] = parse_value(key, value)

        return namelist
    def get_namelist_attr(self, namelist_file: str) -> Any:
        """
        Get an attribute from the namelist.

        Args:
            namelist (Dict[str, Dict[str, Any]]): The namelist dictionary.
            key (str): The key to retrieve.

        Returns:
            Any: The value of the attribute, or None if not found.
        """
        kk = self.read_namelist(namelist_file)
        # Convert all values in kk to lists if they're not already
        for section in kk:
            for key, value in kk[section].items():
                if isinstance(value, bool):
                    kk[section][key] = value
                elif not isinstance(value, list):
                    kk[section][key] = [value]
        #
        
        return kk

    def write_setup_cfg(self, namelist: Dict[str, Dict[str, Any]]) -> None:
        """
        Write SETUP.CFG in namelist format with all required parameters.
        
        Args:
            namelist (Dict[str, Dict[str, Any]]): The namelist dictionary
        """
        # Define default values for all required parameters
        default_setup = {
            'tratio': 0.75,
            'delt': 0.0,
            'mgmin': 10,
            'khmax': 9999,
            'kmixd': 0,
            'kmsl': 0,
            'kagl': 1,
            'k10m': 1,
            'nstr': 0,
            'mhrs': 9999,
            'nver': 0,
            'tout': 60,
            'tm_pres': 1,
            'tm_tpot': 1,
            'tm_tamb': 1,
            'tm_rain': 1,
            'tm_mixd': 1,
            'tm_relh': 1,
            'tm_sphu': 1,
            'tm_mixr': 1,
            'tm_dswf': 1,
            'tm_terr': 1,
            'dxf': 1.00,
            'dyf': 1.00,
            'dzf': 0.01,
            'messg': 'MESSAGE'
        }
        
        # Update defaults with any provided values from namelist
        for key in default_setup:
            if key in namelist:
                value = namelist[key]
                if isinstance(value, list):
                    value = value[0]
                default_setup[key] = value
        
        os.makedirs('../working', exist_ok=True)
        with open('../working/SETUP.CFG', 'w') as f:
            f.write('&SETUP\n')
            for key, value in default_setup.items():
                if isinstance(value, str):
                    f.write(f" {key} = '{value}',\n")
                elif isinstance(value, float):
                    f.write(f" {key} = {value:.2f},\n")
                else:
                    f.write(f" {key} = {value},\n")
            f.write(' /')
