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

class UpdateNamelist(NamelistReader):
    def __init__(self, main_nl: Dict[str, Any]):
        # Initialize with general settings
        self.__dict__.update(main_nl['general'])
        tmp = self._read_source_namelist(ref_nml, evaluation_item, ref_source, 'ref')

        for evaluation_item in evaluation_items:
            self._process_evaluation_item(evaluation_item, sim_nml, ref_nml)

    def _process_evaluation_item(self, evaluation_item: str, sim_nml: Dict[str, Any], ref_nml: Dict[str, Any]):
        """Process a single evaluation item for both reference and simulation data."""
        sim_sources = self._ensure_list(sim_nml['general'][f'{evaluation_item}_sim_source'])
        ref_sources = self._ensure_list(ref_nml['general'][f'{evaluation_item}_ref_source'])

        # Process reference sources
        for ref_source in ref_sources:
            self._process_ref_source(evaluation_item, ref_source, ref_nml)

        # Process simulation sources
        for sim_source in sim_sources:
            self._process_sim_source(evaluation_item, sim_source, sim_nml)

    @staticmethod
    def _ensure_list(value):
        """Ensure the given value is a list."""
        return [value] if isinstance(value, str) else value

    def _process_ref_source(self, evaluation_item: str, ref_source: str, ref_nml: Dict[str, Any]):
        """Process a single reference source for an evaluation item."""
        # Read the namelist for this reference source
        #print all the attributes
        tmp = self._read_source_namelist(ref_nml, evaluation_item, ref_source, 'ref')
        # Initialize the evaluation item dictionary if it doesn't exist
        ref_nml.setdefault(evaluation_item, {})
        # Process each attribute for the reference source
        attributes = [
            'data_type', 'data_groupby', 'tim_res', 'grid_res', 'syear', 'eyear', 'dir',
            'varname', 'varunit', 'suffix', 'prefix'
        ]
        for attr in attributes:
            self._set_attribute(ref_nml, evaluation_item, ref_source, attr, tmp, 'ref')
        # Special handling for station data
        if ref_nml[evaluation_item][f'{ref_source}_data_type'] == 'stn':
            self._set_attribute(ref_nml, evaluation_item, ref_source, 'fulllist', tmp, 'ref')
            try:
                ref_nml[evaluation_item][f'{ref_source}_max_uparea'] = tmp[evaluation_item]['max_uparea']
                ref_nml[evaluation_item][f'{ref_source}_min_uparea'] = tmp[evaluation_item]['min_uparea']
            except KeyError:
                pass

    def _process_sim_source(self, evaluation_item: str, sim_source: str, sim_nml: Dict[str, Any]):
        """Process a single simulation source for an evaluation item."""
        # Read the namelist for this simulation source
        tmp = self._read_source_namelist(sim_nml, evaluation_item, sim_source, 'sim')

        # Initialize the evaluation item dictionary if it doesn't exist
        sim_nml.setdefault(evaluation_item, {})

        # Process each attribute for the simulation source
        attributes = [
            'data_type', 'data_groupby', 'tim_res', 'grid_res', 'syear', 'eyear',
            'suffix', 'prefix', 'model', 'varname', 'varunit', 'dir'
        ]
        for attr in attributes:
            self._set_attribute(sim_nml, evaluation_item, sim_source, attr, tmp, 'sim')

        # Special handling for station data
        if sim_nml[evaluation_item][f'{sim_source}_data_type'] == 'stn':
            self._set_attribute(sim_nml, evaluation_item, sim_source, 'fulllist', tmp, 'sim')

    def _read_source_namelist(self, nml: Dict[str, Any], evaluation_item: str, source: str, source_type: str):
        """Read the namelist for a given source."""
        try:
            return self.read_namelist(nml[evaluation_item][f"{source}"])
        except:
            return self.read_namelist(nml['def_nml'][f"{source}"])

    def _set_attribute(self, nml: Dict[str, Any], evaluation_item: str, source: str, attr: str, tmp: Dict[str, Any], source_type: str):
        """Set an attribute for a source in the namelist."""
        key = f'{source}_{attr}'
        try:
            nml[evaluation_item][key] = tmp[evaluation_item][attr]
        except KeyError:
            try:
                nml[evaluation_item][key] = tmp['general'][attr]
            except KeyError:
                if attr == 'dir':
                    self._set_dir_attribute(nml, evaluation_item, source, tmp, source_type)
                elif attr in ['model', 'varname', 'varunit']:
                    self._set_model_attribute(nml, evaluation_item, source, attr, tmp)
                else:
                    print(f"Warning: {attr} is missing in namelist for {evaluation_item} - {source}")
                    nml[evaluation_item][key] = None  # Set to None if missing

    def _set_dir_attribute(self, nml: Dict[str, Any], evaluation_item: str, source: str, tmp: Dict[str, Any], source_type: str):
        """Set the directory attribute for a source."""
        try:
            root_dir = tmp['general']['root_dir']
            sub_dir = tmp[evaluation_item]['sub_dir']
            nml[evaluation_item][f'{source}_dir'] = os.path.join(root_dir, sub_dir)
        except KeyError:
            try:
                nml[evaluation_item][f'{source}_dir'] = tmp['general']['root_dir']
            except KeyError:
                print("dir is missing in namelist")

    def _set_model_attribute(self, nml: Dict[str, Any], evaluation_item: str, source: str, attr: str, tmp: Dict[str, Any]):
        """Set model-related attributes for a simulation source."""
        try:
            model_nml = self.read_namelist(tmp['general']['model_namelist'])
            try:
                nml[evaluation_item][f'{source}_{attr}'] = model_nml['general'][attr]
            except KeyError:
                try:
                    nml[evaluation_item][f'{source}_{attr}'] = model_nml[evaluation_item][attr]
                except KeyError:
                    print(f"{attr} is missing in namelist")
        except KeyError:
            print(f"{attr} is missing in namelist")




