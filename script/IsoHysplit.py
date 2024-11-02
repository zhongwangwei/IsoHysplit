import os
import sys
import subprocess
import re

from Mod_Namelist import *

def read_setup(filename):
    """Read SETUP section from namelist file and save to SETUP.CFG"""
    with open(filename, 'r') as f:
        content = f.read()
    
    # Extract SETUP section
    setup_match = re.search(r'&SETUP\s*(.*?)\s*/\s*$', content, re.MULTILINE | re.DOTALL)
    if not setup_match:
        raise ValueError("SETUP section not found in namelist file")
    
    setup_text = setup_match.group(1)
    
    # Parse parameters
    params = {}
    for line in setup_text.strip().split('\n'):
        line = line.strip()
        # Remove comments
        line = line.split('#')[0].strip()
        if '=' in line:
            key, value = line.split('=', 1)
            params[key.strip()] = value.strip()
    
    # Write SETUP.CFG
    with open('../working/SETUP.CFG', 'w') as f:
        for key, value in params.items():
            f.write(f"{value.strip(',')}\n")

def main():
    """Main function to process namelist file"""
    if len(sys.argv) != 2:
        print("Usage: python IsoHysplit.py <namelist_file>")
        sys.exit(1)
    
    namelist_file = sys.argv[1]
    if not os.path.exists(namelist_file):
        print(f"Error: File {namelist_file} not found")
        sys.exit(1)
        
    try:
        read_setup(namelist_file)
        print("Successfully created SETUP.CFG")
    except Exception as e:
        print(f"Error processing namelist file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
