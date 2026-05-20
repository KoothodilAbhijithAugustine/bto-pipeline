#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BTO LEVEL-1B CALIBRATION PIPELINE
Takes raw L1a FITS files as input, applies physical calibrations, 
and generates science-ready L1b FITS files.
"""

import sys
import os
import glob
import numpy as np
from astropy.io import fits

# =============================================================================
# 1. MISSION CALIBRATION CONSTANTS
# =============================================================================
CAL_HK = {"EXT": (0.1022, -275.66), "DET": (0.5432, -268.3)}
E_SLOPE, E_INTERCEPT = 0.84098, -1.64736   
MAX_CHANNELS = 4096

# =============================================================================
# 2. UTILITIES
# =============================================================================
def get_ebounds_hdu():
    """Generates the OGIP EBOUNDS extension mapping PI to Energy (keV)."""
    channels = np.arange(MAX_CHANNELS, dtype=np.int16)
    e_min = np.maximum(0, E_SLOPE * (channels - 0.5) + E_INTERCEPT)
    e_max = np.maximum(0, E_SLOPE * (channels + 0.5) + E_INTERCEPT)
    
    cols = [
        fits.Column(name='CHANNEL', format='1I', array=channels), 
        fits.Column(name='E_MIN', format='1E', array=e_min, unit='keV'), 
        fits.Column(name='E_MAX', format='1E', array=e_max, unit='keV')
    ]
    
    hdu = fits.BinTableHDU.from_columns(cols, name='EBOUNDS')
    
    # Aligned with fits_parser.py COSI/BTO instrument standards
    hdu.header.update({
        'EXTNAME': 'EBOUNDS', 
        'TELESCOP': 'COSI', 
        'INSTRUME': 'BTO', 
        'HDUCLASS': 'OGIP', 
        'HDUCLAS1': 'RESPONSE', 
        'HDUCLAS2': 'EBOUNDS', 
        'CHANTYPE': 'PI', 
        'DETCHANS': MAX_CHANNELS
    })
    return hdu

# =============================================================================
# 3. CORE PROCESSING LOGIC
# =============================================================================
def process_l1a_to_l1b(l1a_filepath, output_dir):
    """Parses an L1a FITS file, applies calibrations, and writes an L1b FITS file."""
    
    # Keep the exact same standard filename, just saved into the new L1b output directory
    basename = os.path.basename(l1a_filepath)
    l1b_filepath = os.path.join(output_dir, basename)
    
    os.makedirs(output_dir, exist_ok=True)

    try:
        # memmap=False prevents file-locking issues during read
        with fits.open(l1a_filepath, memmap=False) as hdul:
            primary_hdr = hdul[0].header.copy()
            data_hdu = hdul[1]
            data_hdr = data_hdu.header.copy()
            data = data_hdu.data
            
            ext_name = data_hdr.get('EXTNAME', 'UNKNOWN')
            new_cols = []
            hdu_class = 'DATA'
            
            # Enforce Unix-Time metadata standards on the copied headers
            data_hdr.update({
                'TIMESYS': 'UNIX', 
                'MJDREFI': 40587, 
                'MJDREFF': 0.0
            })

            # ---------------------------------------------------------
            # HOUSEKEEPING (HK)
            # ---------------------------------------------------------
            if ext_name == 'HK':
                time_arr = data['TIME']
                pkt_arr = data['PKT_CNT']
                t_ext_raw = data['t_ext_raw']
                t_det1_raw = data['t_det1_raw']
                
                # Apply Calibrations
                t_ext_cal = (t_ext_raw * CAL_HK['EXT'][0]) + CAL_HK['EXT'][1]
                t_det1_cal = (t_det1_raw * CAL_HK['DET'][0]) + CAL_HK['DET'][1]
                
                # Ensure 1D (64-bit float) for TIME and 1E (32-bit float) for temps
                new_cols = [
                    fits.Column(name='TIME', format='1D', array=time_arr),
                    fits.Column(name='PKT_CNT', format='1I', array=pkt_arr),
                    fits.Column(name='t_ext_raw', format='1I', array=t_ext_raw),
                    fits.Column(name='t_det1_raw', format='1I', array=t_det1_raw),
                    fits.Column(name='t_ext', format='1E', array=t_ext_cal.astype(np.float32)),
                    fits.Column(name='t_det1', format='1E', array=t_det1_cal.astype(np.float32))
                ]
                hdu_class = 'HK'

            # ---------------------------------------------------------
            # PHOTON EVENTS (EVT)
            # ---------------------------------------------------------
            elif ext_name == 'EVENTS':
                time_arr = data['TIME']
                pkt_arr = data['PKT_CNT']
                pha_arr = data['PHA']
                
                # In BTO, PI maps 1:1 to PHA before EBOUNDS calibration
                new_cols = [
                    fits.Column(name='TIME', format='1D', array=time_arr),
                    fits.Column(name='PKT_CNT', format='1I', array=pkt_arr),
                    fits.Column(name='PI', format='1I', array=pha_arr)
                ]
                hdu_class = 'EVENTS'

            # ---------------------------------------------------------
            # HISTOGRAMS / LIGHTCURVES (LC)
            # ---------------------------------------------------------
            elif ext_name == 'SPECTRUM':
                time_arr = data['TIME']
                pkt_arr = data['PKT_CNT']
                bins_arr = data['bins']
                
                # Cast integer bins to 32-bit float for rates
                rates_arr = bins_arr.astype(np.float32) 
                
                new_cols = [
                    fits.Column(name='TIME', format='1D', array=time_arr),
                    fits.Column(name='PKT_CNT', format='1I', array=pkt_arr),
                    fits.Column(name='bins', format='29J', array=bins_arr),
                    fits.Column(name='RATES', format='29E', array=rates_arr)
                ]
                hdu_class = 'LIGHTCURVE'

            else:
                print(f"[-] Skipping {basename}: Unknown extension type '{ext_name}'")
                return

            # Build the new FITS structure
            new_data_hdu = fits.BinTableHDU.from_columns(new_cols, header=data_hdr)
            new_data_hdu.header.update({'HDUCLAS1': hdu_class, 'HDUCLAS2': 'ALL'})
            
            # Inject EBOUNDS if this is an event file
            hdul_new = [fits.PrimaryHDU(header=primary_hdr), new_data_hdu]
            if ext_name == 'EVENTS':
                new_data_hdu.header.update({'CHANTYPE': 'PI', 'DETCHANS': MAX_CHANNELS})
                hdul_new.append(get_ebounds_hdu())
                
            fits.HDUList(hdul_new).writeto(l1b_filepath, overwrite=True)
            print(f"[+] L1b Generated: {basename}")

    except Exception as e:
        print(f"[!] Failed to process {l1a_filepath}: {e}")

# =============================================================================
# 4. EXECUTION
# =============================================================================
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python l1a_to_l1b.py <input_l1a_file_or_dir> <output_l1b_dir>")
    else:
        input_path = sys.argv[1]
        out_dir = sys.argv[2]
        
        if os.path.isfile(input_path):
            process_l1a_to_l1b(input_path, out_dir)
        elif os.path.isdir(input_path):
            # Recursively search for FITS files
            search_pattern = os.path.join(input_path, "**", "*.fits")
            files = glob.glob(search_pattern, recursive=True)
            
            print(f"[*] Found {len(files)} L1a FITS files. Processing...")
            for f in files:
                process_l1a_to_l1b(f, out_dir)
            print("[*] Calibration complete.")