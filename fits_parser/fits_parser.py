#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BTO MASTER PIPELINE
Monitors or batch-processes CCSDS hex logs into L0, L1a, and L1b FITS archives.
Supports live daemon mode and selective processing tiers.
"""

import sys
import re
import datetime
import os
import time
import glob
import argparse
import numpy as np
from astropy.io import fits

# =============================================================================
# 1. MISSION CONSTANTS & OFFSETS
# =============================================================================
GPS_EPOCH = datetime.datetime(1980, 1, 6, 0, 0, 0, tzinfo=datetime.timezone.utc)
MET_EPOCH = datetime.datetime(2025, 1, 1, 0, 0, 0, tzinfo=datetime.timezone.utc)
TICKS_TO_SEC = 1.0 / 40_000_000.0
DWT_TICK_SEC = 1.0 / 300_000_000.0

CAL_HK = {"EXT": (0.1022, -275.66), "DET": (0.5432, -268.3)}
E_SLOPE, E_INTERCEPT = 0.84098, -1.64736   
MAX_CHANNELS = 4096      

PKT_SEQ_OFS, PKT_SEQ_LEN = 2, 2  # CCSDS Sequence Control Word
PKT_SEC_OFS, PKT_SEC_LEN = 6, 4
PKT_TKS_OFS, PKT_TKS_LEN = 10, 4
HK_T_BOARD_OFS, HK_T_DET_OFS, HK_T_LEN = 72, 74, 2
LC_BINS_START, LC_BINS_STEP, LC_NUM_BINS = 15, 2, 29
LC_SEC_OFS, LC_TKS_OFS, LC_TIME_LEN = 77, 81, 4
EVT_DATA_START, EVT_WORD_LEN = 20, 8

FLUSH_THRESHOLD = 1000  # Number of packets to buffer before flushing to disk

# =============================================================================
# 2. TIME & ROUTING UTILITIES
# =============================================================================
def pps_to_utc(pps_sec: int, ticks: int = 0) -> datetime.datetime:
    """
    Converts coarse seconds and fine clock ticks into an absolute UTC datetime.
    
    Args:
        pps_sec (int): Seconds elapsed since the GPS Epoch.
        ticks (int, optional): Sub-second ticks from the 40 MHz FPGA clock. Defaults to 0.
        
    Returns:
        datetime.datetime: The calculated absolute UTC timestamp.
    """
    return GPS_EPOCH + datetime.timedelta(seconds=(float(pps_sec) + ticks * TICKS_TO_SEC))

def get_met(dt_utc: datetime.datetime) -> float:
    """
    Converts a UTC datetime into Mission Elapsed Time (MET) in seconds.
    
    Args:
        dt_utc (datetime.datetime): The absolute UTC timestamp.
        
    Returns:
        float: Total seconds elapsed since the mission MET_EPOCH (Jan 1, 2025).
    """
    return (dt_utc - MET_EPOCH).total_seconds()

def get_day_fraction(dt_utc: datetime.datetime) -> str:
    """
    Calculates the millisecond fraction of the current day for use in file naming.
    
    Args:
        dt_utc (datetime.datetime): The absolute UTC timestamp.
        
    Returns:
        str: A zero-padded, 3-digit string representing the fraction of the day.
    """
    sec_mid = (dt_utc - dt_utc.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()
    return f"{int((sec_mid / 86400.0) * 1000):03d}"

class ArchiveRouter:
    """
    Manages the creation of directories and generation of standard FITS file paths 
    for the BTO Triple-Tier Data Archive based on packet types and timestamps.
    """
    def __init__(self, root="BTO_Data_Archive"):
        """
        Initializes the ArchiveRouter.
        
        Args:
            root (str, optional): The root directory for the generated archive. Defaults to "BTO_Data_Archive".
        """
        self.root = root

    def get_path(self, tier: str, apid: int, dt: datetime.datetime, tid: int = 0) -> str:
        """
        Generates the standard directory structure and filename for a given telemetry packet.
        Creates the target directories on the disk if they do not already exist.
        
        Args:
            tier (str): The data processing tier ("L0", "L1a", "L1b").
            apid (int): The Application Process Identifier (0xD6, 0xD7, 0xD8).
            dt (datetime.datetime): The UTC timestamp associated with the data.
            tid (int, optional): The Deadtime Transaction ID (for event files). Defaults to 0.
            
        Returns:
            str: The full path where the resulting FITS or Binary file should be saved.
        """
        yyyy, mm, dd, hh = dt.strftime("%Y"), dt.strftime("%m"), dt.strftime("%d"), dt.strftime("%H")
        yymmdd = dt.strftime("%y%m%d")
        if tier == "L0":
            folder = {0xD6: "D6_Spectrum", 0xD7: "D7_Event_packages", 0xD8: "D8_House_keeping"}.get(apid, "Other")
            base = os.path.join(self.root, "L0_binary", yyyy, mm, dd, folder)
            fname = f"bto_{yymmdd}_{folder.lower()}.bin"
        else:
            base = os.path.join(self.root, f"{tier}_fits", yyyy, mm, dd, hh)
            if apid == 0xD6:   fname = f"cs{yymmdd}bto_lc.fits"
            elif apid == 0xD8: fname = f"cs{yymmdd}bto_hk.fits"
            elif apid == 0xD7: fname = f"cs{yymmdd}{get_day_fraction(dt)}{tid:05d}_bto_evt.fits"
            else:              fname = f"cs{yymmdd}_misc.fits"
        os.makedirs(base, exist_ok=True)
        return os.path.join(base, fname)

# =============================================================================
# 3. PACKET DECODER
# =============================================================================
def decode_packet(ccsds: bytes, apid: int, d7_state: dict):
    """
    Extracts and decodes binary payloads from CCSDS packets into structured Python dictionaries.
    Handles high-precision Digital Watch Timer (DWT) synchronization for photon events.
    
    Supports:
    - 0x0D8: Housekeeping Data (HK)
    - 0x0D6: Histograms / Lightcurves (LC)
    - 0x0D7: Time-Tagged Photon Events (EVT)
    
    Args:
        ccsds (bytes): The raw payload bytes of the incoming CCSDS packet.
        apid (int): The Application Process Identifier extracted from the packet header.
        d7_state (dict): A mutable dictionary tracking the synchronization state 
                         between the 1 Hz PPS End Frames and the 300 MHz DWT.
                         
    Returns:
        dict or None: A dictionary containing parsed telemetry categorized into 'l1a' 
                      and 'l1b' formats, or None if the packet is invalid/pre-mission.
    """
    if len(ccsds) < 14: return None
    
    seq_ctrl = int.from_bytes(ccsds[PKT_SEQ_OFS : PKT_SEQ_OFS + PKT_SEQ_LEN], "big")
    pkt_count = seq_ctrl & 0x3FFF
    
    sec = int.from_bytes(ccsds[PKT_SEC_OFS : PKT_SEC_OFS + PKT_SEC_LEN], "big")
    tks = int.from_bytes(ccsds[PKT_TKS_OFS : PKT_TKS_OFS + PKT_TKS_LEN], "big")
    if apid == 0x0D6 and len(ccsds) >= 85: 
        sec = int.from_bytes(ccsds[LC_SEC_OFS : LC_SEC_OFS + LC_TIME_LEN], "big")
        tks = int.from_bytes(ccsds[LC_TKS_OFS : LC_TKS_OFS + LC_TIME_LEN], "big")
    if sec < 1420070400: return None 
    
    utc_dt = pps_to_utc(sec, tks)
    met = get_met(utc_dt)

    if apid == 0x0D8:
        t_board = int.from_bytes(ccsds[HK_T_BOARD_OFS : HK_T_BOARD_OFS + HK_T_LEN], "big")
        t_det   = int.from_bytes(ccsds[HK_T_DET_OFS : HK_T_DET_OFS + HK_T_LEN], "big")
        return {"type": "HK", "utc": utc_dt, "met": met,
                "l1a": [{"TIME": met, "PKT_CNT": pkt_count, "t_ext_raw": t_board, "t_det1_raw": t_det}],
                "l1b": [{"TIME": met, "PKT_CNT": pkt_count, "t_ext_raw": t_board, "t_det1_raw": t_det, "t_ext": (t_board * CAL_HK['EXT'][0]) + CAL_HK['EXT'][1], "t_det1": (t_det * CAL_HK['DET'][0]) + CAL_HK['DET'][1]}]}
    
    elif apid == 0x0D6:
        raw_bins = [int.from_bytes(ccsds[LC_BINS_START + (i*LC_BINS_STEP) : LC_BINS_START + (i*LC_BINS_STEP) + LC_BINS_STEP], "big") for i in range(LC_NUM_BINS)]
        return {"type": "LC", "utc": utc_dt, "met": met, 
                "l1a": [{"TIME": met, "PKT_CNT": pkt_count, "bins": raw_bins}], 
                "l1b": [{"TIME": met, "PKT_CNT": pkt_count, "bins": raw_bins, "RATES": [float(b) for b in raw_bins]}]}
    
    elif apid == 0x0D7:
        l1a, l1b, tid = [], [], 0
        ptr = EVT_DATA_START
        while ptr <= len(ccsds) - EVT_WORD_LEN:
            word = int.from_bytes(ccsds[ptr : ptr + EVT_WORD_LEN], "big")
            adc, ts_long = (word >> 48) & 0xFFFF, (word >> 32) & 0xFFFF
            
            if adc == 0xABCD and ts_long == 0xABCD:
                d7_state['pps_time'] = ((word & 0xFFFF) << 16) | ((word >> 16) & 0xFFFF)
                if ptr + (EVT_WORD_LEN * 3) <= len(ccsds):
                    w3 = int.from_bytes(ccsds[ptr + (EVT_WORD_LEN * 2) : ptr + (EVT_WORD_LEN * 3)], "big")
                    d7_state['dwt_at_last_pps'] = ((w3 & 0xFFFF) << 16) | ((w3 >> 16) & 0xFFFF)
                    ptr += (EVT_WORD_LEN * 3); continue
                    
            elif (adc & 0x7FFF) == 0x7FFF: 
                tid = (word & 0xFFFF) & 0x3FFF
                
            else:
                pha, dead_tid = adc & 0x0FFF, word & 0xFFFF 
                photon_dwt_32 = ((word >> 16) & 0xFFFFFF) << 8
                if d7_state['pps_time'] > 0:
                    delta_dwt = (photon_dwt_32 - d7_state['dwt_at_last_pps']) & 0xFFFFFFFF
                    if delta_dwt > 0x7FFFFFFF: delta_dwt -= 0x100000000 
                    abs_met = get_met(pps_to_utc(d7_state['pps_time']) + datetime.timedelta(seconds=delta_dwt * DWT_TICK_SEC))
                else: 
                    abs_met = met 
                
                l1a.append({"TIME": abs_met, "PKT_CNT": pkt_count, "PHA": pha, "DEADTIME": dead_tid})
                l1b.append({"TIME": abs_met, "PKT_CNT": pkt_count, "PI": pha})
                
            ptr += EVT_WORD_LEN
            
        return {"type": "EVT", "utc": utc_dt, "met": met, "l1a": l1a, "l1b": l1b, "tid": tid}
    return None

# =============================================================================
# 4. HEADER & FITS UTILITIES
# =============================================================================
def get_ebounds_hdu():
    """
    Generates a standard OGIP-compliant EBOUNDS FITS extension.
    This dynamically creates an energy calibration table mapping the 16-bit 
    Pulse Invariant (PI) channels to minimum and maximum continuous energy bounds (keV).
    
    Returns:
        fits.BinTableHDU: The constructed EBOUNDS binary table extension.
    """
    channels = np.arange(MAX_CHANNELS, dtype=np.int16)
    e_min = np.maximum(0, E_SLOPE * (channels - 0.5) + E_INTERCEPT)
    e_max = np.maximum(0, E_SLOPE * (channels + 0.5) + E_INTERCEPT)
    cols = [fits.Column(name='CHANNEL', format='1I', array=channels), fits.Column(name='E_MIN', format='1E', array=e_min, unit='keV'), fits.Column(name='E_MAX', format='1E', array=e_max, unit='keV')]
    hdu = fits.BinTableHDU.from_columns(cols, name='EBOUNDS')
    hdu.header.update({'EXTNAME': 'EBOUNDS', 'TELESCOP': 'COSI', 'INSTRUME': 'BTO', 'HDUCLASS': 'OGIP', 'HDUCLAS1': 'RESPONSE', 'HDUCLAS2': 'EBOUNDS', 'CHANTYPE': 'PI', 'DETCHANS': MAX_CHANNELS})
    return hdu

def inject_metadata(hdu, t_start, t_stop, utc_start, utc_stop, is_primary=False):
    """
    Injects mission-compliant metadata headers into a given FITS Header Data Unit (HDU).
    Ensures strict alignment with the BTO-PackageDefinition.xlsx specification.
    
    Args:
        hdu (fits.PrimaryHDU or fits.BinTableHDU): The target FITS HDU to modify in-place.
        t_start (float): The MET start time of the observation window.
        t_stop (float): The MET end time of the observation window.
        utc_start (datetime.datetime): The absolute UTC start time of the observation window.
        utc_stop (datetime.datetime): The absolute UTC end time of the observation window.
        is_primary (bool, optional): Determines if the metadata is formatted for the Primary 
                                     HDU (True) or a standard Extension HDU (False). Defaults to False.
    """
    obs_id = f"{utc_start.strftime('%y%m%d')}000t"
    header_data = {'TELESCOP': ('COSI', 'Telescope'), 'INSTRUME': ('BTO', 'Instrument'), 'OBS_ID': (obs_id, 'Observation ID'), 'DATE-OBS': (utc_start.strftime('%Y-%m-%dT%H:%M:%S'), 'Start'), 'DATE-END': (utc_stop.strftime('%Y-%m-%dT%H:%M:%S'), 'End')}
    if is_primary: header_data.update({'ORIGIN': ('UCB/SSL', 'Origin'), 'DATE': (datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%S'), 'Created'), 'CREATOR': ('BTO_LIVE_V1', 'Software')})
    else: header_data.update({'HDUCLASS': ('OGIP', 'Standard'), 'DATAMODE': ('NORMAL', 'Datamode'), 'OBSERVER': ('BTO_TEAM', 'PI'), 'OBJECT': ('CAL_SOURCE', 'Target'), 'TIMESYS': ('TT', 'Time System'), 'MJDREFI': (60676, 'MJD Ref'), 'MJDREFF': (0.0008007407407407, 'MJD offset'), 'TIMEREF': ('LOCAL', 'Ref Frame'), 'TASSIGN': ('SATELLITE', 'Time clock'), 'TIMEUNIT': ('s', 'Time unit'), 'TSTART': (t_start, 'Start MET'), 'TSTOP': (t_stop, 'Stop MET'), 'CLOCKAPP': ('F', 'Clock corr?')})
    for key, (val, comment) in header_data.items(): hdu.header[key] = (val, comment)

def flush_cache_to_disk(cache, tier):
    """
    Safely transfers decoded telemetry buffered in memory to persistent FITS files on disk.
    Appends new data rows if the file exists, or creates a new file with headers if it does not.
    
    Args:
        cache (dict): The in-memory buffer of decoded packets mapped by target filepath.
        tier (str): The processing tier label ("L1a" or "L1b") for terminal logging.
    """
    for path, data in list(cache.items()):
        if not data['met_list']: 
            del cache[path]; continue
            
        t_s, t_e = min(data['met_list']), max(data['met_list'])
        u_s, u_e = MET_EPOCH + datetime.timedelta(seconds=t_s), MET_EPOCH + datetime.timedelta(seconds=t_e)
        
        cols = []
        for k in data['rows'][0].keys():
            fmt = '1D' if k == 'TIME' else ('29E' if k == 'RATES' else ('29J' if k == 'bins' else ('1I' if k in ['PI', 'PHA', 'PKT_CNT'] else '1J')))
            arr = np.array([r[k] for r in data['rows']])
            if fmt == '1I': arr = arr.astype(np.int16)
            cols.append(fits.Column(name=k, format=fmt, array=arr))
            
        new_hdu = fits.BinTableHDU.from_columns(cols, name='EVENTS' if ('PHA' in data['rows'][0] or 'PI' in data['rows'][0]) else 'DATA')
        inject_metadata(new_hdu, t_s, t_e, u_s, u_e, is_primary=False)
        if tier == 'L1b' and 'PI' in data['rows'][0]: new_hdu.header.update({'HDUCLAS1': 'EVENTS', 'HDUCLAS2': 'ALL', 'CHANTYPE': 'PI', 'DETCHANS': MAX_CHANNELS})

        if os.path.exists(path):
            try:
                with fits.open(path, mode='update') as hdul:
                    old_hdu = hdul[1]
                    new_tstop = max(old_hdu.header.get('TSTOP', t_e), t_e)
                    merged_nrows = old_hdu.data.shape[0] + new_hdu.data.shape[0]
                    merged_hdu = fits.BinTableHDU.from_columns(old_hdu.columns, nrows=merged_nrows, name=old_hdu.name, header=old_hdu.header)
                    for colname in old_hdu.columns.names:
                        merged_hdu.data[colname][:] = np.concatenate([old_hdu.data[colname], new_hdu.data[colname]])
                    merged_hdu.header['TSTOP'] = new_tstop
                    hdul[1] = merged_hdu
                    hdul.flush()
            except Exception as e:
                print(f"Warning: Could not append to {path}. Error: {e}")
        else:
            hdul = [fits.PrimaryHDU(), new_hdu]
            inject_metadata(hdul[0], t_s, t_e, u_s, u_e, is_primary=True)
            if tier == 'L1b' and 'PI' in data['rows'][0]: hdul.append(get_ebounds_hdu())
            fits.HDUList(hdul).writeto(path, overwrite=True)
            
        print(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {tier} Updated: {os.path.basename(path)}")
        del cache[path] 

# =============================================================================
# 5. INGESTION ROUTINES (Live & Offline)
# =============================================================================
def read_input_stream(target_path, is_live):
    """
    Generator that dynamically yields lines from either a single file, an offline directory, 
    or continuously tails an active directory in real-time.
    
    Args:
        target_path (str): The file or directory path to monitor/process.
        is_live (bool): Flag indicating whether to run in continuous live mode.
        
    Yields:
        str: A single line of raw hex telemetry.
    """
    if is_live:
        if not os.path.isdir(target_path):
            raise ValueError("Live mode requires a directory path, not a file.")
        
        current_file = None
        f = None
        print(f"[*] Live Monitor Active on: {target_path}")
        
        while True:
            files = sorted(glob.glob(os.path.join(target_path, "*.txt")), key=os.path.getmtime)
            if not files:
                time.sleep(1)
                continue
                
            newest_file = files[-1]
            if current_file != newest_file:
                if f: f.close()
                current_file = newest_file
                f = open(current_file, "r")
                print(f"\n--- Tailing Log: {os.path.basename(current_file)} ---")
                
            line = f.readline()
            if not line:
                time.sleep(0.1) 
                continue
            yield line
    else:
        files_to_process = []
        if os.path.isfile(target_path):
            files_to_process = [target_path]
        elif os.path.isdir(target_path):
            files_to_process = sorted(glob.glob(os.path.join(target_path, "*.txt")), key=os.path.getmtime)
            
        print(f"[*] Offline Mode: Found {len(files_to_process)} file(s) to process.")
        
        for file in files_to_process:
            print(f"--- Processing {os.path.basename(file)} ---")
            with open(file, 'r') as f:
                for line in f:
                    yield line

# =============================================================================
# 6. MAIN PIPELINE LOOP
# =============================================================================
def run_pipeline(input_path: str, levels: list, is_live: bool):
    """
    Core execution loop applying level filters and managing data ingestion.
    
    Args:
        input_path (str): The target file or directory path.
        levels (list): The list of data tiers to process (e.g., ['l0', 'l1a', 'l1b']).
        is_live (bool): Whether to run as a continuous daemon or execute a batch run.
    """
    router = ArchiveRouter()
    cache_a, cache_b = {}, {}
    d7_timing_state = {'pps_time': 0, 'dwt_at_last_pps': 0}
    packet_count = 0
    
    print(f"[*] Target Levels: {', '.join(levels).upper()}")
    
    try:
        for line in read_input_stream(input_path, is_live):
            if not line.strip() or line.startswith("#"): continue
            try:
                stream = bytes.fromhex(re.sub(r"[^0-9A-Fa-f]", "", line))
                idx = 0
                while idx <= len(stream) - 10:
                    if stream[idx:idx+2] != b'\xeb\x90': 
                        idx += 1; continue
                    apid = (int.from_bytes(stream[idx+2:idx+4], "big") & 0x07FF)
                    blen = 512 if apid in [0xD6, 0xD7] else 2 + (int.from_bytes(stream[idx+6:idx+8], "big") + 7) + 2
                    
                    p = decode_packet(stream[idx+2:idx+blen-2], apid, d7_timing_state)
                    if p:
                        if 'l0' in levels:
                            with open(router.get_path("L0", apid, p['utc']), "ab") as f0: 
                                f0.write(stream[idx:idx+blen])
                        
                        if 'l1a' in levels:
                            pa = router.get_path("L1a", apid, p['utc'], p.get('tid',0))
                            if pa not in cache_a: cache_a[pa] = {'rows': [], 'met_list': []}
                            cache_a[pa]['rows'].extend(p['l1a'])
                            cache_a[pa]['met_list'].append(p['met'])
                        
                        if 'l1b' in levels:
                            pb = router.get_path("L1b", apid, p['utc'], p.get('tid',0))
                            if pb not in cache_b: cache_b[pb] = {'rows': [], 'met_list': []}
                            cache_b[pb]['rows'].extend(p['l1b'])
                            cache_b[pb]['met_list'].append(p['met'])
                        
                        packet_count += 1
                    idx += blen
            except Exception: continue
            
            if packet_count >= FLUSH_THRESHOLD:
                if 'l1a' in levels: flush_cache_to_disk(cache_a, "L1a")
                if 'l1b' in levels: flush_cache_to_disk(cache_b, "L1b")
                packet_count = 0
                
        print("\n[*] Input stream complete. Running final disk flush...")
        if 'l1a' in levels: flush_cache_to_disk(cache_a, "L1a")
        if 'l1b' in levels: flush_cache_to_disk(cache_b, "L1b")
        print("[*] Pipeline completed successfully.")
                
    except KeyboardInterrupt:
        print("\n[*] Keyboard Interrupt received. Flushing remaining data to disk...")
        if 'l1a' in levels: flush_cache_to_disk(cache_a, "L1a")
        if 'l1b' in levels: flush_cache_to_disk(cache_b, "L1b")
        print("[*] Pipeline terminated gracefully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BTO Telemetry Pipeline Tool")
    parser.add_argument("--input", required=True, help="Path to input raw hex file or directory")
    parser.add_argument("--level", default="all", choices=['l0', 'l1a', 'l1b', 'all'], help="Target output level (default: all)")
    parser.add_argument("--live", choices=['yes', 'no', 'true', 'false'], default='no', help="Run as a continuous live daemon (default: no)")

    args = parser.parse_args()

    # Parse Arguments
    is_live = args.live.lower() in ['yes', 'true']
    levels = ['l0', 'l1a', 'l1b'] if args.level.lower() == 'all' else [args.level.lower()]

    # Execute
    run_pipeline(args.input, levels, is_live)