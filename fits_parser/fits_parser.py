#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BTO MASTER PIPELINE (CLI EDITION)
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

GPS_EPOCH = datetime.datetime(1980, 1, 6, 0, 0, 0, tzinfo=datetime.timezone.utc)
MET_EPOCH = datetime.datetime(2025, 1, 1, 0, 0, 0, tzinfo=datetime.timezone.utc)
TICKS_TO_SEC = 1.0 / 40_000_000.0
DWT_TICK_SEC = 1.0 / 300_000_000.0

CAL_HK = {"EXT": (0.1022, -275.66), "DET": (0.5432, -268.3)}
E_SLOPE, E_INTERCEPT = 0.84098, -1.64736   
MAX_CHANNELS = 4096      

PKT_SEQ_OFS, PKT_SEQ_LEN = 2, 2  
PKT_SEC_OFS, PKT_SEC_LEN = 6, 4
PKT_TKS_OFS, PKT_TKS_LEN = 10, 4
HK_T_BOARD_OFS, HK_T_DET_OFS, HK_T_LEN = 72, 74, 2
LC_BINS_START, LC_BINS_STEP, LC_NUM_BINS = 15, 2, 29
LC_SEC_OFS, LC_TKS_OFS, LC_TIME_LEN = 77, 81, 4
EVT_DATA_START, EVT_WORD_LEN = 20, 8

FLUSH_THRESHOLD = 1000  

def pps_to_utc(pps_sec: int, ticks: int = 0) -> datetime.datetime:
    return GPS_EPOCH + datetime.timedelta(seconds=(float(pps_sec) + ticks * TICKS_TO_SEC))

def get_met(dt_utc: datetime.datetime) -> float:
    return (dt_utc - MET_EPOCH).total_seconds()

class ArchiveRouter:
    def __init__(self, root="BTO_Data_Archive"):
        self.root = root

    def get_path(self, tier: str, apid: int, dt: datetime.datetime, tid: int = 0) -> str:
        yyyy, mm, dd = dt.strftime("%Y"), dt.strftime("%m"), dt.strftime("%d")
        yymmdd = dt.strftime("%y%m%d")
        
        if tier == "L0":
            folder = {0xD6: "D6_Spectrum", 0xD7: "D7_Event_packages", 0xD8: "D8_House_keeping"}.get(apid, "Other")
            base = os.path.join(self.root, "L0_binary", yyyy, mm, dd, folder)
            fname = f"bto_{yymmdd}_{folder.lower()}.bin"
        else:
            base = os.path.join(self.root, f"{tier}_fits", yyyy, mm, dd)
            if apid == 0xD6:   fname = f"cs{yymmdd}bto_lc.fits"
            elif apid == 0xD8: fname = f"cs{yymmdd}bto_hk.fits"
            elif apid == 0xD7: fname = f"cs{yymmdd}_{tid:05d}_bto_evt.fits"
            else:              fname = f"cs{yymmdd}_misc.fits"
            
        os.makedirs(base, exist_ok=True)
        return os.path.join(base, fname)

def decode_packet(ccsds: bytes, apid: int, d7_state: dict):
    if len(ccsds) < 14: return None
    
    pkt_count = int.from_bytes(ccsds[PKT_SEQ_OFS : PKT_SEQ_OFS + PKT_SEQ_LEN], "big") & 0x3FFF
    sec = int.from_bytes(ccsds[PKT_SEC_OFS : PKT_SEC_OFS + PKT_SEC_LEN], "big")
    tks = int.from_bytes(ccsds[PKT_TKS_OFS : PKT_TKS_OFS + PKT_TKS_LEN], "big")
    if apid == 0x0D6 and len(ccsds) >= 85: 
        sec = int.from_bytes(ccsds[LC_SEC_OFS : LC_SEC_OFS + LC_TIME_LEN], "big")
        tks = int.from_bytes(ccsds[LC_TKS_OFS : LC_TKS_OFS + LC_TIME_LEN], "big")
    if sec < 1420070400: return None 
    
    utc_dt = pps_to_utc(sec, tks)
    unix_time = utc_dt.timestamp()
    
    base_data = {"apid": apid, "pkt_count": pkt_count, "utc": utc_dt, "met": get_met(utc_dt)}

    if apid == 0x0D8:
        t_board = int.from_bytes(ccsds[HK_T_BOARD_OFS : HK_T_BOARD_OFS + HK_T_LEN], "big")
        t_det   = int.from_bytes(ccsds[HK_T_DET_OFS : HK_T_DET_OFS + HK_T_LEN], "big")
        base_data.update({
            "type": "HK",
            "l1a": [{"TIME": unix_time, "PKT_CNT": pkt_count, "t_ext_raw": t_board, "t_det1_raw": t_det}],
            "l1b": [{"TIME": unix_time, "PKT_CNT": pkt_count, "t_ext_raw": t_board, "t_det1_raw": t_det, "t_ext": (t_board * CAL_HK['EXT'][0]) + CAL_HK['EXT'][1], "t_det1": (t_det * CAL_HK['DET'][0]) + CAL_HK['DET'][1]}]
        })
        return base_data
        
    elif apid == 0x0D6:
        raw_bins = [int.from_bytes(ccsds[LC_BINS_START + (i*LC_BINS_STEP) : LC_BINS_START + (i*LC_BINS_STEP) + 2], "big") for i in range(LC_NUM_BINS)]
        base_data.update({
            "type": "LC",
            "l1a": [{"TIME": unix_time, "PKT_CNT": pkt_count, "bins": raw_bins}], 
            "l1b": [{"TIME": unix_time, "PKT_CNT": pkt_count, "bins": raw_bins, "RATES": [float(b) for b in raw_bins]}]
        })
        return base_data
        
    elif apid == 0x0D7:
        l1a, l1b = [], []
        ptr = EVT_DATA_START
        while ptr <= len(ccsds) - EVT_WORD_LEN:
            word = int.from_bytes(ccsds[ptr : ptr + EVT_WORD_LEN], "big")
            adc = (word >> 48) & 0xFFFF
            
            if adc in (0xABCD, 0x0123):
                d7_state['pps_time'] = ((word & 0xFFFF) << 16) | ((word >> 16) & 0xFFFF)
                if ptr + 24 <= len(ccsds):
                    w3 = int.from_bytes(ccsds[ptr+16 : ptr+24], "big")
                    d7_state['dwt_at_last_pps'] = ((w3 & 0xFFFF) << 16) | ((w3 >> 16) & 0xFFFF)
                    ptr += 24; continue
            elif (adc & 0x7FFF) == 0x7FFF: 
                d7_state['current_tid'] = (word & 0xFFFF) & 0x3FFF
            else:
                pha, dead_tid = adc & 0x0FFF, word & 0xFFFF 
                dwt = ((word >> 16) & 0xFFFFFF) << 8
                dt = (dwt - d7_state['dwt_at_last_pps']) & 0xFFFFFFFF
                abs_utc = pps_to_utc(d7_state['pps_time']) + datetime.timedelta(seconds=(dt if dt < 0x7FFFFFFF else dt - 0x100000000) * DWT_TICK_SEC)
                evt_unix_time = abs_utc.timestamp()
                
                
                active_tid = d7_state.get('current_tid', 0)
                l1a.append({"TIME": evt_unix_time, "PKT_CNT": pkt_count, "PHA": pha, "DEADTIME": dead_tid, "_TID": active_tid})
                l1b.append({"TIME": evt_unix_time, "PKT_CNT": pkt_count, "PI": pha, "_TID": active_tid})
                
            ptr += EVT_WORD_LEN
            
        base_data.update({"type": "EVT", "l1a": l1a, "l1b": l1b})
        return base_data
    return None


def get_ebounds_hdu():
    channels = np.arange(MAX_CHANNELS, dtype=np.int16)
    e_min = np.maximum(0, E_SLOPE * (channels - 0.5) + E_INTERCEPT)
    e_max = np.maximum(0, E_SLOPE * (channels + 0.5) + E_INTERCEPT)
    cols = [fits.Column(name='CHANNEL', format='1I', array=channels), fits.Column(name='E_MIN', format='1E', array=e_min, unit='keV'), fits.Column(name='E_MAX', format='1E', array=e_max, unit='keV')]
    hdu = fits.BinTableHDU.from_columns(cols, name='EBOUNDS')
    hdu.header.update({'EXTNAME': 'EBOUNDS', 'TELESCOP': 'COSI', 'INSTRUME': 'BTO', 'HDUCLASS': 'OGIP', 'HDUCLAS1': 'RESPONSE', 'HDUCLAS2': 'EBOUNDS', 'CHANTYPE': 'PI', 'DETCHANS': MAX_CHANNELS})
    return hdu

def inject_metadata(hdu, utc_start, utc_stop, is_primary=False):
    obs_id = f"{utc_start.strftime('%y%m%d')}000t"
    header_data = {'TELESCOP': ('COSI', 'Telescope'), 'INSTRUME': ('BTO', 'Instrument'), 'OBS_ID': (obs_id, 'Observation ID'), 'DATE-OBS': (utc_start.strftime('%Y-%m-%dT%H:%M:%S'), 'Start'), 'DATE-END': (utc_stop.strftime('%Y-%m-%dT%H:%M:%S'), 'End')}
    if is_primary: 
        header_data.update({'ORIGIN': ('UCB/SSL', 'Origin'), 'DATE': (datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%S'), 'Created'), 'CREATOR': ('BTO_LIVE_V1', 'Software')})
    else: 
        header_data.update({'HDUCLASS': ('OGIP', 'Standard'), 'DATAMODE': ('NORMAL', 'Datamode'), 'OBSERVER': ('BTO_TEAM', 'PI'), 'OBJECT': ('CAL_SOURCE', 'Target'), 'TIMESYS': ('UNIX', 'Time System'), 'MJDREFI': (40587, 'MJD Ref'), 'MJDREFF': (0.0, 'MJD offset')})
    for key, (val, comment) in header_data.items(): hdu.header[key] = (val, comment)

def flush_cache_to_disk(cache, tier):
    for path, data in list(cache.items()):
        if not data['rows']: continue
        try:
            cols = []
            for k in data['rows'][0].keys():
                fmt = '1D' if k=='TIME' else '1E' if k in ['t_ext', 't_det1', 'RATES'] else '29J' if k=='bins' else '1I' if k in ['PI', 'PHA', 'PKT_CNT', 't_ext_raw', 't_det1_raw'] else '1J'
                cols.append(fits.Column(name=k, format=fmt, array=np.array([r[k] for r in data['rows']])))
            
            new_hdu = fits.BinTableHDU.from_columns(cols, name='EVENTS' if 'PHA' in data['rows'][0] or 'PI' in data['rows'][0] else 'DATA')
            
            
            utc_start = datetime.datetime.fromtimestamp(data['rows'][0]['TIME'], tz=datetime.timezone.utc)
            utc_stop = datetime.datetime.fromtimestamp(data['rows'][-1]['TIME'], tz=datetime.timezone.utc)
            inject_metadata(new_hdu, utc_start, utc_stop)
            
            if os.path.exists(path):
                with fits.open(path, memmap=False) as hdul:
                    old_data = hdul[1].data
                    new_table = fits.BinTableHDU.from_columns(hdul[1].columns, nrows=len(old_data) + len(data['rows']))
                    for col in hdul[1].columns.names:
                        new_table.data[col] = np.concatenate([old_data[col], [r[col] for r in data['rows']]])
                    final_hdul = fits.HDUList([hdul[0], new_table])
                    if len(hdul) > 2: final_hdul.append(hdul[2])
                    final_hdul.writeto(path + ".tmp", overwrite=True)
                os.replace(path + ".tmp", path)
            else:
                hdul = [fits.PrimaryHDU(), new_hdu]
                if 'PI' in data['rows'][0]: hdul.append(get_ebounds_hdu())
                fits.HDUList(hdul).writeto(path, overwrite=True)
            print(f"[DISK IO] {tier} Flushed {len(data['rows'])} events to {os.path.basename(path)}")
        except Exception as e: print(f"[ERROR] {e}")
        del cache[path]

def read_input_stream(target_path, is_live):
    if is_live:
        while True:
            files = sorted(glob.glob(os.path.join(target_path, "*.txt")), key=os.path.getmtime)
            if not files: time.sleep(1); continue
            with open(files[-1], "r") as f:
                f.seek(0, 2)
                while True:
                    line = f.readline()
                    if not line: break
                    yield line
    else:
        files = [target_path] if os.path.isfile(target_path) else sorted(glob.glob(os.path.join(target_path, "*.txt")))
        for file in files:
            with open(file, 'r') as f:
                for line in f: yield line

def run_pipeline(input_path: str, levels: list, is_live: bool):
    router = ArchiveRouter(); cache_a, cache_b = {}, {}
    d7_timing_state = {'pps_time': 0, 'dwt_at_last_pps': 0, 'current_tid': 0}
    packet_count = 0
    try:
        for line in read_input_stream(input_path, is_live):
            if not line.strip() or line.startswith("#"): continue
            try:
                stream = bytes.fromhex(re.sub(r"[^0-9A-Fa-f]", "", line))
                idx = 0
                while idx <= len(stream) - 10:
                    if stream[idx:idx+2] != b'\xeb\x90': idx += 1; continue
                    apid = (int.from_bytes(stream[idx+2:idx+4], "big") & 0x07FF)
                    blen = 512 if apid in [0xD6, 0xD7] else 2 + (int.from_bytes(stream[idx+6:idx+8], "big") + 7) + 2
                    p = decode_packet(stream[idx+2:idx+blen-2], apid, d7_timing_state)
                    if p:
                        time_str = p['utc'].strftime('%Y-%m-%d %H:%M:%S')
                        print(f"[{time_str}] APID 0x{p['apid']:03X} ({p['type']:<3}) | Seq: {p['pkt_count']:05d} | Rows: {len(p['l1a'])}")
                        
                        if 'l0' in levels:
                            with open(router.get_path("L0", apid, p['utc']), "ab") as f0: f0.write(stream[idx:idx+blen])
                        if 'l1a' in levels:
                            if p['type'] == 'EVT':
                                for row in p['l1a']:
                                    pa = router.get_path("L1a", apid, p['utc'], row.pop('_TID', 0))
                                    if pa not in cache_a: cache_a[pa] = {'rows': [], 'met_list': []}
                                    cache_a[pa]['rows'].append(row)
                                    cache_a[pa]['met_list'].append(p['met'])
                            else:
                                pa = router.get_path("L1a", apid, p['utc'], p.get('tid',0))
                                if pa not in cache_a: cache_a[pa] = {'rows': [], 'met_list': []}
                                cache_a[pa]['rows'].extend(p['l1a'])
                                cache_a[pa]['met_list'].append(p['met'])
                        if 'l1b' in levels:
                            if p['type'] == 'EVT':
                                for row in p['l1b']:
                                    pb = router.get_path("L1b", apid, p['utc'], row.pop('_TID', 0))
                                    if pb not in cache_b: cache_b[pb] = {'rows': [], 'met_list': []}
                                    cache_b[pb]['rows'].append(row)
                                    cache_b[pb]['met_list'].append(p['met'])
                            else:
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
        if 'l1a' in levels: flush_cache_to_disk(cache_a, "L1a")
        if 'l1b' in levels: flush_cache_to_disk(cache_b, "L1b")
    except KeyboardInterrupt: pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--level", default="all", choices=['l0', 'l1a', 'l1b', 'all'])
    parser.add_argument("--live", choices=['yes', 'no'], default='no')
    args = parser.parse_args()
    run_pipeline(args.input, ['l0', 'l1a', 'l1b'] if args.level == 'all' else [args.level], args.live == 'yes')