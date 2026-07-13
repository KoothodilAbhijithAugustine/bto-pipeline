#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BTO PIPELINE
Processes CCSDS hex logs into L0, L1a, and L1b FITS archives.

"""

import argparse
import datetime
import glob
import os
import re
import time

import numpy as np
from astropy.io import fits

# =============================================================================
# 1. MISSION CONSTANTS & OFFSETS
# =============================================================================
GPS_EPOCH = datetime.datetime(1980, 1, 6, 0, 0, 0, tzinfo=datetime.timezone.utc)
MET_EPOCH = datetime.datetime(2025, 1, 1, 0, 0, 0, tzinfo=datetime.timezone.utc)

GPS_UTC_LEAP_SECONDS = [
    (datetime.datetime(1981, 7, 1, tzinfo=datetime.timezone.utc), 1),
    (datetime.datetime(1982, 7, 1, tzinfo=datetime.timezone.utc), 2),
    (datetime.datetime(1983, 7, 1, tzinfo=datetime.timezone.utc), 3),
    (datetime.datetime(1985, 7, 1, tzinfo=datetime.timezone.utc), 4),
    (datetime.datetime(1988, 1, 1, tzinfo=datetime.timezone.utc), 5),
    (datetime.datetime(1990, 1, 1, tzinfo=datetime.timezone.utc), 6),
    (datetime.datetime(1991, 1, 1, tzinfo=datetime.timezone.utc), 7),
    (datetime.datetime(1992, 7, 1, tzinfo=datetime.timezone.utc), 8),
    (datetime.datetime(1993, 7, 1, tzinfo=datetime.timezone.utc), 9),
    (datetime.datetime(1994, 7, 1, tzinfo=datetime.timezone.utc), 10),
    (datetime.datetime(1996, 1, 1, tzinfo=datetime.timezone.utc), 11),
    (datetime.datetime(1997, 7, 1, tzinfo=datetime.timezone.utc), 12),
    (datetime.datetime(1999, 1, 1, tzinfo=datetime.timezone.utc), 13),
    (datetime.datetime(2006, 1, 1, tzinfo=datetime.timezone.utc), 14),
    (datetime.datetime(2009, 1, 1, tzinfo=datetime.timezone.utc), 15),
    (datetime.datetime(2012, 7, 1, tzinfo=datetime.timezone.utc), 16),
    (datetime.datetime(2015, 7, 1, tzinfo=datetime.timezone.utc), 17),
    (datetime.datetime(2017, 1, 1, tzinfo=datetime.timezone.utc), 18),
]

DWT_24BIT_TICK_SEC = 853.33 / 1_000_000_000.0

CAL_HK = {
    "EXT": (0.1022, -275.66),
    "DET1": (0.1022, -275.66),
    "DET2": (0.1022, -275.66),
    "P12": (0.0037, 0.0),
    "M12": (0.0037, 0.0),
    "IMON5V": (0.0001665, 0.0),
}

E_SLOPE, E_INTERCEPT = 0.84098, -1.64736
MAX_CHANNELS = 4096

# CCSDS primary header offsets after 2-byte sync
PKT_ID_OFS, PKT_ID_LEN = 2, 2
PKT_SEQ_OFS, PKT_SEQ_LEN = 4, 2
PKT_LEN_OFS, PKT_LEN_LEN = 6, 2
PKT_SEC_OFS, PKT_SEC_LEN = 8, 4
PKT_TKS_OFS, PKT_TKS_LEN = 12, 4

# Histogram / lightcurve packet layout
LC_BTO_ID_OFS = 16
LC_BINS_START = 17
LC_BINS_STEP = 2
LC_NUM_BINS = 29

# Event packet layout
EVT_DATA_START = 22
EVT_WORD_LEN = 8

# Housekeeping packet layout
HK_BTO_ID_OFS = 16
HK_MODE_OFS = 17
HK_SUT_OFS = 18
HK_CMD_CNT_OFS = 22
HK_FAIL_CMD_CNT_OFS = 24
HK_CS_ERR_CNT_OFS = 26
HK_FSW_VER_OFS = 28
HK_WD_RESET_OFS = 29
HK_RESET_REASON_OFS = 30
HK_FAULT_STATUS_OFS = 31
HK_DET_PWR_STATUS_OFS = 32
HK_ANALOG_STATUS_OFS = 33
HK_HIST_WPTR_OFS = 34
HK_HIST_RPTR_OFS = 36
HK_HIST_NAND_WPTR_OFS = 38
HK_HIST_NAND_RPTR_OFS = 42
HK_PHOT_NAND_WPTR_OFS = 46
HK_PHOT_NAND_RPTR_OFS = 50
HK_FLASH_ERASE_FAIL_OFS = 54
HK_LAST_GRB_TIME_OFS = 56
HK_LAST_GRB_NAND_OFS = 60
HK_LAST_GRB_ID_OFS = 64
HK_ZC_CNT_OFS = 68
HK_SU_CNT_OFS = 70
HK_TEMP_EXT_OFS = 72
HK_TEMP_DET1_OFS = 74
HK_TEMP_DET2_OFS = 76
HK_VMON_P12_1_OFS = 78
HK_VMON_M12_1_OFS = 80
HK_VMON_P12_2_OFS = 82
HK_VMON_M12_2_OFS = 84
HK_IMON_5V_1_OFS = 86
HK_IMON_5V_2_OFS = 88
HK_SPARE_OFS = 90
HK_IT_CS_OFS = 92

FLUSH_THRESHOLD = 1000
MIN_VALID_UTC = datetime.datetime(2024, 1, 1, tzinfo=datetime.timezone.utc)
MAX_VALID_UTC = datetime.datetime(2030, 1, 1, tzinfo=datetime.timezone.utc)

# =============================================================================
# 2. TIME & ROUTING UTILITIES
# =============================================================================
def gps_utc_offset_for_utc(dt_utc: datetime.datetime) -> int:
    offset = 0
    for effective_utc, leap_count in GPS_UTC_LEAP_SECONDS:
        if dt_utc >= effective_utc:
            offset = leap_count
        else:
            break
    return offset

def gps_seconds_to_utc(gps_sec: int, frac_sec: float = 0.0) -> datetime.datetime:
    gps_dt = GPS_EPOCH + datetime.timedelta(seconds=float(gps_sec) + frac_sec)
    utc_dt = gps_dt
    for _ in range(3):
        offset = gps_utc_offset_for_utc(utc_dt)
        new_utc_dt = gps_dt - datetime.timedelta(seconds=offset)
        if new_utc_dt == utc_dt:
            break
        utc_dt = new_utc_dt
    return utc_dt

def get_met(dt_utc: datetime.datetime) -> float:
    return (dt_utc - MET_EPOCH).total_seconds()

def clean_hex_fragment(text: str) -> str:
    text = text.split("#", 1)[0]
    return re.sub(r"[^0-9A-Fa-f]", "", text)

class ArchiveRouter:
    def __init__(self, root="BTO_Data_Archive"):
        self.root = root

    def get_path(self, tier: str, apid: int, dt: datetime.datetime, tid: int = 0) -> str:
        yyyy, mm, dd = dt.strftime("%Y"), dt.strftime("%m"), dt.strftime("%d")
        yymmdd = dt.strftime("%y%m%d")

        if tier == "L0":
            folder = {
                0xD6: "D6_Spectrum",
                0xD7: "D7_Event_packages",
                0xD8: "D8_House_keeping",
            }.get(apid, "Other")
            base = os.path.join(self.root, "L0_binary", yyyy, mm, dd, folder)
            fname = f"bto_{yymmdd}_{folder.lower()}.bin"
        else:
            base = os.path.join(self.root, f"{tier}_fits", yyyy, mm, dd)
            if apid == 0xD6:
                fname = f"cs{yymmdd}bto_lc.fits"
            elif apid == 0xD8:
                fname = f"cs{yymmdd}bto_hk.fits"
            elif apid == 0xD7:
                fname = f"cs{yymmdd}_{tid:05d}_bto_evt.fits"
            else:
                fname = f"cs{yymmdd}_misc.fits"

        os.makedirs(base, exist_ok=True)
        return os.path.join(base, fname)

# =============================================================================
# 3. PACKET DECODER
# =============================================================================
def decode_packet(ccsds: bytes, apid: int, d7_state: dict):
    if len(ccsds) < 16:
        return None

    pkt_count = int.from_bytes(ccsds[PKT_SEQ_OFS:PKT_SEQ_OFS + PKT_SEQ_LEN], "big") & 0x3FFF
    sec = int.from_bytes(ccsds[PKT_SEC_OFS:PKT_SEC_OFS + PKT_SEC_LEN], "big")
    tks = int.from_bytes(ccsds[PKT_TKS_OFS:PKT_TKS_OFS + PKT_TKS_LEN], "big")

    packet_frac = float(tks) / 300_000_000.0
    utc_dt = gps_seconds_to_utc(sec, packet_frac)

    if utc_dt < MIN_VALID_UTC or utc_dt > MAX_VALID_UTC:
        return None

    packet_met = get_met(utc_dt)
    packet_dwt_24bit = (tks >> 8) & 0xFFFFFF
    
    base_data = {"apid": apid, "pkt_count": pkt_count, "utc": utc_dt, "met": packet_met}

    # ---------------------------------------------------------------------
    # Housekeeping
    # ---------------------------------------------------------------------
    if apid == 0x0D8:
        if len(ccsds) < (HK_IT_CS_OFS + 2):
            return None

        bto_id = ccsds[HK_BTO_ID_OFS]
        bto_mode = ccsds[HK_MODE_OFS]
        sut_time = int.from_bytes(ccsds[HK_SUT_OFS:HK_SUT_OFS + 4], "big")
        
        d7_state["hk_sut"] = sut_time
        
        cmd_cnt = int.from_bytes(ccsds[HK_CMD_CNT_OFS:HK_CMD_CNT_OFS + 2], "big")
        fail_cmd_cnt = int.from_bytes(ccsds[HK_FAIL_CMD_CNT_OFS:HK_FAIL_CMD_CNT_OFS + 2], "big")
        cs_err_cnt = int.from_bytes(ccsds[HK_CS_ERR_CNT_OFS:HK_CS_ERR_CNT_OFS + 2], "big")
        fsw_ver = ccsds[HK_FSW_VER_OFS]
        wd_reset = ccsds[HK_WD_RESET_OFS]
        reset_reason = ccsds[HK_RESET_REASON_OFS]
        fault_status = ccsds[HK_FAULT_STATUS_OFS]
        det_pwr_status = ccsds[HK_DET_PWR_STATUS_OFS]
        analog_status = ccsds[HK_ANALOG_STATUS_OFS]

        hist_wptr = int.from_bytes(ccsds[HK_HIST_WPTR_OFS:HK_HIST_WPTR_OFS + 2], "big")
        hist_rptr = int.from_bytes(ccsds[HK_HIST_RPTR_OFS:HK_HIST_RPTR_OFS + 2], "big")
        hist_nand_wptr = int.from_bytes(ccsds[HK_HIST_NAND_WPTR_OFS:HK_HIST_NAND_WPTR_OFS + 4], "big")
        hist_nand_rptr = int.from_bytes(ccsds[HK_HIST_NAND_RPTR_OFS:HK_HIST_NAND_RPTR_OFS + 4], "big")
        phot_nand_wptr = int.from_bytes(ccsds[HK_PHOT_NAND_WPTR_OFS:HK_PHOT_NAND_WPTR_OFS + 4], "big")
        phot_nand_rptr = int.from_bytes(ccsds[HK_PHOT_NAND_RPTR_OFS:HK_PHOT_NAND_RPTR_OFS + 4], "big")
        flash_erase_fail = int.from_bytes(ccsds[HK_FLASH_ERASE_FAIL_OFS:HK_FLASH_ERASE_FAIL_OFS + 2], "big")
        last_grb_time = int.from_bytes(ccsds[HK_LAST_GRB_TIME_OFS:HK_LAST_GRB_TIME_OFS + 4], "big")
        last_grb_nand = int.from_bytes(ccsds[HK_LAST_GRB_NAND_OFS:HK_LAST_GRB_NAND_OFS + 4], "big")
        last_grb_id = int.from_bytes(ccsds[HK_LAST_GRB_ID_OFS:HK_LAST_GRB_ID_OFS + 4], "big")
        zc_cnt = int.from_bytes(ccsds[HK_ZC_CNT_OFS:HK_ZC_CNT_OFS + 2], "big")
        su_cnt = int.from_bytes(ccsds[HK_SU_CNT_OFS:HK_SU_CNT_OFS + 2], "big")

        t_ext_raw = int.from_bytes(ccsds[HK_TEMP_EXT_OFS:HK_TEMP_EXT_OFS + 2], "big")
        t_det1_raw = int.from_bytes(ccsds[HK_TEMP_DET1_OFS:HK_TEMP_DET1_OFS + 2], "big")
        t_det2_raw = int.from_bytes(ccsds[HK_TEMP_DET2_OFS:HK_TEMP_DET2_OFS + 2], "big")
        vmon_p12_1_raw = int.from_bytes(ccsds[HK_VMON_P12_1_OFS:HK_VMON_P12_1_OFS + 2], "big")
        vmon_m12_1_raw = int.from_bytes(ccsds[HK_VMON_M12_1_OFS:HK_VMON_M12_1_OFS + 2], "big")
        vmon_p12_2_raw = int.from_bytes(ccsds[HK_VMON_P12_2_OFS:HK_VMON_P12_2_OFS + 2], "big")
        vmon_m12_2_raw = int.from_bytes(ccsds[HK_VMON_M12_2_OFS:HK_VMON_M12_2_OFS + 2], "big")
        imon_5v_1_raw = int.from_bytes(ccsds[HK_IMON_5V_1_OFS:HK_IMON_5V_1_OFS + 2], "big")
        imon_5v_2_raw = int.from_bytes(ccsds[HK_IMON_5V_2_OFS:HK_IMON_5V_2_OFS + 2], "big")
        spare_raw = int.from_bytes(ccsds[HK_SPARE_OFS:HK_SPARE_OFS + 2], "big")
        it_cs_raw = int.from_bytes(ccsds[HK_IT_CS_OFS:HK_IT_CS_OFS + 2], "big")

        base_data.update({
            "type": "HK",
            "l1a": [{
                "TIME": packet_met,
                "PKT_CNT": pkt_count,
                "BTO_ID": bto_id,
                "BTO_MODE": bto_mode,
                "SUT": sut_time,
                "CMD_CNT": cmd_cnt,
                "FAIL_CMD_CNT": fail_cmd_cnt,
                "CS_ERR_CNT": cs_err_cnt,
                "FSW_VER": fsw_ver,
                "WD_RESET": wd_reset,
                "RESET_REASON": reset_reason,
                "FAULT_STATUS": fault_status,
                "DET_PWR_STATUS": det_pwr_status,
                "ANALOG_STATUS": analog_status,
                "HIST_WPTR": hist_wptr,
                "HIST_RPTR": hist_rptr,
                "HIST_NAND_WPTR": hist_nand_wptr,
                "HIST_NAND_RPTR": hist_nand_rptr,
                "PHOT_NAND_WPTR": phot_nand_wptr,
                "PHOT_NAND_RPTR": phot_nand_rptr,
                "FLASH_ERASE_FAIL": flash_erase_fail,
                "LAST_GRB_TIME": last_grb_time,
                "LAST_GRB_NAND": last_grb_nand,
                "LAST_GRB_ID": last_grb_id,
                "ZC_CNT": zc_cnt,
                "SU_CNT": su_cnt,
                "t_ext_raw": t_ext_raw,
                "t_det1_raw": t_det1_raw,
                "t_det2_raw": t_det2_raw,
                "vmon_p12_1_raw": vmon_p12_1_raw,
                "vmon_m12_1_raw": vmon_m12_1_raw,
                "vmon_p12_2_raw": vmon_p12_2_raw,
                "vmon_m12_2_raw": vmon_m12_2_raw,
                "imon_5v_1_raw": imon_5v_1_raw,
                "imon_5v_2_raw": imon_5v_2_raw,
                "spare_raw": spare_raw,
                "it_cs_raw": it_cs_raw,
            }],
            "l1b": [{
                "TIME": packet_met,
                "PKT_CNT": pkt_count,
                "BTO_ID": bto_id,
                "BTO_MODE": bto_mode,
                "SUT": sut_time,
                "CMD_CNT": cmd_cnt,
                "FAIL_CMD_CNT": fail_cmd_cnt,
                "CS_ERR_CNT": cs_err_cnt,
                "FSW_VER": fsw_ver,
                "WD_RESET": wd_reset,
                "RESET_REASON": reset_reason,
                "FAULT_STATUS": fault_status,
                "DET_PWR_STATUS": det_pwr_status,
                "ANALOG_STATUS": analog_status,
                "HIST_WPTR": hist_wptr,
                "HIST_RPTR": hist_rptr,
                "HIST_NAND_WPTR": hist_nand_wptr,
                "HIST_NAND_RPTR": hist_nand_rptr,
                "PHOT_NAND_WPTR": phot_nand_wptr,
                "PHOT_NAND_RPTR": phot_nand_rptr,
                "FLASH_ERASE_FAIL": flash_erase_fail,
                "LAST_GRB_TIME": last_grb_time,
                "LAST_GRB_NAND": last_grb_nand,
                "LAST_GRB_ID": last_grb_id,
                "ZC_CNT": zc_cnt,
                "SU_CNT": su_cnt,
                "t_ext_raw": t_ext_raw,
                "t_det1_raw": t_det1_raw,
                "t_det2_raw": t_det2_raw,
                "vmon_p12_1_raw": vmon_p12_1_raw,
                "vmon_m12_1_raw": vmon_m12_1_raw,
                "vmon_p12_2_raw": vmon_p12_2_raw,
                "vmon_m12_2_raw": vmon_m12_2_raw,
                "imon_5v_1_raw": imon_5v_1_raw,
                "imon_5v_2_raw": imon_5v_2_raw,
                "spare_raw": spare_raw,
                "it_cs_raw": it_cs_raw,
                "t_ext": (t_ext_raw * CAL_HK["EXT"][0]) + CAL_HK["EXT"][1],
                "t_det1": (t_det1_raw * CAL_HK["DET1"][0]) + CAL_HK["DET1"][1],
                "t_det2": (t_det2_raw * CAL_HK["DET2"][0]) + CAL_HK["DET2"][1],
            }],
        })
        return base_data

    # ---------------------------------------------------------------------
    # Lightcurve / histogram
    # ---------------------------------------------------------------------
    elif apid == 0x0D6:
        if len(ccsds) < (LC_BINS_START + (LC_NUM_BINS * LC_BINS_STEP)):
            return None

        raw_bins = [
            int.from_bytes(
                ccsds[LC_BINS_START + (i * LC_BINS_STEP): LC_BINS_START + (i * LC_BINS_STEP) + 2],
                "big",
            )
            for i in range(LC_NUM_BINS)
        ]

        base_data.update({
            "type": "LC",
            "l1a": [{"TIME": packet_met, "PKT_CNT": pkt_count, "bins": raw_bins}],
            "l1b": [{"TIME": packet_met, "PKT_CNT": pkt_count, "COUNTS": raw_bins}],
        })
        return base_data

    # ---------------------------------------------------------------------
    # Event-by-event
    # ---------------------------------------------------------------------
    elif apid == 0x0D7:
        l1a, l1b = [], []
        ptr = EVT_DATA_START

        while ptr <= len(ccsds) - EVT_WORD_LEN:
            wb = ccsds[ptr:ptr + EVT_WORD_LEN]

            if wb == b"\x00\x00\x00\x00\x00\x00\x00\x00":
                ptr += EVT_WORD_LEN
                continue

            adc_data = int.from_bytes(wb[0:2], "big")
            ts_long = int.from_bytes(wb[2:4], "big")
            ts_short = int.from_bytes(wb[4:6], "big")
            raw_deadtime = int.from_bytes(wb[6:8], "big")

            # -----------------------------------------------
            # 0. Skip Multi-Word Event Blocks (Cross-Packet Safe)
            # -----------------------------------------------
            if d7_state.get("skip_next", 0) > 0:
                d7_state["skip_next"] -= 1
                ptr += EVT_WORD_LEN
                continue

            # -----------------------------------------------
            # 1. Empty NAND Padding / Heartbeat Evaluator
            # -----------------------------------------------
            if adc_data == 0xFFFF and ts_long == 0xFFFF:
                if ts_short == 0xFFFF:
                    # Pure empty NAND flash padding. Discard.
                    pass
                else:
                    # True Heartbeat (Event 1: Seq num). 
                    # The next 2 events in memory are PPS Time and SUT Time. 
                    # We flag the state machine to explicitly skip them.
                    d7_state["skip_next"] = 2
                ptr += EVT_WORD_LEN
                continue

            # -----------------------------------------------
            # 2. Start Frame Marker
            # -----------------------------------------------
            if adc_data == 0x0123 and ts_long == 0x4567:
                d7_state["current_tid"] += 1
                
                pps_time = (raw_deadtime << 16) | ts_short
                if pps_time >= 0:
                    if pps_time < 1_000_000_000:
                        anchor_utc = MET_EPOCH + datetime.timedelta(seconds=pps_time)
                    else:
                        anchor_utc = gps_seconds_to_utc(pps_time)
                        
                    if MIN_VALID_UTC <= anchor_utc <= MAX_VALID_UTC:
                        d7_state["anchor_met"] = get_met(anchor_utc)
                        d7_state["last_evt_met"] = d7_state["anchor_met"]
                
                d7_state["last_photon_dwt"] = None
                ptr += EVT_WORD_LEN
                continue

            # -----------------------------------------------
            # 3. End Frame Marker
            # -----------------------------------------------
            if adc_data == 0xABCD:
                if ts_long == 0xABCD:
                    # End Frame 1 (PPS Time). Skip the next 2 events.
                    d7_state["skip_next"] = 2
                elif ts_long in [0xEF00, 0xEFFF]:
                    # Failsafe: End Frame 2 (PPS Subsec) for Trigger or Calibration mode.
                    d7_state["skip_next"] = 1
                ptr += EVT_WORD_LEN
                continue

            # -----------------------------------------------
            # 4. Photon Event (Standard & Saturated)
            # -----------------------------------------------
            
            # Extract Hardware Flags
            flag_su = (adc_data >> 15) & 1
            flag_up = (adc_data >> 14) & 1
            flag_ld = (adc_data >> 13) & 1
            flag_pseudo = (adc_data >> 12) & 1

            # 12-Bit Two's Complement Conversion
            if (adc_data & 0x7FFF) == 0x7FFF:
                # Flight software explicit overflow/saturated marker.
                # Force it to the maximum FITS channel to prevent rollover.
                pha = 4095 
            else:
                pha_raw = adc_data & 0x0FFF
                pha_signed = pha_raw - 4096 if pha_raw >= 2048 else pha_raw
                
                # FITS Format Correction:
                # Taking the absolute value rectifies negative signal generator 
                # pulses (-51mV, etc) into positive FITS channels and floors baseline noise.
                pha = abs(pha_signed)
            
            full_timestamp = (ts_long << 16) | ts_short
            photon_dwt_24bit = (full_timestamp >> 8) & 0xFFFFFF

            if d7_state.get("last_photon_dwt") is None:
                d7_state["last_photon_dwt"] = photon_dwt_24bit
                if d7_state.get("last_evt_met") is None:
                    d7_state["last_evt_met"] = packet_met

            tick_diff = photon_dwt_24bit - d7_state["last_photon_dwt"]
            
            if tick_diff < -0x800000: 
                tick_diff += 0x1000000 
            elif tick_diff > 0x800000: 
                tick_diff -= 0x1000000

            evt_met = d7_state["last_evt_met"] + (tick_diff * DWT_24BIT_TICK_SEC)
            
            d7_state["last_photon_dwt"] = photon_dwt_24bit
            d7_state["last_evt_met"] = evt_met

            active_tid = max(1, d7_state.get("current_tid", 1))
            
            l1a.append({
                "TIME": evt_met,
                "PKT_CNT": pkt_count,
                "PHA": pha,
                "FLAG_SU": flag_su,
                "FLAG_UP": flag_up,
                "FLAG_LD": flag_ld,
                "FLAG_PSEUDO": flag_pseudo,
                "DEADTIME": raw_deadtime,
                "_TID": active_tid,
            })
            l1b.append({
                "TIME": evt_met,
                "PKT_CNT": pkt_count,
                "PI": pha,
                "_TID": active_tid,
            })

            ptr += EVT_WORD_LEN

        base_data.update({"type": "EVT", "l1a": l1a, "l1b": l1b})
        return base_data

    return None

# =============================================================================
# 4. HEADER & FITS UTILITIES
# =============================================================================
def get_ebounds_hdu():
    channels = np.arange(MAX_CHANNELS, dtype=np.int16)
    e_min = np.maximum(0, E_SLOPE * (channels - 0.5) + E_INTERCEPT)
    e_max = np.maximum(0, E_SLOPE * (channels + 0.5) + E_INTERCEPT)

    cols = [
        fits.Column(name="CHANNEL", format="1I", array=channels),
        fits.Column(name="E_MIN", format="1E", array=e_min, unit="keV"),
        fits.Column(name="E_MAX", format="1E", array=e_max, unit="keV"),
    ]
    hdu = fits.BinTableHDU.from_columns(cols, name="EBOUNDS")
    hdu.header.update({
        "EXTNAME": "EBOUNDS",
        "TELESCOP": "COSI",
        "INSTRUME": "BTO",
        "HDUCLASS": "OGIP",
        "HDUCLAS1": "RESPONSE",
        "HDUCLAS2": "EBOUNDS",
        "CHANTYPE": "PI",
        "DETCHANS": MAX_CHANNELS,
    })
    return hdu

def get_eneband_hdu():
    bin_width = MAX_CHANNELS // LC_NUM_BINS
    minchan = np.arange(LC_NUM_BINS, dtype=np.int16) * bin_width
    maxchan = minchan + bin_width - 1
    maxchan[-1] = MAX_CHANNELS - 1

    e_min = np.maximum(0, E_SLOPE * (minchan - 0.5) + E_INTERCEPT)
    e_max = np.maximum(0, E_SLOPE * (maxchan + 0.5) + E_INTERCEPT)

    cols = [
        fits.Column(name="MINCHAN", format="1I", array=minchan),
        fits.Column(name="MAXCHAN", format="1I", array=maxchan),
        fits.Column(name="E_MIN", format="1E", array=e_min.astype(np.float32), unit="keV"),
        fits.Column(name="E_MAX", format="1E", array=e_max.astype(np.float32), unit="keV"),
    ]
    hdu = fits.BinTableHDU.from_columns(cols, name="ENEBAND")
    hdu.header.update({
        "EXTNAME": "ENEBAND",
        "TELESCOP": "COSI",
        "INSTRUME": "BTO",
        "HDUCLASS": "OGIP",
        "HDUCLAS1": "RESPONSE",
        "HDUCLAS2": "EBOUNDS",
    })
    return hdu

def make_obs_id(utc_start: datetime.datetime, is_event: bool) -> str:
    if is_event:
        return f"{utc_start.strftime('%y%m%d')}000t"
    return utc_start.strftime("%Y%m%d")

def inject_metadata(hdu, t_start, t_stop, utc_start, utc_stop, is_primary=False, obs_id=None):
    if obs_id is None:
        obs_id = utc_start.strftime("%y%m%d") + "000t"

    header_data = {
        "TELESCOP": ("COSI", "Telescope"),
        "INSTRUME": ("BTO", "Instrument"),
        "OBS_ID": (obs_id, "Observation ID"),
        "DATE-OBS": (utc_start.strftime("%Y-%m-%dT%H:%M:%S"), "Start"),
        "DATE-END": (utc_stop.strftime("%Y-%m-%dT%H:%M:%S"), "End"),
    }

    if is_primary:
        header_data.update({
            "ORIGIN": ("UCB/SSL", "Origin"),
            "DATE": (datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S"), "Created"),
            "CREATOR": ("BTO_LIVE_V5.23", "Software"),
        })
    else:
        header_data.update({
            "HDUCLASS": ("OGIP", "Standard"),
            "DATAMODE": ("NORMAL", "Datamode"),
            "OBSERVER": ("BTO_TEAM", "PI"),
            "OBJECT": ("CAL_SOURCE", "Target"),
            "TIMESYS": ("TT", "Time System"),
            "MJDREFI": (60676, "MJD Ref"),
            "MJDREFF": (0.0008007407407407, "MJD offset"),
            "TIMEREF": ("LOCAL", "Ref Frame"),
            "TASSIGN": ("SATELLITE", "Time clock"),
            "TIMEUNIT": ("s", "Time unit"),
            "TSTART": (t_start, "Start MET"),
            "TSTOP": (t_stop, "Stop MET"),
            "CLOCKAPP": ("F", "Clock corr?"),
        })

    for key, (val, comment) in header_data.items():
        hdu.header[key] = (val, comment)

def _make_column(name, values):
    arr = np.array(values)

    if name == "TIME":
        fmt = "1D"
    elif name in ["t_ext", "t_det1", "t_det2"]:
        fmt = "1E"
    elif name in ["bins", "COUNTS"]:
        fmt = f"{LC_NUM_BINS}J"
    elif name in ["PKT_CNT", "BTO_ID", "BTO_MODE", "FSW_VER", "WD_RESET", "RESET_REASON", "FAULT_STATUS", "DET_PWR_STATUS", "ANALOG_STATUS", "FLAG_SU", "FLAG_UP", "FLAG_LD", "FLAG_PSEUDO"]:
        fmt = "1I"
    elif name in ["SUT", "CMD_CNT", "FAIL_CMD_CNT", "CS_ERR_CNT", "HIST_WPTR", "HIST_RPTR",
                  "FLASH_ERASE_FAIL", "ZC_CNT", "SU_CNT",
                  "t_ext_raw", "t_det1_raw", "t_det2_raw",
                  "vmon_p12_1_raw", "vmon_m12_1_raw", "vmon_p12_2_raw", "vmon_m12_2_raw",
                  "imon_5v_1_raw", "imon_5v_2_raw", "spare_raw", "it_cs_raw",
                  "HIST_NAND_WPTR", "HIST_NAND_RPTR", "PHOT_NAND_WPTR", "PHOT_NAND_RPTR",
                  "LAST_GRB_TIME", "LAST_GRB_NAND", "LAST_GRB_ID", "DEADTIME"]:
        fmt = "1J"
    elif name in ["PI", "PHA"]:
        fmt = "1I"  # FITS 1I is inherently a SIGNED 16-bit integer
    else:
        fmt = "1J"

    return fits.Column(name=name, format=fmt, array=arr)

def flush_cache_to_disk(cache, tier):
    for path, data in list(cache.items()):
        if not data["rows"]:
            continue

        try:
            rows = sorted(data["rows"], key=lambda r: float(r["TIME"]))
            first_row = rows[0]

            if "PI" in first_row or "PHA" in first_row:
                ext_name, hdu_class = "EVENTS", "EVENTS"
            elif "t_ext_raw" in first_row:
                ext_name, hdu_class = "HK", "HK"
            elif "bins" in first_row or "COUNTS" in first_row:
                ext_name, hdu_class = "RATE", "LIGHTCURVE"
            else:
                ext_name, hdu_class = "DATA", "DATA"

            cols = []
            for k in first_row.keys():
                values = [r[k] for r in rows]
                cols.append(_make_column(k, values))

            new_hdu = fits.BinTableHDU.from_columns(cols, name=ext_name)

            t_s = float(rows[0]["TIME"])
            t_e = float(rows[-1]["TIME"])
            utc_start = MET_EPOCH + datetime.timedelta(seconds=t_s)
            utc_stop = MET_EPOCH + datetime.timedelta(seconds=t_e)

            obs_id = utc_start.strftime("%y%m%d") + ("000t" if ext_name == "EVENTS" else "")
            inject_metadata(new_hdu, t_s, t_e, utc_start, utc_stop, obs_id=obs_id)
            new_hdu.header.update({"HDUCLAS1": hdu_class, "HDUCLAS2": "ALL"})

            if ext_name == "EVENTS" and tier == "L1b":
                new_hdu.header.update({"CHANTYPE": "PI", "DETCHANS": MAX_CHANNELS})

            if os.path.exists(path):
                with fits.open(path, memmap=False) as hdul:
                    old_data = hdul[1].data
                    new_tstop = max(hdul[1].header.get("TSTOP", t_e), t_e)

                    new_table = fits.BinTableHDU.from_columns(
                        hdul[1].columns,
                        nrows=len(old_data) + len(rows),
                        name=hdul[1].name,
                        header=hdul[1].header,
                    )

                    for col in hdul[1].columns.names:
                        incoming = np.array([r[col] for r in rows])
                        new_table.data[col] = np.concatenate([old_data[col], incoming])

                    new_table.header["TSTOP"] = new_tstop
                    final_hdul = fits.HDUList([hdul[0], new_table])
                    for ext in hdul[2:]:
                        final_hdul.append(ext.copy())

                    final_hdul[0].header["DATE-END"] = utc_stop.strftime("%Y-%m-%dT%H:%M:%S")
                    final_hdul.writeto(path + ".tmp", overwrite=True)

                os.replace(path + ".tmp", path)
            else:
                hdul = [fits.PrimaryHDU(), new_hdu]
                inject_metadata(hdul[0], t_s, t_e, utc_start, utc_stop, is_primary=True, obs_id=obs_id)

                if tier == "L1b":
                    if "PI" in first_row:
                        hdul.append(get_ebounds_hdu())
                    elif "COUNTS" in first_row:
                        hdul.append(get_eneband_hdu())

                fits.HDUList(hdul).writeto(path, overwrite=True)

            print(f"[DISK IO] {tier} Flushed {len(rows)} rows to {os.path.basename(path)}")
        except Exception as e:
            print(f"[ERROR] {e}")

        del cache[path]

# =============================================================================
# 5. EXECUTION
# =============================================================================
def read_input_stream(target_path, is_live):
    if is_live:
        while True:
            files = sorted(glob.glob(os.path.join(target_path, "*.txt")), key=os.path.getmtime)
            if not files:
                time.sleep(1)
                continue

            latest = files[-1]
            try:
                with open(latest, "r") as f:
                    f.seek(0, os.SEEK_END)
                    while True:
                        line = f.readline()
                        if not line:
                            break
                        clean = clean_hex_fragment(line)
                        if len(clean) % 2 != 0:
                            clean = clean[:-1]
                        if clean:
                            yield clean
            except FileNotFoundError:
                time.sleep(1)

            time.sleep(0.25)
    else:
        files = (
            [target_path]
            if os.path.isfile(target_path)
            else sorted(glob.glob(os.path.join(target_path, "*.txt")))
        )

        for file in files:
            hex_lines = []
            with open(file, "r") as f:
                for line in f:
                    clean = clean_hex_fragment(line.strip())
                    if clean:
                        hex_lines.append(clean)

            full_hex = "".join(hex_lines)

            if len(full_hex) % 2 != 0:
                print(f"[WARNING] File {file} ended with incomplete byte. Truncating last char.")
                full_hex = full_hex[:-1]

            if full_hex:
                yield full_hex

def run_pipeline(input_path: str, levels: list, is_live: bool):
    router = ArchiveRouter()
    cache_a, cache_b = {}, {}
    
    global_d7_state = {
        "current_tid": 0, 
        "last_sut_val": 0,
        "last_photon_dwt": None,
        "last_evt_met": None,
        "hk_sut": None,
        "skip_next": 0   # Track multi-word events across packet boundaries
    }
     
    packet_count = 0
    total_parsed = 0
    
    for hex_chunk in read_input_stream(input_path, is_live):
        if not hex_chunk:
            continue

        try:
            stream = bytes.fromhex(hex_chunk)
        except ValueError as e:
            print(f"[WARNING] Skipping invalid hex chunk: {e}")
            continue

        idx = 0
        while idx <= len(stream) - 16:
            if stream[idx:idx + 2] != b"\xeb\x90":
                idx += 1
                continue

            apid = int.from_bytes(stream[idx + 2:idx + 4], "big") & 0x07FF

            if apid not in [0x0D6, 0x0D7, 0x0D8]:
                idx += 1
                continue

            pkt_data_len = int.from_bytes(stream[idx + 6:idx + 8], "big")
            blen = pkt_data_len + 9

            if blen > 65536 or idx + blen > len(stream):
                idx += 1
                continue

            packet = stream[idx:idx + blen]
            p = decode_packet(packet, apid, global_d7_state)

            if p:
                total_parsed += 1
                if "l0" in levels:
                    with open(router.get_path("L0", apid, p["utc"]), "ab") as f0:
                        f0.write(packet)

                if "l1a" in levels and p.get("type") in ["EVT", "HK", "LC"]:
                    if p["type"] == "EVT":
                        for row in p["l1a"]:
                            pa = router.get_path("L1a", apid, p["utc"], row.pop("_TID", 0))
                            if pa not in cache_a:
                                cache_a[pa] = {"rows": [], "met_list": []}
                            cache_a[pa]["rows"].append(row)
                            cache_a[pa]["met_list"].append(p["met"])
                    else:
                        pa = router.get_path("L1a", apid, p["utc"], p.get("tid", 0))
                        if pa not in cache_a:
                            cache_a[pa] = {"rows": [], "met_list": []}
                        cache_a[pa]["rows"].extend(p["l1a"])
                        cache_a[pa]["met_list"].append(p["met"])

                if "l1b" in levels and p.get("type") in ["EVT", "HK", "LC"]:
                    if p["type"] == "EVT":
                        for row in p["l1b"]:
                            pb = router.get_path("L1b", apid, p["utc"], row.pop("_TID", 0))
                            if pb not in cache_b:
                                cache_b[pb] = {"rows": [], "met_list": []}
                            cache_b[pb]["rows"].append(row)
                            cache_b[pb]["met_list"].append(p["met"])
                    else:
                        pb = router.get_path("L1b", apid, p["utc"], p.get("tid", 0))
                        if pb not in cache_b:
                            cache_b[pb] = {"rows": [], "met_list": []}
                        cache_b[pb]["rows"].extend(p["l1b"])
                        cache_b[pb]["met_list"].append(p["met"])

                packet_count += 1
                idx += blen
            else:
                idx += 2

        if packet_count >= FLUSH_THRESHOLD:
            if "l1a" in levels:
                flush_cache_to_disk(cache_a, "L1a")
            if "l1b" in levels:
                flush_cache_to_disk(cache_b, "L1b")
            packet_count = 0

    if "l1a" in levels:
        flush_cache_to_disk(cache_a, "L1a")
    if "l1b" in levels:
        flush_cache_to_disk(cache_b, "L1b")
        
    print(f"\n[PIPELINE COMPLETE] Valid Packets Parsed: {total_parsed}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--level", default="all", choices=["l0", "l1a", "l1b", "all"])
    parser.add_argument("--live", choices=["yes", "no"], default="no")
    args = parser.parse_args()

    run_pipeline(
        args.input,
        ["l0", "l1a", "l1b"] if args.level == "all" else [args.level],
        args.live == "yes",
    )
