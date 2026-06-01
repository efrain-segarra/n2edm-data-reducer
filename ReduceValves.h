#pragma once

#include <string>
#include <fstream>
#include <regex>
#include <cstdint>
#include <cstdio>
#include <iostream>

#include "constants.h"

// Binary format: 13 × uint64_t per row + uint64_t sentinel = 112 bytes/row
static constexpr int      VALVES_NCOLS    = 13;
static constexpr int      VALVES_ROW_BYTES = (VALVES_NCOLS + 1) * 8;
static constexpr uint64_t VALVES_SENTINEL  = 0xDEADBEEFDEADBEEFULL;

static constexpr int COL_TIMESTAMP = 0;
static constexpr int COL_DOUTA     = 1;
static constexpr int COL_DOUTB     = 2;

// Bit positions within DOUTA (A_n = bit n-1)
static constexpr int BIT_A19 = 18;  // Hg reservoir
static constexpr int BIT_A27 = 26;  // neutron shutter (N.O.)
static constexpr int BIT_A29 = 28;  // Hg shutter (N.O.)

struct ValvesResult {
    double T_fill_end      = DUMMY_VAL;  // A27 first closes  — filling stops
    double T_fill_hg_start = DUMMY_VAL;  // A29 first opens   — Hg released into chamber
    double T_fill_hg_stop  = DUMMY_VAL;  // A29 first closes  — Hg shutter closed
    double T_count_start   = DUMMY_VAL;  // A27 opens         — counting begins
    double T_hg_pol_start  = DUMMY_VAL;  // A19 opens         — Hg polarisation fill
    double T_hg_pol_stop   = DUMMY_VAL;  // A19 closes        — polarisation done
};

// Parse firstTimeStamp from .hd header (line: firstTimeStamp = <value>L;)
static uint64_t ParseFirstTimestamp(const std::string& hd_path) {
    std::ifstream f(hd_path);
    std::string line;
    std::regex re(R"(firstTimeStamp\s*=\s*(\d+))");
    std::smatch m;
    while (std::getline(f, line)) {
        if (std::regex_search(line, m, re))
            return std::stoull(m[1].str());
    }
    return 0;
}

ValvesResult ReduceValves(const std::string& hd_path) {
    ValvesResult result;

    // Derive .EDMdat path: strip trailing _000 suffix before extension
    std::string dat_path = hd_path;
    const std::string old_suffix = "_000.hd";
    auto pos = dat_path.rfind(old_suffix);
    if (pos != std::string::npos)
        dat_path.replace(pos, old_suffix.size(), ".EDMdat");
    else {
        pos = dat_path.rfind(".hd");
        if (pos != std::string::npos)
            dat_path.replace(pos, 3, ".EDMdat");
    }

    // Use livecount firstTimeStamp as cycle t=0 — the neutron counter opens at the
    // true start of the cycle, while the valve file opens ~35 s later.
    std::string livecount_hd = hd_path;
    auto vpos = livecount_hd.rfind("valves_000.hd");
    if (vpos != std::string::npos)
        livecount_hd.replace(vpos, 13, "livecount_000.hd");
    uint64_t t0 = ParseFirstTimestamp(livecount_hd);
    if (t0 == 0)
        t0 = ParseFirstTimestamp(hd_path);  // fallback to valve file itself
    if (t0 == 0) {
        std::cerr << "ReduceValves: could not parse firstTimeStamp from " << hd_path << "\n";
        return result;
    }

    FILE* fp = fopen(dat_path.c_str(), "rb");
    if (!fp) return result;  // file absent for this cycle — leave fields as DUMMY_VAL

    uint64_t row[VALVES_NCOLS + 1];
    uint64_t prev_douta = 0;
    bool first_row = true;
    bool a27_closed = false;  // track whether we've seen A27 close (to distinguish open event)
    bool a29_closed_once = false;

    while (fread(row, sizeof(uint64_t), VALVES_NCOLS + 1, fp) == (size_t)(VALVES_NCOLS + 1)) {
        if (row[VALVES_NCOLS] != VALVES_SENTINEL) continue;

        uint64_t douta = row[COL_DOUTA];

        if (first_row) {
            prev_douta = douta;
            first_row = false;
            continue;
        }

        if (douta == prev_douta) continue;

        double t = (double)(row[COL_TIMESTAMP] - t0) / 1e9;

        auto changed = [](uint64_t p, uint64_t c, int b) { return ((p >> b) & 1) != ((c >> b) & 1); };
        auto bit     = [](uint64_t w, int b) -> int { return (int)((w >> b) & 1); };

        // A27: neutron shutter (N.O.) — 0=open, 1=closed
        if (changed(prev_douta, douta, BIT_A27)) {
            if (bit(douta, BIT_A27) == 1 && result.T_fill_end == DUMMY_VAL) {
                result.T_fill_end = t;   // first close = fill ends
                a27_closed = true;
            } else if (bit(douta, BIT_A27) == 0 && a27_closed && result.T_count_start == DUMMY_VAL) {
                result.T_count_start = t;  // first open after close = counting begins
            }
        }

        // A29: Hg shutter (N.O.) — 0=open, 1=closed (starts CLOSED, opens briefly)
        if (changed(prev_douta, douta, BIT_A29)) {
            if (bit(douta, BIT_A29) == 0 && result.T_fill_hg_start == DUMMY_VAL) {
                result.T_fill_hg_start = t;   // first open pulse
                a29_closed_once = false;
            } else if (bit(douta, BIT_A29) == 1 && !a29_closed_once && result.T_fill_hg_stop == DUMMY_VAL) {
                result.T_fill_hg_stop = t;    // close after first open
                a29_closed_once = true;
            }
        }

        // A19: Hg reservoir (normally closed) — 0=closed, 1=open
        if (changed(prev_douta, douta, BIT_A19)) {
            if (bit(douta, BIT_A19) == 1 && result.T_hg_pol_start == DUMMY_VAL)
                result.T_hg_pol_start = t;
            else if (bit(douta, BIT_A19) == 0 && result.T_hg_pol_stop == DUMMY_VAL)
                result.T_hg_pol_stop = t;
        }

        prev_douta = douta;
    }

    fclose(fp);
    return result;
}
