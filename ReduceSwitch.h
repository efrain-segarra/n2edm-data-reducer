#pragma once

#include <string>
#include <fstream>
#include <regex>
#include <cstdint>
#include <cstdio>
#include <iostream>

#include "constants.h"

// Binary format: 17 × uint64_t per row + uint64_t sentinel = 144 bytes/row
static constexpr int      SWITCH_NCOLS     = 17;
static constexpr int      SWITCH_ROW_BYTES = (SWITCH_NCOLS + 1) * 8;
static constexpr uint64_t SWITCH_SENTINEL  = 0xDEADBEEFDEADBEEFULL;

// Column indices
static constexpr int SW_COL_TIMESTAMP = 0;
static constexpr int SW_COL_STATUS1   = 9;
static constexpr int SW_COL_EXECUTE1  = 13;

// status1 values
static constexpr uint64_t SW_STATUS_IDLE      = 0;
static constexpr uint64_t SW_STATUS_RETURNING = 2;
static constexpr uint64_t SW_STATUS_MOVING    = 5;

struct SwitchResult {
    double Switch_count_start    = DUMMY_VAL;  // switch begins moving to count position
    double Switch_count_duration = DUMMY_VAL;  // seconds until settled at count position
    double Switch_fill_start     = DUMMY_VAL;  // switch begins returning to fill position
    double Switch_fill_duration  = DUMMY_VAL;  // seconds until settled at fill position
};

static uint64_t ParseFirstTimestamp_sw(const std::string& hd_path) {
    std::ifstream f(hd_path);
    std::string line;
    std::regex re(R"(firstTimeStamp\s*=\s*(\d+))");
    std::smatch m;
    while (std::getline(f, line))
        if (std::regex_search(line, m, re))
            return std::stoull(m[1].str());
    return 0;
}

SwitchResult ReduceSwitch(const std::string& hd_path) {
    SwitchResult result;

    // Derive .EDMdat path: strip trailing _000 from stem
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

    // Use livecount firstTimeStamp as cycle t=0
    std::string livecount_hd = hd_path;
    auto vpos = livecount_hd.rfind("switch_000.hd");
    if (vpos != std::string::npos)
        livecount_hd.replace(vpos, 13, "livecount_000.hd");
    uint64_t t0 = ParseFirstTimestamp_sw(livecount_hd);
    if (t0 == 0)
        t0 = ParseFirstTimestamp_sw(hd_path);
    if (t0 == 0) {
        std::cerr << "ReduceSwitch: could not parse firstTimeStamp from " << hd_path << "\n";
        return result;
    }

    FILE* fp = fopen(dat_path.c_str(), "rb");
    if (!fp) return result;

    uint64_t row[SWITCH_NCOLS + 1];
    uint64_t prev_execute1 = 0;
    bool first_row = true;

    while (fread(row, sizeof(uint64_t), SWITCH_NCOLS + 1, fp) == (size_t)(SWITCH_NCOLS + 1)) {
        if (row[SWITCH_NCOLS] != SWITCH_SENTINEL) continue;

        uint64_t execute1 = row[SW_COL_EXECUTE1];
        uint64_t status1  = row[SW_COL_STATUS1];

        if (first_row) {
            prev_execute1 = execute1;
            first_row = false;
            continue;
        }

        if (execute1 == prev_execute1) continue;

        double t = (double)(row[SW_COL_TIMESTAMP] - t0) / 1e9;

        if (execute1 == 1) {
            // Switch started a move — classify by status
            if (status1 == SW_STATUS_MOVING && result.Switch_count_start == DUMMY_VAL)
                result.Switch_count_start = t;
            else if (status1 == SW_STATUS_RETURNING && result.Switch_fill_start == DUMMY_VAL)
                result.Switch_fill_start = t;
        } else {
            // Switch settled (execute1 went 1→0)
            if (result.Switch_count_start != DUMMY_VAL && result.Switch_count_duration == DUMMY_VAL)
                result.Switch_count_duration = t - result.Switch_count_start;
            else if (result.Switch_fill_start != DUMMY_VAL && result.Switch_fill_duration == DUMMY_VAL)
                result.Switch_fill_duration = t - result.Switch_fill_start;
        }

        prev_execute1 = execute1;
    }

    fclose(fp);
    return result;
}
