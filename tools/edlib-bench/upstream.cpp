#include "edlib.h"

static bool valid_path(const char* query, int query_len,
                       const char* target, int target_len,
                       EdlibAlignMode mode, const EdlibAlignResult& result) {
    if (!result.alignment || result.numLocations < 1 ||
        !result.startLocations || !result.endLocations) return false;
    const int start = result.startLocations[0];
    const int end = result.endLocations[0];
    if (start < 0 || end < start - 1 || end >= target_len ||
        (mode != EDLIB_MODE_HW && start != 0) ||
        (mode == EDLIB_MODE_NW && end != target_len - 1)) return false;
    const int span_len = end + 1 - start;
    target += start;
    int q = 0, t = 0, cost = 0;
    for (int i = 0; i < result.alignmentLength; ++i) {
        const auto op = result.alignment[i];
        switch (op) {
            case EDLIB_EDOP_MATCH:
            case EDLIB_EDOP_MISMATCH:
                if (q == query_len || t == span_len ||
                    (query[q] == target[t]) != (op == EDLIB_EDOP_MATCH)) return false;
                ++q;
                ++t;
                break;
            case EDLIB_EDOP_INSERT:
                if (q == query_len) return false;
                ++q;
                break;
            case EDLIB_EDOP_DELETE:
                if (t == span_len) return false;
                ++t;
                break;
            default:
                return false;
        }
        cost += op != EDLIB_EDOP_MATCH;
    }
    return q == query_len && t == span_len && cost == result.editDistance;
}

extern "C" int upstream_align(const char* query, int query_len,
                              const char* target, int target_len, int k,
                              int mode_id, bool path, bool check_path) {
    static_assert(EDLIB_MODE_NW == 0 && EDLIB_MODE_SHW == 1 && EDLIB_MODE_HW == 2,
                  "mode IDs must match the Rust benchmark adapter");
    const auto mode = static_cast<EdlibAlignMode>(mode_id);
    const auto task = path ? EDLIB_TASK_PATH : EDLIB_TASK_DISTANCE;
    const EdlibAlignConfig config{k, mode, task, nullptr, 0};
    const auto result = edlibAlign(query, query_len, target, target_len, config);
    if (result.status != EDLIB_STATUS_OK) {
        return -2;
    }
    const int distance = check_path && !valid_path(query, query_len, target, target_len, mode, result)
        ? -3 : result.editDistance;
    edlibFreeAlignResult(result);
    return distance;
}
