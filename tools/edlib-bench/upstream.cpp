#include "edlib.h"

extern "C" int upstream_distance(const char* query, int query_len,
                                 const char* target, int target_len, int k) {
    // The benchmark uses unbounded global and bounded prefix searches.
    const auto mode = k < 0 ? EDLIB_MODE_NW : EDLIB_MODE_SHW;
    const EdlibAlignConfig config{k, mode, EDLIB_TASK_DISTANCE, nullptr, 0};
    const auto result = edlibAlign(query, query_len, target, target_len, config);
    if (result.status != EDLIB_STATUS_OK) {
        return -2;
    }
    const int distance = result.editDistance;
    edlibFreeAlignResult(result);
    return distance;
}
