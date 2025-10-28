// clock_utils.h
#ifndef CLOCK_UTILS_H
#define CLOCK_UTILS_H

#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>

static inline uint64_t __cu_ts_to_ns(const struct timespec *ts) {
    return (uint64_t)ts->tv_sec * 1000000000ull + (uint64_t)ts->tv_nsec;
}

// Single active timer (per TU/thread) – START then END.
static struct timespec __cu_start_ts;

#define START_CLOCK do { \
    clock_gettime(CLOCK_MONOTONIC, &__cu_start_ts); \
} while (0)

#define END_CLOCK(name_literal) do { \
    struct timespec __cu_end_ts; \
    clock_gettime(CLOCK_MONOTONIC, &__cu_end_ts); \
    uint64_t __cu_start_ns = __cu_ts_to_ns(&__cu_start_ts); \
    uint64_t __cu_end_ns   = __cu_ts_to_ns(&__cu_end_ts); \
    uint64_t __cu_delta_ns = __cu_end_ns - __cu_start_ns; \
    double   __cu_ms       = __cu_delta_ns / 1e6; \
    printf("%s: %.3f ms (%" PRIu64 " ns)\n", (name_literal), __cu_ms, __cu_delta_ns); \
} while (0)


// ------------------------------------------------------------
// Additional macros for getting elapsed time directly
// ------------------------------------------------------------
static struct timespec __cu_tclock_start_ts;

#define START_TCLOCK do { \
    clock_gettime(CLOCK_MONOTONIC, &__cu_tclock_start_ts); \
} while (0)

#define GET_TCLOCK ( \
    ({ \
        struct timespec __cu_tclock_end_ts; \
        clock_gettime(CLOCK_MONOTONIC, &__cu_tclock_end_ts); \
        uint64_t __cu_tclock_start_ns = __cu_ts_to_ns(&__cu_tclock_start_ts); \
        uint64_t __cu_tclock_end_ns   = __cu_ts_to_ns(&__cu_tclock_end_ts); \
        uint64_t __cu_tclock_delta_ns = __cu_tclock_end_ns - __cu_tclock_start_ns; \
        (double)__cu_tclock_delta_ns / 1e6; \
    }) \
)

#endif // CLOCK_UTILS_H
