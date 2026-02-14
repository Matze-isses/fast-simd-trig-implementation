CC      := gcc
CFLAGS  := -O2 -Wextra
LDLIBS  := -lm -lflint

SRC := \
  ./cmeasure/cbind_to_hw_thread.c \
  ./cmeasure/cmeasure.c \
  ./cmeasure/CrystalClockInC.c \
  ./tests/test_interface.c \
  ./tests/value_generation.c \
  ./tests/arb_comparison.c \
  ./util/bit_printing.c \
  ./tan.c \
  ./sin.c

# AVX2 flags (matches your command; -mavx is redundant if -mavx2 is present but kept)
AVX2FLAGS   := -mavx -mavx2 -mfma

# AVX-512 flags (common, practical subset)
# -mavx512f is the foundation; BW/DQ/VL are widely useful and common on Intel AVX-512 parts.
AVX512FLAGS := -mavx512f -mavx512bw -mavx512dq -mavx512vl -mfma

.PHONY: all avx2 avx512 clean

all: avx2 avx512

avx2: trig_simd_avx2

avx512: trig_simd_avx512

trig_simd_avx2: $(SRC)
	$(CC) $(CFLAGS) $(AVX2FLAGS) $(SRC) -o $@ $(LDLIBS)

trig_simd_avx512: $(SRC)
	$(CC) $(CFLAGS) $(AVX512FLAGS) $(SRC) -o $@ $(LDLIBS)

clean:
	rm -f trig_simd_avx2 trig_simd_avx512
