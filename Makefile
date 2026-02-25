CC      := gcc
CFLAGS  := -O2 -Wextra
LDLIBS  := -lm -lflint

# Common sources (NO sin/tan implementations here)
SRC_COMMON := \
  ./cmeasure/cbind_to_hw_thread.c \
  ./cmeasure/cmeasure.c \
  ./cmeasure/CrystalClockInC.c \
  ./tests/test_interface.c \
  ./tests/value_generation.c \
  ./tests/arb_comparison.c \
  ./util/bit_printing.c

# SIMD impl sources
SRC_AVX2   := ./vfast_sin_avx2.c ./vfast_tan_avx2.c
SRC_AVX512 := ./vfast_sin_avx512.c ./vfast_tan_avx512.c

AVX2FLAGS   := -mavx2 -mfma
AVX512FLAGS := -mavx512f -mavx512bw -mavx512dq -mavx512vl -mfma

BUILDDIR := build
OBJDIR_AVX2   := $(BUILDDIR)/avx2
OBJDIR_AVX512 := $(BUILDDIR)/avx512

OBJ_COMMON_AVX2   := $(patsubst ./%.c,$(OBJDIR_AVX2)/%.o,$(SRC_COMMON))
OBJ_COMMON_AVX512 := $(patsubst ./%.c,$(OBJDIR_AVX512)/%.o,$(SRC_COMMON))

OBJ_AVX2   := $(patsubst ./%.c,$(OBJDIR_AVX2)/%.o,$(SRC_AVX2))
OBJ_AVX512 := $(patsubst ./%.c,$(OBJDIR_AVX512)/%.o,$(SRC_AVX512))

.PHONY: all avx2 avx512 clean

all: avx2 avx512

avx2: trig_simd_avx2
avx512: trig_simd_avx512

trig_simd_avx2: $(OBJ_COMMON_AVX2) $(OBJ_AVX2)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

trig_simd_avx512: $(OBJ_COMMON_AVX512) $(OBJ_AVX512)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

# Compile rules (same source path, different obj dirs + flags)
$(OBJDIR_AVX2)/%.o: ./%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(AVX2FLAGS) -c $< -o $@

$(OBJDIR_AVX512)/%.o: ./%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(AVX512FLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR) trig_simd_avx2 trig_simd_avx512
