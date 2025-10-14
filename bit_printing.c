#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdint.h> 
#include "bit_printing.h"

void print_double_bits(double value) {
  uint64_t bits;
  union {
    double d;
    uint64_t u;
  } conv;

  conv.d = value;
  bits = conv.u;
  printf("%d ", (int)((bits >> 63) & 1));

  for (int i = 62; i >= 52; i--) {
    printf("%d", (int)((bits >> i) & 1));
  }

  printf(" ");

  for (int i = 51; i >= 0; i--) {
    printf("%d", (int)((bits >> i) & 1));
  }

  printf("\n");
}

void print_bits_ulong(unsigned long x) {
    int nbits = sizeof(unsigned long) * 8;
    for (int i = nbits - 1; i >= 0; i--) {
        putchar( (x & (1ul << i)) ? '1' : '0' );
        // optional: add a space or separator every 8 bits
        if (i % 8 == 0 && i != 0) putchar(' ');
    }
    printf("\n");
}

void print_bits_u8(uint8_t x) {
    for (int i = 7; i >= 0; i--) {
        putchar( (x & (1u << i)) ? '1' : '0' );
    }
    printf("\n");
}
