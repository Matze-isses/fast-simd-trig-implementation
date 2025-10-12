// Intrinsics Macro Header

#ifndef J512_INTRINSICS_H
#define J512_INTRINSICS_H
#ifdef __cplusplus
extern "C" {
#endif

// }
#define INT64S __m512i 
#define INT32S __m256i
#define DOUBLES __m512d 
#define MASK __mask8

#define LOAD_DOU _mm512_loadu_pd
#define LOAD_INT _mm512_loadu_si512
#define LOAD_32 _mm256_loadu_epi32

#define ADD _m512_add_pd
#define SUB _m512_sub_pd
#define MUL _mm512_mul_pd
#define FMADD _mm512_fmadd_pd
#define FMSUB _mm512_fmsub_pd
#define FNMADD _mm512_fnmadd_pd
#define STORE_DOU _mm512_storeu_pd
#define ZERO_DOU _mm512_setzero_pd
#define SET_DOU _mm512_set1_pd

#define MASK_MOV _mm512_mask_mov_pd
#define MASKZ_MOV _mm512_maskz_mov_pd
#define MASK_BLEND _mm512_mask_blend_pd
#define MASKOADDDOU _mm512_maskz_add_pd
#define MASKMULDOU(A, B, C, D) _mm512_mask_mul_pd(A, B, C, D)
#define MASKADDDOU(A, B, C, D) _mm512_mask_add_pd(A, B, C, D)

#define OR _mm512_or_pd
#define AND _mm512_and_pd
#define ANDNOT _mm512_and_pd
#define AND64 _mm512_and_si512
#define AND32 _mm256_and_si256
#define OR64 _mm512_or_si512
#define XOR64 _mm512_xor_si512

#define SHR32 _mm512_srli_epi32
#define SHL32 _mm512_slli_epi32

#define SHR64 _mm512_srli_epi64
#define SHL64 _mm512_slli_epi64

#define SHR25632 _mm256_srli_epi32
#define ADD32 _mm256_add_epi32

#define MAX32 _mm256_max_epi32
#define SET16 _mm512_set1_epi16
#define SET32 _mm512_set1_epi32
#define SET64 _mm512_set1_epi64

#define ZERO32 _mm256_setzero_si256
#define MADD16 _mm512_madd_epi16
#define ADD64 _mm512_add_epi64

#define MULT32 _mm512_mullo_epi32
#define MULT64 _mm512_mullo_epi64

#define SUB32 _mm256_sub_epi32
#define SUB64 _mm512_sub_epi64

#define CAST_DOU_INT _mm512_castpd_si512
#define CAST_INT_DOU _mm512_casts512_pd
#define CUT32 _mm512_cvtepi64_epi32
#define CUTD32 _mm512_cvttpd_epi32
#define CUT64 _mm512_cvtepi32_epi64
#define TESTAND32 _mm256_test_epi32_mask
#define TESTNAND32 _mm256_testn_epi32_mask
#define BLENDvDOUBLE _mm512_blendv_pd

#define MAND _kand_mask8
#define MXOR _kxor_mask8
#define MNOT _knot_mask8
#define MOR _kor_mask8

#ifdef __cplusplus
}
#endif

#endif










