/* Determine the highest SIMD capability  */

/* Intel/AMD core SIMD with 8 double float width  */
#if defined(__AVX512F__)
  #include <immintrin.h>
  #define _simd_type __m512d
  #define _simd_loadu_pd _mm512_loadu_pd
  #define _simd_set1_pd _mm512_set1_pd
  #define _simd_add_pd _mm512_add_pd
  #define _simd_sub_pd _mm512_sub_pd
  #define _simd_mul_pd _mm512_mul_pd
  #define _simd_div_pd _mm512_div_pd
  #define _simd_storeu_pd _mm512_storeu_pd
  #define _simd_setzero_pd _mm512_setzero_pd
  #define SIMD_WIDTH 8

/* Intel/AMD core SIMD with 4 double float width  */
#elif defined(__AVX__)
  #include <immintrin.h>
  #define _simd_type __m256d
  #define _simd_loadu_pd _mm256_loadu_pd
  #define _simd_set1_pd _mm256_set1_pd
  #define _simd_add_pd _mm256_add_pd
  #define _simd_sub_pd _mm256_sub_pd
  #define _simd_mul_pd _mm256_mul_pd
  #define _simd_div_pd _mm256_div_pd
  #define _simd_storeu_pd _mm256_storeu_pd
  #define _simd_setzero_pd _mm256_setzero_pd
  #define SIMD_WIDTH 4

/* Intel/AMD core SIMD with 2 double float width  */
#elif defined(__SSE2__)
  #include <emmintrin.h>
  #define _simd_type __m128d
  #define _simd_loadu_pd _mm_loadu_pd
  #define _simd_set1_pd _mm_set1_pd
  #define _simd_add_pd _mm_add_pd
  #define _simd_sub_pd _mm_sub_pd
  #define _simd_mul_pd _mm_mul_pd
  #define _simd_div_pd _mm_div_pd
  #define _simd_storeu_pd _mm_storeu_pd
  #define _simd_setzero_pd _mm_setzero_pd
  #define SIMD_WIDTH 2

/* No SIMD detected */
#else
  #define SIMD_WIDTH 1
#endif
