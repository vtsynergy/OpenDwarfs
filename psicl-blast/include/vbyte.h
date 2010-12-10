/* macro to read a vbyte number from cmem (which is incremented over the number)
 * into n.  Note that cmem must be a uchar * */

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

#define vbyte_getVbyte(cmem, n)                                                  \
    if (1) {                                                                  \
        int B2_GETVBYTE_count = 0;                                   \
                                                                              \
        *(n) = 0;                                                             \
        while (((unsigned char)(*(cmem)) >= 128)) {                             \
            *(n) |= (*(cmem) & 0x7f) << B2_GETVBYTE_count;                    \
            B2_GETVBYTE_count += 7;                                           \
            (cmem)++;                                                         \
        }                                                                     \
        *(n) |= (unsigned char)(*(cmem)) << B2_GETVBYTE_count;                  \
        (cmem)++;                                                             \
    } else

/* macro to put a vbyte number to cmem (which is incremented over the number)
 * from n (which is destroyed by the operation).  Note that cmem must be
 * a uchar * */
#define vbyte_putVbyte(cmem, n)                                                  \
    if (1) {                                                                  \
    while ((n) >= 128) {                                                      \
            *(cmem) = ((n) & 0x7f) | 0x80;                                    \
            (n) >>= 7;                                                        \
            (cmem)++;                                                         \
        }                                                                     \
        *(cmem) = ((unsigned char) n);                                        \
        (cmem)++;                                                             \
    } else

// Same as vbyte_putVbyte but does not destroy n
#define vbyte_safePutVbyte(cmem, originalN)                                   \
    if (1) {                                                                  \
    uint8 n = originalN;                                                        \
    while ((n) >= 128) {                                                      \
            *(cmem) = ((n) & 0x7f) | 0x80;                                    \
            (n) >>= 7;                                                        \
            (cmem)++;                                                         \
        }                                                                     \
        *(cmem) = ((unsigned char) n);                                        \
        (cmem)++;                                                             \
    } else

// Reads a vbyte into a 64-bit integer
unsigned char* vbyte_get64vbyte(unsigned char* cmem, uint8* n);

#ifdef __MY_EXTERN_C__
}
#endif

