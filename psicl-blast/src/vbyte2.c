/* vbyte.c implements variable byte encoding to stdio files
 * 
 * written Hugh Williams, August 1995 
 * updated Dirk, March 2000
 * update nml, 2003-03-07
 *
 */

/*#include "firstinclude.h"*/

#include "vbyte2.h"

#include <errno.h>

#include <stdint.h>
#include <limits.h>

/* XXX: names for bitmasks: form bits.h */
unsigned int vbyte_read(FILE* fp, unsigned long int* n) {
    unsigned int ret = 1;
    unsigned int count = 0; 
    uint8_t byte;

    *n = 0;
    while (fread(&byte, 1, 1, fp)) {
        if (byte & 0x80) {
            /* For each time we get a group of 7 bits, need to left-shift the 
               latest group by an additional 7 places, since the groups were 
               effectively stored in reverse order.*/
            *n |= ((byte & 0x7f) << count);
            ret++;
            count += 7;
        } else if (ret > (sizeof(*n))) {
            if (ret == (sizeof(*n) + 1) 
              && (byte < (1 << (sizeof(*n) + 1)))) {
                /* its a large, but valid number */
                *n |= byte << count;
                return ret;
            } else {
                /* we've overflowed */

                /* pretend we detected it as it happened */
                fseek(fp, - (int) (ret - (sizeof(*n) + 1)), SEEK_CUR);
                errno = EOVERFLOW;
                return 0;
            }
        } else {
            *n |= byte << count;
            return ret;
        }
    }

    /* got a read error */
    return 0;
}

unsigned int vbyte_write(FILE* fp, unsigned long int n) {
    unsigned int ret = 1;
    uint8_t byte;

    while (n >= 128) {
        /* write the bytes out least significant to most significant */
        byte = (uint8_t) (n & 0x7f) | 0x80;
        if (!fwrite(&byte, 1, 1, fp)) {
            return 0;
        }
        n >>= 7;
        ret++;
    }

    byte = (uint8_t) n;
    if (!fwrite(&byte, 1, 1, fp)) {
        return 0;
    }
    
    return ret;
}

unsigned int vbyte_scan(FILE* fp) {
    unsigned int ret = 1;
    uint8_t byte;

    while (fread(&byte, 1, 1, fp)) {
        if (byte & 0x80) {
            /* For each time we get a group of 7 bits, need to left-shift the 
               latest group by an additional 7 places, since the groups were 
               effectively stored in reverse order.*/
            ret++;
        } else if (ret > sizeof(unsigned long int)) {
            if (ret == (sizeof(unsigned long int) + 1) 
              && (byte < (1 << (sizeof(unsigned long int) + 1)))) {
                /* its a large, but valid number */
                return ret;
            } else {
                /* we've overflowed */

                /* pretend we detected it as it happened */
                fseek(fp, - (int) (ret - (sizeof(unsigned long int) + 1)), SEEK_CUR);
                errno = EOVERFLOW;
                return 0;
            }
        } else {
            return ret;
        }
    }

    /* got a read error */
    return 0;
}

unsigned int vbyte_len(unsigned long int n) {
    unsigned int ret = 1;
    while (n >= 128) {
        n >>= 7;
        ret++;
    }
    return ret;
}

