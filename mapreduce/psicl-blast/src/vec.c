/* vec.c implements the byte vector described in vec.h 
 *
 * based on code by Hugh Williams
 *
 * written nml 2003-02-28
 *
 */

#include "_vec.h"  

#include "vec.h"  
#include "vbyte.h" 

#include <stdint.h>
#include <limits.h>
#include <fcntl.h>
#include <unistd.h>

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

int vec_init_inplace(unsigned long int initsize, struct vec *v) {
    if ((v->vector = malloc(initsize))) {
        v->size = initsize;
        v->pos = v->len = 0;
        v->err = 0;
        return 1;
    } else {
        return 0;
    }
}

struct vec* vec_init(unsigned long int initsize) {
    struct vec* v = malloc(sizeof(*v));

    if (v && (v->vector = malloc(initsize))) {
        v->size = initsize;
        v->pos = v->len = 0;
        v->err = 0;
        return v;
    } else if (v) {
        free(v);
    }

    return NULL;
}

struct vec* vec_reinit(struct vec* v, unsigned long int initsize) {
    void* newmem;

    if (!v) {
        /* need to allocate a new vector */
        if ((v = malloc(sizeof(*v))) 
          && (v->vector = malloc(initsize))) {
            v->size = initsize;
            v->pos = v->len = 0;
            v->err = 0;
        } else if (v) {
            free(v);
            return NULL;
        } else {
            return NULL;
        }
    } else {
        /* reinit old one */
        if ((newmem = realloc(v->vector, initsize))) {
            v->size = initsize;
            v->vector = newmem;
            v->pos = v->len = 0;
            v->err = 0;
        } else {
            return NULL;
        }
    }

    return v;
}

int vec_free(struct vec* v) {
    free(v->vector);
    free(v);
    return 1;
}

void vec_resetpos(struct vec* v) {
    v->pos = 0;
}

void vec_reseterr(struct vec* v) {
    v->err = 0;
}

unsigned int vec_append(struct vec* dst, const struct vec* src) {
    unsigned int size = src->len - src->pos;

    /* expand vector if required */
    while ((dst->size < dst->pos + size) && vec_expand(dst, dst->pos + size)) ;

    if (!dst->err && (dst != src)) {
        /* append stuff */
        memcpy(dst->vector + dst->pos, src->vector + src->pos, size);
        dst->pos += size;
        dst->len = dst->pos;
        return size;
    } else {
        return 0;
    }
}

/* XXX: names for bitmasks: form bits.h */
unsigned long int vec_putvbyte(struct vec* v, unsigned long int n) {
    unsigned long int ret = 1;

    while (n >= 128) {
        if ((v->pos >= v->size) && !vec_expand(v, GROWTHFACTOR * v->size + 1)) {
            v->err = ENOSPC;
            return 0;
        }
        
        v->vector[v->pos++] = (char) ((n & 0x7f) | 0x80);
        n >>= 7;
        ret++;
    }

    if ((v->pos >= v->size) && !vec_expand(v, GROWTHFACTOR * v->size + 1)) {
        v->err = ENOSPC;
        return 0;
    }
    
    v->vector[v->pos++] = (char) n;
    v->len = v->pos;                 /* truncate vector */
    
    return ret;
}

unsigned long int vec_getvbyte(struct vec* v, unsigned long int* n) {
    unsigned long int ret = 1, 
                      count = 0, 
                      get = 0;

    *n = 0;
    while (v->pos < v->len) {
        if (((get = ((uint8_t) v->vector[v->pos++])) & 0x80) == 0x80) {
            /* For each time we get a group of 7 bits, need to left-shift the 
               latest group by an additional 7 places, since the groups were 
               effectively stored in reverse order.*/
            *n |= ((get & 0x7f) << count);
            ret++;
            count += 7;
        } else if (ret > (sizeof(*n))) {
            if (ret == ((sizeof(*n) + 1)) 
              && (get < (1 << (sizeof(*n) + 1)))) {
                /* its a large, but valid number */
                *n |= get << count;
                return ret;
            } else {
                /* we've overflowed, return error */
                return 0;
            }
        } else {
            /* Now get the final 7 bits, which need to be left-shifted by a 
               factor of 7. */
            *n |= get << count;
            return ret;
        }
    }

    v->err = ENOSPC;
    return 0;
}

unsigned long int vec_scanvbyte(struct vec* v) {
    unsigned long int ret = 1;

    while ((v->pos < v->len)) {
        if (!(v->vector[v->pos++] & 0x80)) { 
            return ret;
        }
        ret++;
    }

    v->err = ENOSPC;
    return 0;
}

int vec_expand(struct vec* v, unsigned int size) {
    /* grow the buffer (+ 1 is to step around case where size is 0) */
    void* newmem = realloc(v->vector, size);

    if (newmem) {
        v->size = size;
        v->vector = newmem;
        return 1;
    } else {
        v->err = ENOMEM;
        return 0;
    }
}

int vec_err(const struct vec* v) {
    return v->err;
}

int vec_eof(const struct vec* v) {
    return (v->pos == v->len);
}

unsigned long int vec_position(const struct vec* v) {
    return v->pos;
}

unsigned long int vec_length(const struct vec* v) {
    return v->len;
}

struct vec* vec_copy(const struct vec* v) {
    struct vec* nv = malloc(sizeof(*nv));

    if (nv && (nv->vector = malloc(v->len))) {
        nv->pos = v->pos;
        nv->len = v->len;
        nv->size = v->len;
        nv->err = 0;
        memcpy(nv->vector, v->vector, nv->size);
    } else if (nv) {
        free(nv);
        return NULL;
    }

    return nv;
}

int vec_read_inplace(int srcfd, struct vec *v, unsigned int len) {
    if ((v->vector = malloc(len))) {
        unsigned int rlen;
        char *dst = v->vector;
        unsigned int size = len;

        while (len && ((rlen = read(srcfd, dst, len)) > 0)) {
            len -= rlen;
            dst += rlen;
        }

        if (!len) {
            v->len = v->size = size;
            v->pos = 0;
            v->err = 0;
            return 1;
        } else {
            free(v->vector);
            return 0;
        } 
    } else {
        return 0;
    }
}

struct vec* vec_read(int srcfd, unsigned int len) {
    struct vec* v;

    if ((v = malloc(sizeof(*v))) && vec_read_inplace(srcfd, v, len)) {
        return v;
    } else {
        if (v) {
            free(v);
        }
        return NULL;
    }
}

unsigned int vec_write(struct vec* src, int dstfd) {
    unsigned int wlen,
                 len = src->len - src->pos;
    char *srcbuf = src->vector;


    while (len && ((wlen = write(dstfd, srcbuf, len)) >= 0)) {
        len -= wlen;
        srcbuf += wlen;
    }

    if (!len) {
        return src->len - src->pos;
    } else {
        src->err = errno;
        return 0;
    }
}

int vec_cmp(const struct vec *one, const struct vec *two) {
    if (one->len != two->len) {
        return one->len - two->len;
    } else {
        return (memcmp(one->vector, two->vector, one->len));
    }
}

unsigned int vec_vbyte_len(unsigned long int n) {
    unsigned int ret = 1;
    while (n >= 128) {
        n >>= 7;
        ret++;
    }
    return ret;
}

