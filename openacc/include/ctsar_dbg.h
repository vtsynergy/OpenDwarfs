#ifndef __CTSAR_DBG
#define __CTSAR_DBG

enum CTSAR_DBG_LEVEL{
    DBG_ALWAYS=0,
    DBG_INFO=1,
    DBG_TEST=2,
    DBG_DEBUG=3
};
extern int sdebug;

#define DBG(level, ...) if(sdebug > level) fprintf(stderr, "ctsar: " __VA_ARGS__);

#endif // __CTSAR_DBG
