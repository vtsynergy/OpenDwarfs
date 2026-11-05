/* vim: set sts=4 sw=4 expandtab:*/
/*
 * =====================================================================================
 *
 *       Filename:  cmain.c
 *
 *    Description:  wrapper to make gem have a c based main... don't ask
 *
 *        Version:  1.0
 *        Created:  04/14/12 21:19:57
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tom Scogland (), 
 *        Company:  
 *
 * =====================================================================================
 */
#include "ctsar.h"
#include "gem.h"

int main(int argc, char **argv){

    ctsar_pre_init();
    return cpp_main(argc, argv);
}
