/* Main file for CRC (combinational-logic dwarf) -- OpenACC translation of
 * opencl/combinational-logic/crc/src/crc_algo.c + crc_kernel.cl.
 *
 * Computes a 32-bit Ethernet CRC on a number of input pages with the "Slice-by-8" algorithm
 * (Stephan Brumme / Intel). Same input reader (read_crc), same blocking of pages into blocks of
 * `-p` parallel CRCs (default 16), one work-item per page, same kernel body, same serial
 * slice-by-8 check (-a). CLI: -i file [-a] [-v] [-p n] [-r n] [-w n] [-k ignored].
 * As in the OpenCL program, the device buffers are created outside the timed region; the timed
 * region is, for every block, H2D of its pages, the kernel and D2H of its CRCs, queued
 * asynchronously (the non-blocking enqueues) and followed by one wait (clFinish).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <openacc.h>
#include "ocd_compat.h"
#include "../inc/crc_formats.h"
#include "dwarf_timer.h"
#include "dwarf_acc.h" /* [FIX-OCL-34] pinned launch configuration */

/* the lookup table lives on the device as well (OpenCL: __constant) */
#include "../inc/eth_crc32_lut.h"
#pragma acc declare copyin(crc32Lookup)

unsigned char verbosity = 0;

/* Slice-by-8 CRC on the host (serial check), altered from Stephan Brumme's crc32 page */
static uint32_t crc32_8bytes(const void *data, size_t length) {
  uint32_t *current = (uint32_t *)data;
  uint32_t crc = 0xFFFFFFFF;
  while (length >= 8) {
    uint32_t one = *current++ ^ crc;
    uint32_t two = *current++;
    crc = crc32Lookup[7][one & 0xFF] ^ crc32Lookup[6][(one >> 8) & 0xFF] ^ crc32Lookup[5][(one >> 16) & 0xFF] ^
          crc32Lookup[4][one >> 24] ^ crc32Lookup[3][two & 0xFF] ^ crc32Lookup[2][(two >> 8) & 0xFF] ^
          crc32Lookup[1][(two >> 16) & 0xFF] ^ crc32Lookup[0][two >> 24];
    length -= 8;
  }
  unsigned char *currentChar = (unsigned char *)current;
  while (length--) crc = (crc >> 8) ^ crc32Lookup[0][(crc & 0xFF) ^ *currentChar++];
  return ~crc;
}

/* One work-item of crc_kernel.cl (crc32_slice8): CRC of page `gid` of this block. */
#pragma acc routine seq
static unsigned int crc_item(const unsigned int *data, unsigned int length_bytes, unsigned int length_ints,
                             size_t gid) {
  unsigned int crc = 0xFFFFFFFF, one, two;
  size_t i = gid * length_ints, j;
  while (length_bytes >= 8) {  /* process eight bytes at once */
    one = data[i++] ^ crc;
    two = data[i++];
    crc = crc32Lookup[7][one & 0xFF] ^ crc32Lookup[6][(one >> 8) & 0xFF] ^ crc32Lookup[5][(one >> 16) & 0xFF] ^
          crc32Lookup[4][one >> 24] ^ crc32Lookup[3][two & 0xFF] ^ crc32Lookup[2][(two >> 8) & 0xFF] ^
          crc32Lookup[1][(two >> 16) & 0xFF] ^ crc32Lookup[0][two >> 24];
    length_bytes -= 8;
  }
  while (length_bytes) {  /* remaining 1 to 7 bytes */
    one = data[i++];
    const unsigned char *currentChar = (const unsigned char *)&one;
    j = 0;
    while (length_bytes && j < 4) {
      length_bytes = length_bytes - 1;
      crc = (crc >> 8) ^ crc32Lookup[0][(crc & 0xFF) ^ currentChar[j]];
      j = j + 1;
    }
  }
  return ~crc;
}

/* enqueueCRCDevice(): write the block's pages, run the kernel over global_size work-items, read the
 * CRCs back, all on one async queue. crc_kernel.cl uses neither barriers nor local memory; the
 * pinned shape is local 128, global roundup(pages) (the loop bound replaces the gid >= num_items guard). */
static void crc_block(unsigned int *h_num, unsigned int *h_answer, size_t global_size,
                      unsigned int page_size, unsigned int num_words) {
  size_t nbytes_words = (size_t)num_words * global_size;
#pragma acc update device(h_num[0 : nbytes_words]) async(1)
#pragma acc parallel loop gang vector ACC_LAUNCH(global_size) present(h_num[0 : nbytes_words], h_answer[0 : global_size]) async(1)
  for (size_t gid = 0; gid < global_size; gid++) h_answer[gid] = crc_item(h_num, page_size, num_words, gid);
#pragma acc update self(h_answer[0 : global_size]) async(1)
}

static void usage(void) {
  printf("crc -i <input_file> [-a] [-v] [-p <parallel crcs>] [-r <num_execs>] [-w <wg_size>]\n");
  exit(0);
}

int main(int argc, char **argv) {
  size_t global_size;
  unsigned int *h_num;
  unsigned int run_serial = 0, num_pages = 1, num_execs = 1, page_size = 0;
  unsigned int num_wg_sizes = 0, num_block_sizes = 0;
  unsigned int *num_parallel_crcs = NULL;
  size_t *wg_sizes = NULL;
  char *file = NULL;
  int c;
  void *tmp;

  while ((c = getopt(argc, argv, "avn:s:i:p:w:k:hr:")) != -1) {
    switch (c) {
      case 'h': usage(); break;
      case 'v': verbosity++; break;
      case 'a': run_serial = 1; break;
      case 'i': file = optarg; printf("Reading Input from '%s'\n", file); break;
      case 'r': num_execs = atoi(optarg); printf("Executing %d times\n", num_execs); break;
      case 'p':
        num_block_sizes++;
        tmp = realloc(num_parallel_crcs, sizeof(unsigned int) * num_block_sizes);
        check(tmp != NULL, "Heap Overflow! num_parallel_crcs");
        num_parallel_crcs = tmp;
        num_parallel_crcs[num_block_sizes - 1] = atoi(optarg);
        break;
      case 'w':
        num_wg_sizes++;
        tmp = realloc(wg_sizes, sizeof(size_t) * num_wg_sizes);
        check(tmp != NULL, "Heap Overflow! wg_sizes");
        wg_sizes = tmp;
        wg_sizes[num_wg_sizes - 1] = atoi(optarg);
        break;
      case 'k': printf("Note: -k ignored (single OpenACC kernel)\n"); break;
      default: fprintf(stderr, "Invalid argument\n\n"); usage();
    }
  }
  check(file != NULL, "-i option must be supplied!");
  h_num = read_crc(&num_pages, &page_size, file);
  if (!num_block_sizes) { num_block_sizes = 1; num_parallel_crcs = malloc(sizeof(unsigned int)); num_parallel_crcs[0] = 16; }
  if (!num_wg_sizes) { num_wg_sizes = 1; wg_sizes = malloc(sizeof(size_t)); wg_sizes[0] = 1; }
  unsigned int num_words = page_size / 4;
  if (verbosity) printf("num_words = %u\n", num_words);

  unsigned int *remainders = (unsigned int *)calloc(num_pages, sizeof(unsigned int));
  unsigned int bad = 0;
  acc_init(acc_get_device_type()); /* device discovery/initialisation, outside the timed region */
  dwarf_init("crc", "combinational-logic", "openacc");
  dwarf_set_size(num_pages, page_size);
  if (num_execs > 1) dwarf_set_repeats((int)num_execs);

  for (unsigned int h = 0; h < num_block_sizes; h++) {
    unsigned int per = num_parallel_crcs[h], num_blocks = num_pages / per, num_pages_last_block;
    if (num_pages % per != 0) { num_blocks++; num_pages_last_block = num_pages % per; }
    else num_pages_last_block = per;
    if (verbosity) printf("Num Pages: %u - Num Parallel CRCs: %u - Num blocks = %u\n", num_pages, per, num_blocks);
    /* one input and one output buffer per block (OpenCL dev_input[i]/dev_output[i]), sized for the pages that exist */
    for (unsigned int i = 0; i < num_blocks; i++) {
      size_t off = (size_t)i * per, cnt = (i == num_blocks - 1) ? num_pages_last_block : per;
      unsigned int *in = &h_num[off * num_words], *out = &remainders[off];
#pragma acc enter data create(in[0 : cnt * num_words], out[0 : cnt])
    }
    for (unsigned int k = 0; k < num_wg_sizes; k++) {
      for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
        dwarf_region_begin();
        for (unsigned int i = 0; i < num_blocks; i++) {
          global_size = (i == num_blocks - 1) ? num_pages_last_block : per;
          crc_block(&h_num[(size_t)i * per * num_words], &remainders[(size_t)i * per], global_size, page_size, num_words);
        }
#pragma acc wait
        dwarf_region_end(rep);
        if (verbosity >= 2) /* same lines as the OpenCL host (the first num_blocks remainders) */
          for (unsigned int i = 0; i < num_blocks; i++) printf("Parallel Computation: '%X'\n", remainders[i]);
        if (run_serial) {
          printf("Validating results with serial CRC...\n");
          for (unsigned int i = 0; i < num_pages; i++) {
            unsigned int cpu_remainder = crc32_8bytes(&h_num[(size_t)i * num_words], page_size);
            if (cpu_remainder != remainders[i]) {
              bad++;
              if (bad <= 10)
                fprintf(stderr, "ERROR: OpenACC and CPU Slice-by-8 remainders for page %u differ [OpenACC: '%X', CPU: '%X']\n", i + 1, remainders[i], cpu_remainder);
            }
          }
        }
      }
    }
    for (unsigned int i = 0; i < num_blocks; i++) {
      size_t off = (size_t)i * per, cnt = (i == num_blocks - 1) ? num_pages_last_block : per;
      unsigned int *in = &h_num[off * num_words], *out = &remainders[off];
#pragma acc exit data delete(in[0 : cnt * num_words], out[0 : cnt])
    }
  }
  if (run_serial) printf("%s: %u pages differ from serial slice-by-8\n", bad ? "FAILURE" : "SUCCESS", bad);
  free(h_num); free(remainders);
  return dwarf_report(bad == 0);
}
