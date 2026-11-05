/* Main file for CRC (combinational-logic dwarf) -- OpenMP target-offload translation of
 * opencl/combinational-logic/crc/src/crc_algo.c.
 *
 * Computes a 32-bit Ethernet CRC on a number of input pages with the "Slice-by-8" algorithm
 * (Stephan Brumme / Intel). Same input reader (read_crc), same blocking of pages into blocks of
 * `-p` parallel CRCs (default 16), one work-item per page, same kernel body as crc_kernel.cl, same
 * serial slice-by-8 check (-a). CLI: -i file [-a] [-v] [-p n] [-r n] [-w n] [-k ignored].
 * Timed region: for every block, H2D of its pages, the kernel, D2H of its CRCs.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "ocd_compat.h"
#include "../inc/crc_formats.h"
#include "dwarf_timer.h"

/* the lookup table must exist on the device as well (OpenCL: __constant) */
#pragma omp declare target
#include "../inc/eth_crc32_lut.h"
#pragma omp end declare target

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
#pragma omp declare target
static inline unsigned int crc_item(const unsigned int *data, unsigned int length_bytes, unsigned int length_ints,
                                    size_t gid) {
  unsigned int crc = 0xFFFFFFFF, one, two;
  size_t i = gid * length_ints, j;
  while (length_bytes >= 8) {  // process eight bytes at once
    one = data[i++] ^ crc;
    two = data[i++];
    crc = crc32Lookup[7][one & 0xFF] ^ crc32Lookup[6][(one >> 8) & 0xFF] ^ crc32Lookup[5][(one >> 16) & 0xFF] ^
          crc32Lookup[4][one >> 24] ^ crc32Lookup[3][two & 0xFF] ^ crc32Lookup[2][(two >> 8) & 0xFF] ^
          crc32Lookup[1][(two >> 16) & 0xFF] ^ crc32Lookup[0][two >> 24];
    length_bytes -= 8;
  }
  while (length_bytes) {  // remaining 1 to 7 bytes
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
#pragma omp end declare target

/* enqueueCRCDevice(): write the block's pages, run the kernel over global_size work-items, read the
 * CRCs back. crc_kernel.cl uses neither barriers nor local memory, so the work-group partition
 * (local_size, -w) has no semantic effect and the kernel is expressed as one flat parallel loop. */
static void crc_block(const unsigned int *h_num, unsigned int *h_answer, size_t global_size, size_t local_size,
                      unsigned int page_size, unsigned int num_words) {
  size_t nbytes_words = (size_t)num_words * global_size;
  (void)local_size;
#pragma omp target teams distribute parallel for map(to : h_num[0 : nbytes_words]) map(from : h_answer[0 : global_size])
  for (size_t gid = 0; gid < global_size; gid++) h_answer[gid] = crc_item(h_num, page_size, num_words, gid);
}

static void usage(void) {
  printf("crc -i <input_file> [-a] [-v] [-p <parallel crcs>] [-r <num_execs>] [-w <wg_size>]\n");
  exit(0);
}

int main(int argc, char **argv) {
  size_t global_size, local_size;
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
      case 'k': printf("Note: -k ignored (single OpenMP kernel)\n"); break;
      default: fprintf(stderr, "Invalid argument\n\n"); usage();
    }
  }
  check(file != NULL, "-i option must be supplied!");
  h_num = read_crc(&num_pages, &page_size, file);
  if (!num_block_sizes) { num_block_sizes = 1; num_parallel_crcs = malloc(sizeof(unsigned int)); num_parallel_crcs[0] = 16; }
  if (!num_wg_sizes) { num_wg_sizes = 1; wg_sizes = malloc(sizeof(size_t)); wg_sizes[0] = 1; }
  unsigned int num_words = page_size / 4;
  if (verbosity) printf("num_words = %u\n", num_words);

  unsigned int *remainders = (unsigned int *)malloc(sizeof(unsigned int) * num_pages);
  unsigned int bad = 0;
  dwarf_init("crc", "combinational-logic", "omp");
  dwarf_set_size(num_pages, page_size);

  for (unsigned int h = 0; h < num_block_sizes; h++) {
    unsigned int num_blocks = num_pages / num_parallel_crcs[h], num_pages_last_block;
    if (num_pages % num_parallel_crcs[h] != 0) { num_blocks++; num_pages_last_block = num_pages % num_parallel_crcs[h]; }
    else num_pages_last_block = num_parallel_crcs[h];
    if (verbosity) printf("Num Pages: %u - Num Parallel CRCs: %u - Num blocks = %u\n", num_pages, num_parallel_crcs[h], num_blocks);
    for (unsigned int k = 0; k < num_wg_sizes; k++) {
      if (num_execs > 1) dwarf_set_repeats((int)num_execs);
      for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
        dwarf_region_begin();
        for (unsigned int i = 0; i < num_blocks; i++) {
          if (i == num_blocks - 1) {  // last block
            global_size = num_pages_last_block;
            local_size = wg_sizes[k];
            if ((global_size % local_size) != 0) {  // same fix-up as the OpenCL host code
              local_size = 1;
              while ((global_size % local_size) == 0) local_size = local_size << 1;
              local_size = local_size >> 1;
            }
          } else {
            global_size = num_parallel_crcs[h];
            local_size = wg_sizes[k];
          }
          crc_block(&h_num[(size_t)i * num_parallel_crcs[h] * num_words], &remainders[(size_t)i * num_parallel_crcs[h]],
                    global_size, local_size, page_size, num_words);
        }
        dwarf_region_end(rep);
        if (run_serial) {
          printf("Validating results with serial CRC...\n");
          for (unsigned int i = 0; i < num_pages; i++) {
            unsigned int cpu_remainder = crc32_8bytes(&h_num[(size_t)i * num_words], page_size);
            if (cpu_remainder != remainders[i]) {
              bad++;
              if (bad <= 10)
                fprintf(stderr, "ERROR: OMP and CPU Slice-by-8 remainders for page %u differ [OMP: '%X', CPU: '%X']\n", i + 1, remainders[i], cpu_remainder);
            }
          }
        }
      }
    }
  }
  if (run_serial) printf("%s: %u pages differ from serial slice-by-8\n", bad ? "FAILURE" : "SUCCESS", bad);
  free(h_num); free(remainders);
  return dwarf_report(bad == 0);
}
