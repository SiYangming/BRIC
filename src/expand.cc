/************************************************************************/
/* Author: Qin Ma <maqin@uga.edu>, Step. 19, 2013
 * Biclustering expansion, greedy add the possible genes and negatively
 * regulated genes from outside of the current bicluster in a given background
 */

#include "expand.h"
#include "block.h"
#include "cluster.h"
#include "read_array.h"
#include "write_block.h"

discrete **another_arr_c;
char **another_genes;
char **another_conds;
int another_rows;
int another_cols;

/**********************************************************************************/
static char *atom = NULL;
static char delims[] = " \t\r\n";

static int intersect_row(const std::vector<discrete> &colcand, discrete *g2,
                         const int cols) {
  int cnt = 0;
  for (int col = 0; col < cols; col++)
    if (colcand[col] != 0 && colcand[col] == g2[col])
      cnt++;
  return cnt;
}
static int reverse_row(const std::vector<discrete> &colcand, discrete *g2,
                       const int cols) {
  int cnt = 0;
  for (int col = 0; col < cols; col++) {
    if (colcand[col] != 0 && symbols[colcand[col]] == -symbols[g2[col]])
      cnt++;
  }
  return cnt;
}

void store_block(std::unique_ptr<Block> &b_ptr, const std::vector<int> &ge,
                 const std::vector<int> &co) {
  b_ptr->genes = ge;
  b_ptr->conds = co;
}
static void init_expand() {
  another_genes = genes_n;
  another_conds = conds_n;
  another_arr_c = arr_c;
  another_rows = rows;
  another_cols = cols;
}

