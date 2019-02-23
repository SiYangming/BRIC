#ifndef BLOCK_H
#define BLOCK_H

#include <algorithm>    // std::find
#include <vector>
#include <cstddef> // std::size_t

/* biclustering block */
struct BlockBase {
  std::vector<int> genes;
  std::vector<int> conds;
  double score;
  int block_rows_pre;
  int core_rownum;
  int core_colnum;
  
  std::size_t block_rows() const {
    return genes.size();
  }
  
  std::size_t block_cols() const {
    return conds.size();
  }
  
  bool contains(int gene) const {
    // https://stackoverflow.com/a/3451045 If searching for an element is important, I'd recommend std::set instead of std::vector.
    return std::find(genes.begin(), genes.end(), gene) != genes.end();
  }
  int cond_low_bound;
  long double pvalue;
};

/* KL biclustering block */
struct Block : BlockBase {
  double significance;
};

/* biclustering block in version 1*/
struct Block1 : BlockBase {};

#endif
