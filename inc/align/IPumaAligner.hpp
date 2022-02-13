#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "PWAlign.hpp"
#include "../../ipuma-lib/src/driver.hpp"
#include "../../ipuma-lib/src/swatlib/vector.hpp"

namespace pastis {

// uses light representation
class
    IPumaAligner : public PWAlign {
 private:
  std::vector<std::string> rseqs_;
  std::vector<std::string> cseqs_;
  std::tuple<int, int> gaps_;
  ipu::batchaffine::SWAlgorithm *driver_algo;
  uint32_t g_batch_sz_;

 public:
  IPumaAligner(int gap_open, int gap_ext, uint32_t bsz = 1e6) : PWAlign(), g_batch_sz_(bsz) {
    this->gaps_ = std::make_tuple<>(gap_open, gap_ext);

        static const ipu::SWConfig SW_CONFIGURATION = {
            .gapInit = -gap_open,
            .gapExtend = -gap_ext,
            .matchValue = ALN_MATCH_SCORE,
            .mismatchValue = -ALN_MISMATCH_COST,
            .ambiguityValue = -ALN_AMBIGUITY_COST,
            .similarity = swatlib::Similarity::blosum62,
            .datatype = swatlib::DataType::aminoAcid,
    };

static const ipu::batchaffine::IPUAlgoConfig ALGO_CONFIGURATION = {
            KLIGN_IPU_TILES,
            KLIGN_IPU_MAXAB_SIZE,
            KLIGN_IPU_MAX_BATCHES,
            KLIGN_IPU_BUFSIZE,
            ipu::batchaffine::VertexType::cpp,
            ipu::partition::Algorithm::roundRobin
};

    init_single_ipu(SW_CONFIGURATION, ALGO_CONFIGURATION);
    this->driver_algo = getDriver();
  }

  ~IPumaAligner() {
  }

  void
  construct_seqs(std::shared_ptr<DistFastaData> dfd) override;

  void
  construct_seqs_bl(std::shared_ptr<DistFastaData> dfd) override;

  void
  aln_batch(std::tuple<uint64_t, uint64_t, CommonKmerLight *> *mattuples,
            uint64_t beg, uint64_t end,
            uint64_t bl_roffset, uint64_t bl_coffset,
            const params_t &params) override;

  size_t
  rseq_len(uint64_t i) {
    return rseqs_[i].size();
  }

  size_t
  cseq_len(uint64_t i) {
    return cseqs_[i].size();
  }
};

}  // namespace pastis
