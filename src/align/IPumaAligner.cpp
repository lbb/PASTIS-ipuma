// Created by Saliya Ekanayake on 2019-09-03.

#include "../../inc/align/IPumaAligner.hpp"


#include <iostream>
#include <fstream>


#include <algorithm>
#include <utility>

#include <set>

#include "../../inc/util.hpp"

using std::max;
using std::min;
using std::string;
using std::to_string;
using std::vector;

extern shared_ptr<pastis::ParallelOps> parops;

namespace pastis {

void IPumaAligner::construct_seqs(std::shared_ptr<DistFastaData> dfd) {
  uint64_t rseq_cnt = dfd->get_nrseqs();
  uint64_t cseq_cnt = dfd->get_ncseqs();

  rseqs_.resize(rseq_cnt);
  cseqs_.resize(cseq_cnt);

  int cur_nbr = 0;
  uint64_t cur_rseq = 0;
  uint64_t cur_cseq = 0;
  for (auto &nbr : dfd->my_nbrs) {
    uint64_t nbr_seqs_count = nbr.nbr_seq_end_idx - nbr.nbr_seq_start_idx + 1;
    FastaData *curfd = dfd->fd;
    uint64_t sidx = nbr.nbr_seq_start_idx;
    if (nbr.nbr_rank != parops->g_rank)  // not myself
    {
      curfd = dfd->recv_fds_[cur_nbr];
      sidx = 0;
      ++cur_nbr;
    }

    vector<string> &cur_seqs = (nbr.rc_flag == 1 ? rseqs_ : cseqs_);
    uint64_t seq_beg = (nbr.rc_flag == 1 ? cur_rseq : cur_cseq);

#pragma omp parallel
    {
      ushort len;
      uint64_t start_offset, end_offset_inclusive;

#pragma omp for
      for (uint64_t i = 0; i < nbr_seqs_count; ++i) {
        char *buff = curfd->get_sequence(sidx + i, len, start_offset, end_offset_inclusive);
        cur_seqs[seq_beg + i] = std::move(string(buff + start_offset, len));
      }
    }

    if (nbr.rc_flag == 1)
      cur_rseq += nbr_seqs_count;
    else
      cur_cseq += nbr_seqs_count;
  }

  // diag cell row and col seqs are the same
  if (parops->grid->GetRankInProcRow() == parops->grid->GetRankInProcCol()) {
    assert(cseqs_.empty());
    cseqs_.assign(rseqs_.begin(), rseqs_.end());
  }
}

std::tuple<int16_t, int16_t> convertPackedRange(int32_t packed) {
  int16_t begin = packed & 0xffff;
  int16_t end = packed >> 16;
  return {begin, end};
}

void IPumaAligner::aln_batch(std::tuple<uint64_t, uint64_t, CommonKmerLight *> *mattuples, uint64_t beg, uint64_t end, uint64_t bl_roffset, uint64_t bl_coffset, const params_t &params) {
  parops->info("IPUMA Aligner aln_batch");

  parops->tp->start_timer("sim:align_pre");

  uint64_t npairs = end - beg;
  vector<string> seqs_q(npairs);  // queries - shorter seqs
  vector<string> seqs_r(npairs);  // refs - longer seqs
  uint64_t max_rlen = 0;
  uint64_t max_qlen = 0;

  int numThreads = 1;
#ifdef THREADED
  #pragma omp parallel
  {
    numThreads = omp_get_num_threads();
  }
#endif

#pragma omp parallel for reduction(max : max_rlen, max_qlen)
  for (uint64_t i = beg; i < end; ++i) {
    uint64_t lr = std::get<0>(mattuples[i]);
    uint64_t lc = std::get<1>(mattuples[i]);
    assert(lr + bl_roffset < rseqs_.size() && lc + bl_coffset < cseqs_.size());
    string &rseq = rseqs_[lr + bl_roffset];
    string &cseq = cseqs_[lc + bl_coffset];
    if (rseq.size() < cseq.size()) {
      seqs_q[i - beg] = rseq;
      seqs_r[i - beg] = cseq;
    } else {
      seqs_q[i - beg] = cseq;
      seqs_r[i - beg] = rseq;
    }

    max_rlen = max(max_rlen, seqs_r[i - beg].size());
    max_qlen = max(max_qlen, seqs_q[i - beg].size());
  }

  std::vector<ipu::Comparison> pairs(npairs);
  vector<string> seq_pool;
  std:unordered_map<int, int> seqid_map;
  int inserts = 0;

  for (uint64_t i = beg; i < end; ++i) {
    int j = i - beg;
    uint64_t lr = std::get<0>(mattuples[i]);
    uint64_t lc = std::get<1>(mattuples[i]);
    assert(lr + bl_roffset < rseqs_.size() && lc + bl_coffset < cseqs_.size());
    int lrI = 0;
    auto iter = seqid_map.find(lr);
    if (iter == seqid_map.end()) {
      seq_pool.emplace_back(rseqs_[lr + bl_roffset]);
      auto pos = inserts++;
      seqid_map[lr] = pos;
      lrI = pos;
    } else {
       lrI = iter->second;
    }

    int lcI = 0;
     iter = seqid_map.find(lc);
    if (iter == seqid_map.end()) {
      seq_pool.emplace_back(cseqs_[lc + bl_roffset]);
      auto pos = inserts++;
      seqid_map[lc] = pos;
      lcI = pos;
    } else {
       lcI = iter->second;
    }
    pairs[j] = {lrI, lcI};

    max_rlen = max(max_rlen, rseqs_[lr + bl_roffset].size());
    max_qlen = max(max_qlen, cseqs_[lc + bl_roffset].size());
  }

  printf("Map has %d entries, while we have N-pairs %d\n", seqid_map.size(), npairs);

  parops->tp->stop_timer("sim:align_pre");

  parops->tp->start_timer("sim:align");


// {
//   ofstream myfile;
//   myfile.open ("As.txt");
//   for (auto &&i : seqs_r) {
//     myfile << i << '\n';
//   }
//   myfile.close();
// }
// {
//   ofstream myfile;
//   myfile.open ("Bs.txt");
//   for (auto &&i : seqs_q) {
//     myfile << i << '\n';
//   }
//   myfile.close();
// }

  driver_algo->compare_mn_local(seq_pool, pairs);
  // driver_algo->compare_local(seqs_r, seqs_q);
  auto [scores, r_range_result, q_range_result] = driver_algo->get_result();

  parops->tp->stop_timer("sim:align");

  // This seems to be the same block ignore for now
  omp_set_num_threads(numThreads);

  parops->tp->start_timer("sim:align_post");


#pragma omp for
    for (uint64_t i = 0; i < seqs_r.size(); ++i) {
      int len_r = seqs_r[i].size();
      int len_q = seqs_q[i].size();

      auto [ref_begin, ref_end] = convertPackedRange(r_range_result[i]);
      auto [query_begin, query_end] = convertPackedRange(q_range_result[i]);

      int alen_r = ref_end - ref_begin;
      int alen_q = query_end - query_begin;

      double cov_r = (double)(alen_r) / len_r;
      double cov_q = (double)(alen_q) / len_q;

      if (max(cov_r, cov_q) >= params.aln_cov_thr)  // coverage constraint
      {
        CommonKmerLight *ckl = std::get<2>(mattuples[beg + i]);
        ckl->score_aln = (float)(scores[i]) /
                         (float)min(len_r, len_q);
      }
    }

  parops->tp->stop_timer("sim:align_post");
}

void IPumaAligner::construct_seqs_bl(std::shared_ptr<DistFastaData> dfd) {
  parops->logger->log("construct_seqs_bl not implemented");
  exit(1);
}

}  // namespace pastis
