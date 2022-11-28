/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Wavefront Alignment Algorithms benchmarking tool
 */
#include <omp.h>

#define EXTERNAL_BENCHMARKS

#include "utils/commons.h"
#include "utils/sequence_buffer.h"
#include "system/profiler_timer.h"

#include "alignment/score_matrix.h"
#include "edit/edit_dp.h"
#include "gap_linear/nw.h"
#include "gap_affine/swg.h"
#include "gap_affine2p/affine2p_matrix.h"
#include "gap_affine2p/affine2p_dp.h"
#include "wavefront/wavefront_align.h"

#include "benchmark/benchmark_indel.h"
#include "benchmark/benchmark_edit.h"
#include "benchmark/benchmark_gap_linear.h"
#include "benchmark/benchmark_gap_affine.h"
#include "benchmark/benchmark_gap_affine2p.h"

#ifdef EXTERNAL_BENCHMARKS
#include "benchmark/external/benchmark_bitpal.h"
#include "benchmark/external/benchmark_daligner.h"
#include "benchmark/external/benchmark_diffutils.h"
#include "benchmark/external/benchmark_edlib.h"
#include "benchmark/external/benchmark_gaba.h"
#include "benchmark/external/benchmark_ksw2.h"
#include "benchmark/external/benchmark_lh3.h"
#include "benchmark/external/benchmark_wfalm.h"
#endif

/*
 * Algorithms
 */
typedef enum {
  // Test
  alignment_test,
  // Indel
  alignment_indel_wavefront,
  // Edit
  alignment_edit_bpm,
  alignment_edit_dp,
  alignment_edit_dp_banded,
  alignment_edit_wavefront,
  // Gap-linear
  alignment_gap_linear_nw,
  alignment_gap_linear_wavefront,
  // Gap-affine
  alignment_gap_affine_swg,
  alignment_gap_affine_swg_endsfree,
  alignment_gap_affine_swg_banded,
  alignment_gap_affine_wavefront,
  // Gap-affine dual-cost
  alignment_gap_affine2p_dp,
  alignment_gap_affine2p_wavefront,
#ifdef EXTERNAL_BENCHMARKS
  // External algorithms
  alignment_bitpal_edit,
  alignment_bitpal_scored,
  alignment_daligner,
  alignment_diffutils,
  alignment_edlib,
  alignment_gaba_aband,
  alignment_ksw2_extz2_sse,
  alignment_ksw2_extd2_sse,
  alignment_lv89,
  alignment_wfalm,
  alignment_wfalm_lowmem,
  alignment_wfalm_rec
#endif
} alignment_algorithm_type;
bool align_benchmark_is_wavefront(
    const alignment_algorithm_type algorithm) {
  return algorithm == alignment_indel_wavefront ||
         algorithm == alignment_edit_wavefront ||
         algorithm == alignment_gap_linear_wavefront ||
         algorithm == alignment_gap_affine_wavefront ||
         algorithm == alignment_gap_affine2p_wavefront;
}
/*
 * Generic parameters
 */
typedef struct {
  // Algorithm
  alignment_algorithm_type algorithm;
  // I/O
  char *input_filename;
  char *output_filename;
  bool output_full;
  int target_bases_aligned;
  // I/O internals
  FILE* input_file;
  char* line1;
  char* line2;
  size_t line1_allocated;
  size_t line2_allocated;
  FILE* output_file;
  // Penalties
  linear_penalties_t linear_penalties;
  affine_penalties_t affine_penalties;
  affine2p_penalties_t affine2p_penalties;
  // Wavefront parameters
  bool wfa_score_only;
  wf_heuristic_strategy wfa_heuristic;
  int wfa_heuristic_p1;
  int wfa_heuristic_p2;
  int wfa_heuristic_p3;
  wavefront_memory_t wfa_memory_mode;
  uint64_t wfa_max_memory;
  int wfa_max_threads;
  // Other algorithms parameters
  int bandwidth;
#ifdef EXTERNAL_BENCHMARKS
  bool ksw2_approx_max__drop;
  int ksw2_bandwidth;
  int ksw2_zdrop;
#endif
  // Misc
  bool check_display;
  bool check_correct;
  bool check_score;
  bool check_alignments;
  int check_metric;
  int check_bandwidth;
  int plot;
  // Profile
  profiler_timer_t timer_global;
  // System
  int num_threads;
  int batch_size;
  int progress;
  int verbose;
} benchmark_args;
benchmark_args parameters = {
  // Algorithm
  .algorithm = alignment_test,
  // I/O
  .input_filename = NULL,
  .output_filename = NULL,
  .output_full = false,
  .output_file = NULL,
  .target_bases_aligned = 0,
  // I/O internals
  .input_file = NULL,
  .line1 = NULL,
  .line2 = NULL,
  .line1_allocated = 0,
  .line2_allocated = 0,
  // Penalties
  .linear_penalties = {
      .match = 0,
      .mismatch = 4,
      .indel = 2,
  },
  .affine_penalties = {
      .match = 0,
      .mismatch = 4,
      .gap_opening = 6,
      .gap_extension = 2,
  },
  .affine2p_penalties = {
      .match = 0,
      .mismatch = 4,
      .gap_opening1 = 6,
      .gap_extension1 = 2,
      .gap_opening2 = 24,
      .gap_extension2 = 1,
  },
  // Wavefront parameters
  .wfa_score_only = false,
  .wfa_heuristic = wf_heuristic_none,
  .wfa_heuristic_p1 = -1,
  .wfa_heuristic_p2 = -1,
  .wfa_heuristic_p3 = -1,
  .wfa_memory_mode = wavefront_memory_ultralow,
  .wfa_max_memory = UINT64_MAX,
  .wfa_max_threads = 1,
  // Other algorithms parameters
  .bandwidth = -1,
#ifdef EXTERNAL_BENCHMARKS
  .ksw2_approx_max__drop = false,
  .ksw2_bandwidth = -1,
  .ksw2_zdrop = -1,
#endif
  // Misc
  .check_bandwidth = -1,
  .check_display = false,
  .check_correct = false,
  .check_score = false,
  .check_alignments = false,
  .check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE,
  .plot = 0,
  // System
  .num_threads = 1,
  .batch_size = 10000,
  .progress = 100000,
  .verbose = 0,
};

/*
 * Benchmark UTest
 */
void align_pairwise_test() { }

/*
 * Configuration
 */
wavefront_aligner_t* align_input_configure_wavefront(
    align_input_t* const align_input) {
  // Set attributes
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.memory_mode = parameters.wfa_memory_mode;
  if (parameters.wfa_score_only) {
    attributes.alignment_scope = compute_score;
  }
  // WF-Heuristic
  switch (parameters.wfa_heuristic) {
    case wf_heuristic_none:
      attributes.heuristic.strategy = wf_heuristic_none;
      break;
    case wf_heuristic_banded_static:
      attributes.heuristic.strategy = wf_heuristic_banded_static;
      attributes.heuristic.min_k = parameters.wfa_heuristic_p1;
      attributes.heuristic.max_k = parameters.wfa_heuristic_p2;
      break;
    case wf_heuristic_banded_adaptive:
      attributes.heuristic.strategy = wf_heuristic_banded_adaptive;
      attributes.heuristic.min_k = parameters.wfa_heuristic_p1;
      attributes.heuristic.max_k = parameters.wfa_heuristic_p2;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p3;
      break;
    case wf_heuristic_wfadaptive:
      attributes.heuristic.strategy = wf_heuristic_wfadaptive;
      attributes.heuristic.min_wavefront_length = parameters.wfa_heuristic_p1;
      attributes.heuristic.max_distance_threshold = parameters.wfa_heuristic_p2;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p3;
      break;
    case wf_heuristic_xdrop:
      attributes.heuristic.strategy = wf_heuristic_xdrop;
      attributes.heuristic.xdrop = parameters.wfa_heuristic_p1;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p2;
      break;
    case wf_heuristic_zdrop:
      attributes.heuristic.strategy = wf_heuristic_zdrop;
      attributes.heuristic.zdrop = parameters.wfa_heuristic_p1;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p2;
      break;
    default:
      break;
  }
  // Select flavor
  switch (parameters.algorithm) {
    case alignment_indel_wavefront:
      attributes.distance_metric = indel;
      break;
    case alignment_edit_wavefront:
      attributes.distance_metric = edit;
      break;
    case alignment_gap_linear_wavefront:
      attributes.distance_metric = gap_linear;
      attributes.linear_penalties = parameters.linear_penalties;
      break;
    case alignment_gap_affine_wavefront:
      attributes.distance_metric = gap_affine;
      attributes.affine_penalties = parameters.affine_penalties;
      break;
    case alignment_gap_affine2p_wavefront:
      attributes.distance_metric = gap_affine_2p;
      attributes.affine2p_penalties = parameters.affine2p_penalties;
      break;
    default:
      return NULL; // No WF selected
      break;
  }
  // Select alignment form
  attributes.alignment_form.span = alignment_end2end;
  // Misc
  attributes.plot.enabled = (parameters.plot > 0);
  attributes.plot.resolution_points = parameters.plot;
  attributes.system.verbose = parameters.verbose;
  attributes.system.max_memory_abort = parameters.wfa_max_memory;
  attributes.system.max_num_threads = parameters.wfa_max_threads;
  // Allocate
  return wavefront_aligner_new(&attributes);
}
void align_input_configure_global(
    align_input_t* const align_input) {
  // Clear
  benchmark_align_input_clear(align_input);
  // Penalties
  align_input->linear_penalties = parameters.linear_penalties;
  align_input->affine_penalties = parameters.affine_penalties;
  align_input->affine2p_penalties = parameters.affine2p_penalties;
  // Alignment form
  align_input->ends_free = false;
  // Output
  align_input->output_file = parameters.output_file;
  align_input->output_full = parameters.output_full;
  // MM
  align_input->mm_allocator = mm_allocator_new(BUFFER_SIZE_1M);
  // WFA
  if (align_benchmark_is_wavefront(parameters.algorithm)) {
    align_input->wf_aligner = align_input_configure_wavefront(align_input);
  } else {
    align_input->wf_aligner = NULL;
  }
  // PROFILE/STATS
  timer_reset(&align_input->timer);
  // DEBUG
  align_input->debug_flags = 0;
  align_input->debug_flags |= parameters.check_metric;
  if (parameters.check_display) align_input->debug_flags |= ALIGN_DEBUG_DISPLAY_INFO;
  if (parameters.check_correct) align_input->debug_flags |= ALIGN_DEBUG_CHECK_CORRECT;
  if (parameters.check_score) align_input->debug_flags |= ALIGN_DEBUG_CHECK_SCORE;
  if (parameters.check_alignments) align_input->debug_flags |= ALIGN_DEBUG_CHECK_ALIGNMENT;
  align_input->check_linear_penalties = &parameters.linear_penalties;
  align_input->check_affine_penalties = &parameters.affine_penalties;
  align_input->check_affine2p_penalties = &parameters.affine2p_penalties;
  align_input->check_bandwidth = parameters.check_bandwidth;
  align_input->verbose = parameters.verbose;
}
void align_benchmark_free(
    align_input_t* const align_input) {
  if (align_input->wf_aligner) wavefront_aligner_delete(align_input->wf_aligner);
  mm_allocator_delete(align_input->mm_allocator);
}
/*
 * I/O
 */
bool align_benchmark_read_input(
    FILE* input_file,
    char** line1,
    char** line2,
    size_t* line1_allocated,
    size_t* line2_allocated,
    const int seqs_processed,
    align_input_t* const align_input) {
  // Parameters
  int line1_length=0, line2_length=0;
  // Read queries
  line1_length = getline(line1,line1_allocated,input_file);
  if (line1_length==-1) return false;
  line2_length = getline(line2,line2_allocated,input_file);
  if (line1_length==-1) return false;
  // Configure input
  align_input->sequence_id = seqs_processed;
  align_input->pattern = *line1 + 1;
  align_input->pattern_length = line1_length - 2;
  align_input->pattern[align_input->pattern_length] = '\0';
  align_input->text = *line2 + 1;
  align_input->text_length = line2_length - 2;
  align_input->text[align_input->text_length] = '\0';
  return true;
}
/*
 * Display
 */
void align_benchmark_print_progress(
    const int seqs_processed) {
  const uint64_t time_elapsed_alg = timer_get_current_total_ns(&parameters.timer_global);
  const float rate_alg = (float)seqs_processed/(float)TIMER_CONVERT_NS_TO_S(time_elapsed_alg);
  fprintf(stderr,"...processed %d reads (alignment = %2.3f seq/s)\n",seqs_processed,rate_alg);
}
void align_benchmark_print_results(
    align_input_t* const align_input,
    const int seqs_processed,
    const bool print_stats) {
  // Print benchmark results
  fprintf(stderr,"[Benchmark]\n");
  fprintf(stderr,"=> Total.reads            %d\n",seqs_processed);
  fprintf(stderr,"=> Time.Benchmark      ");
  timer_print(stderr,&parameters.timer_global,NULL);
  fprintf(stderr,"  => Time.Alignment    %2.3f (s)\n",
      TIMER_CONVERT_NS_TO_S(timer_get_total_ns(&align_input->timer)));
  //timer_print(stderr,&align_input->timer,&parameters.timer_global);
  // Print Stats
  const bool checks_enabled =
      parameters.check_display || parameters.check_correct ||
      parameters.check_score || parameters.check_alignments;
  if (checks_enabled) {
    benchmark_print_stats(stderr,align_input,true);
  }
}
/*
 * Benchmark
 */
void align_benchmark_run_algorithm(
    align_input_t* const align_input) {
  // Select algorithm
  switch (parameters.algorithm) {
    // Indel
    case alignment_indel_wavefront:
      benchmark_indel_wavefront(align_input);
      break;
    // Edit
    case alignment_edit_bpm:
      benchmark_edit_bpm(align_input);
      break;
    case alignment_edit_dp:
      benchmark_edit_dp(align_input);
      break;
    case alignment_edit_dp_banded:
      benchmark_edit_dp_banded(align_input,parameters.bandwidth);
      break;
    case alignment_edit_wavefront:
      benchmark_edit_wavefront(align_input);
      break;
    // Gap-linear
    case alignment_gap_linear_nw:
      benchmark_gap_linear_nw(align_input,&parameters.linear_penalties);
      break;
    case alignment_gap_linear_wavefront:
      benchmark_gap_linear_wavefront(align_input,&parameters.linear_penalties);
      break;
    // Gap-affine
    case alignment_gap_affine_swg:
      benchmark_gap_affine_swg(align_input,&parameters.affine_penalties);
      break;
    case alignment_gap_affine_swg_endsfree:
      benchmark_gap_affine_swg_endsfree(
          align_input,&parameters.affine_penalties);
      break;
    case alignment_gap_affine_swg_banded:
      benchmark_gap_affine_swg_banded(align_input,
          &parameters.affine_penalties,parameters.bandwidth);
      break;
    case alignment_gap_affine_wavefront:
      benchmark_gap_affine_wavefront(align_input,&parameters.affine_penalties);
      break;
    // Gap-affine 2p
    case alignment_gap_affine2p_dp:
      benchmark_gap_affine2p_dp(align_input,&parameters.affine2p_penalties);
      break;
    case alignment_gap_affine2p_wavefront:
      benchmark_gap_affine2p_wavefront(align_input,&parameters.affine2p_penalties);
      break;
#ifdef EXTERNAL_BENCHMARKS
      /*
       * External Algorithms
       */
      case alignment_bitpal_edit:
        benchmark_bitpal_m0_x1_g1(align_input);
        break;
      case alignment_bitpal_scored:
        benchmark_bitpal_m1_x4_g2(align_input);
        break;
      case alignment_daligner:
        benchmark_daligner(align_input);
        break;
      case alignment_diffutils:
        benchmark_diffutils(align_input,true);
        break;
      case alignment_edlib:
        benchmark_edlib(align_input);
        break;
      case alignment_gaba_aband:
        benchmark_gaba_aband(align_input,&parameters.affine_penalties);
        break;
      case alignment_ksw2_extz2_sse:
        benchmark_ksw2_extz2_sse(
            align_input,&parameters.affine_penalties,
            parameters.ksw2_approx_max__drop,
            parameters.ksw2_bandwidth,parameters.ksw2_zdrop);
        break;
      case alignment_ksw2_extd2_sse:
        benchmark_ksw2_extd2_sse(
            align_input,&parameters.affine2p_penalties,
            parameters.ksw2_approx_max__drop,
            parameters.ksw2_bandwidth,parameters.ksw2_zdrop);
        break;
      case alignment_lv89:
        benchmark_lv89(align_input);
        break;
      case alignment_wfalm:
        benchmark_wfalm_global_affine(align_input,&parameters.affine_penalties);
        break;
      case alignment_wfalm_lowmem:
        benchmark_wfalm_global_affine_lowmem(align_input,&parameters.affine_penalties);
        break;
      case alignment_wfalm_rec:
        benchmark_wfalm_global_affine_rec(align_input,&parameters.affine_penalties);
        break;
#endif
    default:
      fprintf(stderr,"Algorithm not implemented\n");
      exit(1);
      break;
  }
}
void align_benchmark_sequential() {
  // PROFILE
  timer_reset(&parameters.timer_global);
  timer_start(&parameters.timer_global);
  // I/O files
  parameters.input_file = fopen(parameters.input_filename, "r");
  if (parameters.input_file == NULL) {
    fprintf(stderr,"Input file '%s' couldn't be opened\n",parameters.input_filename);
    exit(1);
  }
  if (parameters.output_filename != NULL) {
    parameters.output_file = fopen(parameters.output_filename, "w");
  }
  // Global configuration
  align_input_t align_input;
  align_input_configure_global(&align_input);
  // Read-align loop
  const bool aggregated_stats = (parameters.target_bases_aligned == 0);
  int seqs_processed = 0, progress = 0;
  while (true) {
    // Read input sequence-pair
    const bool input_read = align_benchmark_read_input(
        parameters.input_file,&parameters.line1,&parameters.line2,
        &parameters.line1_allocated,&parameters.line2_allocated,
        seqs_processed,&align_input);
    if (!input_read) break;
    // Execute the selected algorithm
    if (aggregated_stats) {
      align_benchmark_run_algorithm(&align_input);
    } else {
      timer_reset(&align_input.timer);
      int bases_aligned = 0, times_aligned = 0;
      while (bases_aligned < parameters.target_bases_aligned) {
        align_benchmark_run_algorithm(&align_input);
        bases_aligned += align_input.pattern_length;
        times_aligned++;
      }
      /*
       * Print stats
       * <LENGTH>    <TIME(ms)>    <SCORE>    <ERROR(%)>
       */
      const double time_ms = TIMER_CONVERT_NS_TO_MS(timer_get_total_ns(&align_input.timer));
      if (align_input.wf_aligner == NULL) {
        fprintf(stderr,"%d\t%2.6f\n",
            align_input.pattern_length,
            time_ms/(double)times_aligned);
      } else {
        fprintf(stderr,"%d\t%2.6f\t%d\t%2.2f\n",
            align_input.pattern_length,
            time_ms/(double)times_aligned,
            align_input.wf_aligner->cigar->score,
            100.0*(double)cigar_score_edit(align_input.wf_aligner->cigar)/(double)align_input.pattern_length);
      }
    }
    // Update progress
    ++seqs_processed;
    if (aggregated_stats && ++progress == parameters.progress) {
      progress = 0;
      align_benchmark_print_progress(seqs_processed);
    }
  }
  // Print benchmark results
  if (aggregated_stats) {
    timer_stop(&parameters.timer_global);
    align_benchmark_print_results(&align_input,seqs_processed,true);
  }
  // Free
  align_benchmark_free(&align_input);
  fclose(parameters.input_file);
  if (parameters.output_file) fclose(parameters.output_file);
  free(parameters.line1);
  free(parameters.line2);
}
/*
 * Generic Menu
 */
void usage() {
  fprintf(stderr,
      "USE: ./align_benchmark -a <algorithm> -i <input>                        \n"
      "      Options::                                                         \n"
      "        [Algorithm]                                                     \n"
      "          --algorithm|a <algorithm>                                     \n"
      "            [Indel (Longest Common Subsequence)]                        \n"
      "              indel-wfa                                                 \n"
      "            [Edit (Levenshtein)]                                        \n"
      "              edit-bpm                                                  \n"
      "              edit-dp                                                   \n"
      "              edit-dp-banded                                            \n"
      "              edit-wfa                                                  \n"
      "            [Gap-linear (Needleman-Wunsch)]                             \n"
      "              gap-linear-nw                                             \n"
      "              gap-linear-wfa                                            \n"
      "            [Gap-affine (Smith-Waterman-Gotoh)]                         \n"
      "              gap-affine-swg                                            \n"
      "              gap-affine-swg-banded                                     \n"
      "              gap-affine-wfa                                            \n"
      "            [Gap-affine-2pieces (Concave 2-pieces)]                     \n"
      "              gap-affine2p-dp                                           \n"
      "              gap-affine2p-wfa                                          \n"
#ifdef EXTERNAL_BENCHMARKS
      "            [External/BitPal]                                           \n"
      "              bitpal-edit          (Edit)[score-only]                   \n"
      "              bitpal-scored        (Gap-linear)[score-only]             \n"
      "            [External/Daligner]                                         \n"
      "              daligner             (Edit)                               \n"
      "            [External/Diffutils]                                        \n"
      "              diffutils            (Edit)                               \n"
      "            [External/Edlib]                                            \n"
      "              edlib                (Edit)                               \n"
      "            [External/GABA]                                             \n"
      "              gaba-aband           (Gap-affine)                         \n"
      "            [External/KSW2]                                             \n"
      "              ksw2-extz2-sse       (Gap-affine)                         \n"
      "              ksw2-extd2-sse       (Gap-affine-2pieces)                 \n"
      "            [External/LV89]                                             \n"
      "              lv89                 (Edit)[score-only]                   \n"
      "            [External/WFAlm]                                            \n"
      "              wfalm                (Gap-affine)                         \n"
      "              wfalm-lowmem         (Gap-affine)[low-mem]                \n"
      "              wfalm-rec            (Gap-affine)[ultralow-mem]           \n"
#endif
      "        [Input & Output]                                                \n"
      "          --input|i <File>                                              \n"
      "          --output|o <File>                                             \n"
      "          --output-full <File>                                          \n"
      "          --repeat|-r <TargetBasesAligned>                              \n"
      "        [Penalties & Span]                                              \n"
      "          --linear-penalties|p M,X,I                                    \n"
      "          --affine-penalties|g M,X,O,E                                  \n"
      "          --affine2p-penalties M,X,O1,E1,O2,E2                          \n"
      "          --ends-free P0,Pf,T0,Tf                                       \n"
      "        [Wavefront parameters]                                          \n"
      "          --wfa-score-only                                              \n"
      "          --wfa-memory-mode 'high'|'med'|'low'|'ultralow'               \n"
      "          --wfa-heuristic <Strategy>                                    \n"
      "          --wfa-heuristic-parameters  <P1>,<P2>[,<P3>]                  \n"
      "            [Strategy='banded-static']                                  \n"
      "              P1 = minimum-diagonal-band (e.g., -100)                   \n"
      "              P2 = maximum-diagonal-band (e.g., +100)                   \n"
      "            [Strategy='banded-adaptive']                                \n"
      "              P1 = minimum-diagonal-band (e.g., -100)                   \n"
      "              P2 = maximum-diagonal-band (e.g., +100)                   \n"
      "              P3 = steps-between-cutoffs                                \n"
      "            [Strategy='wfa-adaptive']                                   \n"
      "              P1 = minimum-wavefront-length                             \n"
      "              P2 = maximum-difference-distance                          \n"
      "              P3 = steps-between-cutoffs                                \n"
      "            [Strategy='xdrop']                                          \n"
      "              P1 = x-drop                                               \n"
      "              P2 = steps-between-cutoffs                                \n"
      "            [Strategy='zdrop']                                          \n"
      "              P1 = z-drop                                               \n"
      "              P2 = steps-between-cutoffs                                \n"
      "          --wfa-max-memory <Bytes>                                      \n"
      "          --wfa-bidirectional                                           \n"
      "          --wfa-max-threads <INT> (intra-parallelism; default=1)        \n"
      "        [Other Parameters]                                              \n"
      "          --bandwidth <INT>                                             \n"
#ifdef EXTERNAL_BENCHMARKS
      "          --ba-block-size <INT>                                         \n"
      "          --ksw2-approx-max-drop                                        \n"
      "          --ksw2-bandwidth <INT>                                        \n"
      "          --ksw2-zdrop <INT>                                            \n"
#endif
      "        [Misc]                                                          \n"
      "          --check|c 'correct'|'score'|'alignment'                       \n"
      "          --check-distance 'indel'|'edit'|'linear'|'affine'|'affine2p'  \n"
      "          --check-bandwidth <INT>                                       \n"
      "          --plot                                                        \n"
      "        [System]                                                        \n"
    //"          --num-threads|t <INT>                                         \n"
    //"          --batch-size <INT>                                            \n"
    //"          --progress|P <INT>                                            \n"
      "          --help|h                                                      \n");
}
void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    /* Input */
    { "algorithm", required_argument, 0, 'a' },
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "output-full", required_argument, 0, 800 },
    { "repeat", required_argument, 0, 'r' },
    /* Penalties */
    { "linear-penalties", required_argument, 0, 'p' },
    { "affine-penalties", required_argument, 0, 'g' },
    { "affine2p-penalties", required_argument, 0, 900 },
    /* Wavefront parameters */
    { "wfa-score-only", no_argument, 0, 1000 },
    { "wfa-memory-mode", required_argument, 0, 1001 },
    { "wfa-heuristic", required_argument, 0, 1002 },
    { "wfa-heuristic-parameters", required_argument, 0, 1003 },
    { "wfa-max-memory", required_argument, 0, 1005 },
    { "wfa-max-threads", required_argument, 0, 1007 },
    /* Other alignment parameters */
    { "bandwidth", required_argument, 0, 2000 },
#ifdef EXTERNAL_BENCHMARKS
    { "ba-block-size", required_argument, 0, 2001 },
    { "ksw2-approx-max-drop", no_argument, 0, 2002 },
    { "ksw2-bandwidth", required_argument, 0, 2003 },
    { "ksw2-zdrop", required_argument, 0, 2004 },
#endif
    /* Misc */
    { "check", required_argument, 0, 'c' },
    { "check-distance", required_argument, 0, 3001 },
    { "check-bandwidth", required_argument, 0, 3002 },
    { "plot", optional_argument, 0, 3003 },
    /* System */
    { "num-threads", required_argument, 0, 't' },
    { "batch-size", required_argument, 0, 4000 },
    { "progress", required_argument, 0, 'P' },
    { "verbose", optional_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  if (argc <= 1) {
    usage();
    exit(0);
  }
  while (1) {
    c=getopt_long(argc,argv,"a:i:o:r:p:g:P:c:v::t:h",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    /*
     * Algorithm
     */
    case 'a': {
      /* Test bench */
      if (strcmp(optarg,"test")==0) {
        parameters.algorithm = alignment_test;
      // Indel
      } else if (strcmp(optarg,"indel-wfa")==0) {
        parameters.algorithm = alignment_indel_wavefront;
      // Edit
      } else if (strcmp(optarg,"edit-bpm")==0) {
        parameters.algorithm = alignment_edit_bpm;
      } else if (strcmp(optarg,"edit-dp")==0) {
        parameters.algorithm = alignment_edit_dp;
      } else if (strcmp(optarg,"edit-dp-banded")==0) {
        parameters.algorithm = alignment_edit_dp_banded;
      } else if (strcmp(optarg,"edit-wfa")==0) {
        parameters.algorithm = alignment_edit_wavefront;
      // Gap-Linear
      } else if (strcmp(optarg,"gap-linear-nw")==0 ||
                 strcmp(optarg,"gap-linear-dp")==0) {
        parameters.algorithm = alignment_gap_linear_nw;
      } else if (strcmp(optarg,"gap-linear-wfa")==0) {
        parameters.algorithm = alignment_gap_linear_wavefront;
      // Gap-Affine
      } else if (strcmp(optarg,"gap-affine-swg")==0 ||
                 strcmp(optarg,"gap-affine-dp")==0) {
        parameters.algorithm = alignment_gap_affine_swg;
      } else if (strcmp(optarg,"gap-affine-swg-banded")==0 ||
                 strcmp(optarg,"gap-affine-dp-banded")==0) {
        parameters.algorithm = alignment_gap_affine_swg_banded;
      } else if (strcmp(optarg,"gap-affine-wfa")==0) {
        parameters.algorithm = alignment_gap_affine_wavefront;
      // Gap-Affine 2-Pieces
      } else if (strcmp(optarg,"gap-affine2p-dp")==0) {
        parameters.algorithm = alignment_gap_affine2p_dp;
      } else if (strcmp(optarg,"gap-affine2p-wfa")==0) {
        parameters.algorithm = alignment_gap_affine2p_wavefront;
#ifdef EXTERNAL_BENCHMARKS
      /*
       * External Algorithm
       */
      // External (BitPal)
      } else if (strcmp(optarg,"bitpal-edit")==0) {
        parameters.algorithm = alignment_bitpal_edit;
      } else if (strcmp(optarg,"bitpal-scored")==0) {
        parameters.algorithm = alignment_bitpal_scored;
      // External (Daligner)
      } else if (strcmp(optarg,"daligner")==0) {
        parameters.algorithm = alignment_daligner;
      // External (Diffutils)
      } else if (strcmp(optarg,"diffutils")==0) {
        parameters.algorithm = alignment_diffutils;
      // External (Edlib)
      } else if (strcmp(optarg,"edlib")==0) {
        parameters.algorithm = alignment_edlib;
      // External (Gaba)
      } else if (strcmp(optarg,"gaba-aband")==0) {
        parameters.algorithm = alignment_gaba_aband;
      // External (KSW2)
      } else if (strcmp(optarg,"ksw2-extz2-sse")==0) {
        parameters.algorithm = alignment_ksw2_extz2_sse;
      } else if (strcmp(optarg,"ksw2-extd2-sse")==0) {
        parameters.algorithm = alignment_ksw2_extd2_sse;
      // External (LV89)
      } else if (strcmp(optarg,"lv89")==0) {
        parameters.algorithm = alignment_lv89;
      // External (wfalm)
      } else if (strcmp(optarg,"wfalm")==0) {
        parameters.algorithm = alignment_wfalm;
      } else if (strcmp(optarg,"wfalm-lowmem")==0) {
        parameters.algorithm = alignment_wfalm_lowmem;
      } else if (strcmp(optarg,"wfalm-rec")==0) {
        parameters.algorithm = alignment_wfalm_rec;
#endif
      } else {
        fprintf(stderr,"Algorithm '%s' not recognized\n",optarg);
        exit(1);
      }
      break;
    }
    /*
     * Input/Output
     */
    case 'i':
      parameters.input_filename = optarg;
      break;
    case 'o':
      parameters.output_filename = optarg;
      break;
    case 800: // --output-full
      parameters.output_filename = optarg;
      parameters.output_full = true;
      break;
    case 'r': // --repeat
      parameters.target_bases_aligned = atoi(optarg);
      break;
    /*
     * Penalties
     */
    case 'p': { // --linear-penalties M,X,I
      char* sentinel = strtok(optarg,",");
      parameters.linear_penalties.match = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.linear_penalties.mismatch = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.linear_penalties.indel = atoi(sentinel);
      break;
    }
    case 'g': { // --affine-penalties M,X,O,E
      char* sentinel = strtok(optarg,",");
      parameters.affine_penalties.match = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.mismatch = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.gap_opening = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.gap_extension = atoi(sentinel);
      break;
    }
    case 900: { // --affine2p-penalties M,X,O1,E1,O2,E2
      char* sentinel = strtok(optarg,",");
      parameters.affine2p_penalties.match = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.mismatch = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.gap_opening1 = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.gap_extension1 = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.gap_opening2 = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine2p_penalties.gap_extension2 = atoi(sentinel);
      break;
    }
    /*
     * Wavefront parameters
     */
    case 1000: // --wfa-score-only
      parameters.wfa_score_only = true;
      break;
    case 1001: // --wfa-memory-mode in {'high','med','low','ultralow'}
      if (strcmp(optarg,"high")==0) {
        parameters.wfa_memory_mode = wavefront_memory_high;
      } else if (strcmp(optarg,"med")==0) {
        parameters.wfa_memory_mode = wavefront_memory_med;
      } else if (strcmp(optarg,"low")==0) {
        parameters.wfa_memory_mode = wavefront_memory_low;
      } else if (strcmp(optarg,"low")==0 || strcmp(optarg,"biwfa")==0) {
        parameters.wfa_memory_mode = wavefront_memory_ultralow;
      } else {
        fprintf(stderr,"Option '--wfa-memory-mode' must be in {'high','med','low','ultralow'}\n");
        exit(1);
      }
      break;
    case 1002: // --wfa-heuristic in {'none'|'banded-static'|'banded-adaptive'|'wfa-adaptive'|'xdrop'|'zdrop'}
      if (strcmp(optarg,"none")==0) {
        parameters.wfa_heuristic = wf_heuristic_none;
      } else if (strcmp(optarg,"banded-static")==0 || strcmp(optarg,"banded")==0) {
        parameters.wfa_heuristic = wf_heuristic_banded_static;
      } else if (strcmp(optarg,"banded-adaptive")==0) {
        parameters.wfa_heuristic = wf_heuristic_banded_adaptive;
      } else if (strcmp(optarg,"wfa-adaptive")==0) {
        parameters.wfa_heuristic = wf_heuristic_wfadaptive;
      } else if (strcmp(optarg,"xdrop")==0) {
        parameters.wfa_heuristic = wf_heuristic_xdrop;
      } else if (strcmp(optarg,"zdrop")==0) {
        parameters.wfa_heuristic = wf_heuristic_zdrop;
      } else {
        fprintf(stderr,"Option '--wf-heuristic' must be in {'none'|'banded-static'|'banded-adaptive'|'wfa-adaptive'|'xdrop'|'zdrop'}\n");
        exit(1);
      }
      break;
    case 1003: { // --wfa-heuristic-parameters  <P1>,<P2>[,<P3>]
      char* sentinel = strtok(optarg,",");
      const int p1 = atoi(sentinel);
      parameters.wfa_heuristic_p1 = p1;
      sentinel = strtok(NULL,",");
      const int p2 = atoi(sentinel);
      parameters.wfa_heuristic_p2 = p2;
      sentinel = strtok(NULL,",");
      if (sentinel != NULL) {
        const int p3 = atoi(sentinel);
        parameters.wfa_heuristic_p3 = p3;
      }
      break;
    }
    case 1005: // --wfa-max-memory
      parameters.wfa_max_memory = atol(optarg);
      break;
    case 1007: // --wfa-max-threads
      parameters.wfa_max_threads = atoi(optarg);
      break;
    /*
     * Other alignment parameters
     */
    case 2000: // --bandwidth
      parameters.bandwidth = atoi(optarg);
      break;
#ifdef EXTERNAL_BENCHMARKS
    case 2002: // --ksw2-approx-max-drop
      parameters.ksw2_approx_max__drop = true;
      break;
    case 2003: // --ksw2-bandwidth
      parameters.ksw2_bandwidth = atoi(optarg);
      break;
    case 2004: // --ksw2-zdrop
      parameters.ksw2_zdrop = atoi(optarg);
      break;
#endif
    /*
     * Misc
     */
    case 'c':
      if (optarg ==  NULL) { // default = score
        parameters.check_correct = true;
        parameters.check_score = true;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"display")==0) {
        parameters.check_display = true;
      } else if (strcasecmp(optarg,"correct")==0) {
        parameters.check_correct = true;
        parameters.check_score = false;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"score")==0) {
        parameters.check_correct = true;
        parameters.check_score = true;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"alignment")==0) {
        parameters.check_correct = true;
        parameters.check_score = true;
        parameters.check_alignments = true;
      } else {
        fprintf(stderr,"Option '--check' must be in {'correct','score','alignment'}\n");
        exit(1);
      }
      break;
    case 3001: // --check-distance in {'indel','edit','linear','affine','affine2p'}
      if (strcasecmp(optarg,"indel")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_INDEL;
      } else if (strcasecmp(optarg,"edit")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_EDIT;
      } else if (strcasecmp(optarg,"linear")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_LINEAR;
      } else if (strcasecmp(optarg,"affine")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE;
      } else if (strcasecmp(optarg,"affine2p")==0) {
        parameters.check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE2P;
      } else {
        fprintf(stderr,"Option '--check-distance' must be in {'indel','edit','linear','affine','affine2p'}\n");
        exit(1);
      }
      break;
    case 3002: // --check-bandwidth
      parameters.check_bandwidth = atoi(optarg);
      break;
    case 3003: // --plot
      parameters.plot = (optarg==NULL) ? 1000 : atoi(optarg);
      break;
    /*
     * System
     */
    case 't': // --num-threads
      parameters.num_threads = atoi(optarg);
      break;
    case 4000: // --batch-size
      parameters.batch_size = atoi(optarg);
      break;
    case 'P':
      parameters.progress = atoi(optarg);
      break;
    case 'v':
      if (optarg==NULL) {
        parameters.verbose = 1;
      } else {
        parameters.verbose = atoi(optarg);
        if (parameters.verbose < 0 || parameters.verbose > 3) {
          fprintf(stderr,"Option '--verbose' must be in {0,1,2,3}\n");
          exit(1);
        }
      }
      break;
    case 'h':
      usage();
      exit(1);
    // Other
    case '?': default:
      fprintf(stderr,"Option not recognized \n");
      exit(1);
    }
  }
  // Checks general
  if (parameters.algorithm!=alignment_test && parameters.input_filename==NULL) {
    fprintf(stderr,"Option --input is required \n");
    exit(1);
  }
  // Check 'bandwidth' parameter
  switch (parameters.algorithm) {
    case alignment_edit_dp_banded:
    case alignment_gap_affine_swg_banded:
      if (parameters.bandwidth == -1) {
        fprintf(stderr,"Parameter 'bandwidth' has to be provided for banded algorithms\n");
        exit(1);
      }
      break;
    default:
      if (parameters.bandwidth != -1) {
        fprintf(stderr,"Parameter 'bandwidth' has no effect with the selected algorithm\n");
        exit(1);
      }
      break;
  }
  // Check 'wfa-heuristic'
  switch (parameters.wfa_heuristic) {
    case wf_heuristic_banded_static:
    case wf_heuristic_xdrop:
    case wf_heuristic_zdrop:
      if (parameters.wfa_heuristic_p1 == -1 ||
          parameters.wfa_heuristic_p2 == -1) {
        fprintf(stderr,"Heuristic requires parameters '--wfa-heuristic-parameters' <P1>,<P2>\n");
        exit(1);
      }
      break;
    case wf_heuristic_banded_adaptive:
    case wf_heuristic_wfadaptive:
      if (parameters.wfa_heuristic_p1 == -1 ||
          parameters.wfa_heuristic_p2 == -1 ||
          parameters.wfa_heuristic_p3 == -1) {
        fprintf(stderr,"Heuristic requires parameters '--wfa-heuristic-parameters' <P1>,<P2>,<P3>\n");
        exit(1);
      }
      break;
    default:
      break;
  }
  // Checks parallel
  if (parameters.num_threads > 1) {
    if (parameters.plot > 0) {
      fprintf(stderr,"Parameter 'plot' disabled for parallel executions\n");
      parameters.plot = 0;
    }
  }
}
int main(int argc,char* argv[]) {
  // Parsing command-line options
  parse_arguments(argc,argv);
  // Select option
  if (parameters.algorithm == alignment_test) {
    align_pairwise_test();
  } else {
    // Execute benchmark
    align_benchmark_sequential();
  }
}
