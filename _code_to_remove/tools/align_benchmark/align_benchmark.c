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

/*
 * Generic parameters
 */
typedef struct {
  // I/O
  char *input_filename;
  char *output_filename;
  bool output_full;
  // I/O internals
  FILE* input_file;
  char* line1;
  char* line2;
  size_t line1_allocated;
  size_t line2_allocated;
  FILE* output_file;
  // Penalties
  affine_penalties_t affine_penalties;
  // Wavefront parameters
  bool wfa_score_only;
  bool wfa_bidirectional;
  int wfa_max_threads;
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
  int progress;
  int verbose;
} benchmark_args;
benchmark_args parameters = {
  // I/O
  .input_filename = NULL,
  .output_filename = NULL,
  .output_full = false,
  .output_file = NULL,
  // I/O internals
  .input_file = NULL,
  .line1 = NULL,
  .line2 = NULL,
  .line1_allocated = 0,
  .line2_allocated = 0,
  // Penalties
  .affine_penalties = {
      .match = 0,
      .mismatch = 4,
      .gap_opening = 6,
      .gap_extension = 2,
  },
  // Wavefront parameters
  .wfa_score_only = false,
  .wfa_bidirectional = true,
  .wfa_max_threads = 1,
  // Misc
  .check_bandwidth = -1,
  .check_display = false,
  .check_correct = false,
  .check_score = false,
  .check_alignments = false,
  .check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE,
  .plot = 0,
  // System
  .progress = 100000,
  .verbose = 0,
};

/*
 * Configuration
 */
wavefront_aligner_t* align_input_configure_wavefront(
    align_input_t* const align_input) {
  // Set attributes
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  if (parameters.wfa_score_only) {
    attributes.alignment_scope = compute_score;
  }
  // WF-Heuristic
  attributes.heuristic.strategy = wf_heuristic_none;
  // Select flavor
  attributes.affine_penalties = parameters.affine_penalties;
  // Misc
  attributes.bidirectional_alignment = true; //parameters.wfa_bidirectional;
  attributes.plot_params.plot_enabled = (parameters.plot > 0);
  attributes.plot_params.resolution_points = parameters.plot;
  attributes.system.verbose = parameters.verbose;
  attributes.system.max_num_threads = parameters.wfa_max_threads;
  // Allocate
  return wavefront_aligner_new(&attributes);
}
void align_input_configure_global(
    align_input_t* const align_input) {
  // Clear
  benchmark_align_input_clear(align_input);
  // Penalties
  align_input->affine_penalties = parameters.affine_penalties;
  // Output
  align_input->output_file = parameters.output_file;
  align_input->output_full = parameters.output_full;
  // MM
  align_input->mm_allocator = mm_allocator_new(BUFFER_SIZE_1M);
  // WFA
  align_input->wf_aligner = align_input_configure_wavefront(align_input);
  // PROFILE/STATS
  timer_reset(&align_input->timer);
  // DEBUG
  align_input->debug_flags = 0;
  align_input->debug_flags |= parameters.check_metric;
  if (parameters.check_display) align_input->debug_flags |= ALIGN_DEBUG_DISPLAY_INFO;
  if (parameters.check_correct) align_input->debug_flags |= ALIGN_DEBUG_CHECK_CORRECT;
  if (parameters.check_score) align_input->debug_flags |= ALIGN_DEBUG_CHECK_SCORE;
  if (parameters.check_alignments) align_input->debug_flags |= ALIGN_DEBUG_CHECK_ALIGNMENT;
  align_input->check_affine_penalties = &parameters.affine_penalties;
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
void align_benchmark_plot_wf(
    align_input_t* const align_input,
    const int seq_id) {
  // Setup filename
  char filename[500];
  if (parameters.output_filename != NULL) {
    sprintf(filename,"%s.%03d.wfa",parameters.output_filename,seq_id);
  } else {
    sprintf(filename,"%s.%03d.wfa",parameters.input_filename,seq_id);
  }
  // Open file
  FILE* const wf_plot = fopen(filename,"w");
  wavefront_plot_print(wf_plot,align_input->wf_aligner);
  fclose(wf_plot);
}
/*
 * Benchmark
 */
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
  int seqs_processed = 0, progress = 0;
  while (true) {
    // Read input sequence-pair
    const bool input_read = align_benchmark_read_input(
        parameters.input_file,&parameters.line1,&parameters.line2,
        &parameters.line1_allocated,&parameters.line2_allocated,
        seqs_processed,&align_input);
    if (!input_read) break;
    // Execute the algorithm
    benchmark_gap_affine_wavefront(&align_input,&parameters.affine_penalties);
    // Update progress
    ++seqs_processed;
    if (++progress == parameters.progress) {
      progress = 0;
      align_benchmark_print_progress(seqs_processed);
    }
    // Plot
    if (parameters.plot > 0) align_benchmark_plot_wf(&align_input,seqs_processed);
  }
  // Print benchmark results
  timer_stop(&parameters.timer_global);
  align_benchmark_print_results(&align_input,seqs_processed,true);
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
      "USE: ./align_benchmark -i <input>                                       \n"
      "        [Input & Output]                                                \n"
      "          --input|i <File>                                              \n"
      "          --output|o <File>                                             \n"
      "          --output-full <File>                                          \n"
      "        [Penalties]                                                     \n"
      "          --affine-penalties|g M,X,O,E                                  \n"
      "        [Wavefront parameters]                                          \n"
      "          --wfa-bidirectional                                           \n"
      "          --wfa-max-threads <INT> (intra-parallelism; default=1)        \n"
      "        [Misc]                                                          \n"
      "          --check|c 'correct'|'score'|'alignment'                       \n"
      "          --check-distance 'indel'|'edit'|'linear'|'affine'|'affine2p'  \n"
      "          --check-bandwidth <INT>                                       \n"
      "        [System]                                                        \n"
      "          --help|h                                                      \n");
}
void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    /* Input */
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "output-full", required_argument, 0, 800 },
    /* Penalties */
    { "affine-penalties", required_argument, 0, 'g' },
    /* Wavefront parameters */
    { "wfa-bidirectional", no_argument, 0, 1006 },
    { "wfa-max-threads", required_argument, 0, 1007 },
    /* Misc */
    { "check", required_argument, 0, 'c' },
    { "check-distance", required_argument, 0, 3001 },
    { "check-bandwidth", required_argument, 0, 3002 },
    { "plot", optional_argument, 0, 3003 },
    /* System */
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
    c=getopt_long(argc,argv,"i:o:g:c:P:v::h",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
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
    /*
     * Penalties
     */
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
    /*
     * Wavefront parameters
     */
    case 1000: // --wfa-score-only
      parameters.wfa_score_only = true;
      break;
    case 1006: // --wfa-bidirectional
      parameters.wfa_bidirectional = true;
      break;
    case 1007: // --wfa-max-threads
      parameters.wfa_max_threads = atoi(optarg);
      break;
    /*
     * Misc
     */
    case 'c':
      if (strcasecmp(optarg,"display")==0) {
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
  if (parameters.input_filename==NULL) {
    fprintf(stderr,"Option --input is required \n");
    exit(1);
  }
}
int main(int argc,char* argv[]) {
  // Parsing command-line options
  parse_arguments(argc,argv);
  // Run
  align_benchmark_sequential();
}
