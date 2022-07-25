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
 */

#include "wavefront_bialign.h"

#include "wavefront_align.h"
#include "wavefront_extend.h"
#include "wavefront_compute.h"
#include "wavefront_compute_edit.h"
#include "wavefront_compute_linear.h"
#include "wavefront_compute_affine.h"
#include "wavefront_compute_affine2p.h"
#include "wavefront_backtrace.h"

/*
 * Setup
 */
void wavefront_bialign_junction_configure(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const distance_metric_t distance_metric,
    alignment_form_t* const form,
    const affine_matrix_type component_begin,
    const affine_matrix_type component_end) {
  // Resize wavefront aligner
  wavefront_aligner_resize(wf_aligner_forward,pattern,pattern_length,text,text_length,false);
  wavefront_aligner_resize(wf_aligner_reverse,pattern,pattern_length,text,text_length,true);
  // Configure form forward and reverse
  alignment_span_t span_forward =
      (form->pattern_begin_free > 0 || form->text_begin_free > 0) ? alignment_endsfree : alignment_end2end;
  alignment_form_t form_forward = {
      .span = span_forward,
      .pattern_begin_free = form->pattern_begin_free,
      .pattern_end_free = 0,
      .text_begin_free = form->text_begin_free,
      .text_end_free = 0,
  };
  alignment_span_t span_reverse =
      (form->pattern_end_free > 0 || form->text_end_free > 0) ? alignment_endsfree : alignment_end2end;
  alignment_form_t form_reverse = {
      .span = span_reverse,
      .pattern_begin_free = form->pattern_end_free,
      .pattern_end_free = 0,
      .text_begin_free = form->text_end_free,
      .text_end_free = 0,
  };
  // Configure WF-compute function (global)
  switch (distance_metric) {
    case indel:
    case edit:
      wf_aligner_forward->align_status.wf_align_compute = &wavefront_compute_edit;
      break;
    case gap_linear:
      wf_aligner_forward->align_status.wf_align_compute = &wavefront_compute_linear;
      break;
    case gap_affine:
      wf_aligner_forward->align_status.wf_align_compute = &wavefront_compute_affine;
      break;
    case gap_affine_2p:
      wf_aligner_forward->align_status.wf_align_compute = &wavefront_compute_affine2p;
      break;
    default:
      fprintf(stderr,"[WFA] Distance function not implemented\n");
      exit(1);
      break;
  }
  // Initialize wavefront (forward)
  wf_aligner_forward->alignment_form = form_forward;
  wf_aligner_forward->component_begin = component_begin;
  if (span_forward == alignment_end2end) {
    wavefront_align_end2end_initialize(wf_aligner_forward);
  } else {
    wavefront_align_endsfree_initialize(wf_aligner_forward,pattern_length,text_length);
  }
  // Initialize wavefront (reverse)
  wf_aligner_reverse->alignment_form = form_reverse;
  wf_aligner_reverse->component_begin = component_end;
  if (span_reverse == alignment_end2end) {
    wavefront_align_end2end_initialize(wf_aligner_reverse);
  } else {
    wavefront_align_endsfree_initialize(wf_aligner_reverse,pattern_length,text_length);
  }
}
/*
 * Wavefronts junction check
 */
void wavefront_bialign_junction_i2i(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    wavefront_t* const iwf_forward,
    wavefront_t* const iwf_reverse,
    const int k_forward,
    const int k_reverse,
    wf_bialign_junction_t* const junction) {
  // Parameters
  const int text_length = wf_aligner_forward->text_length;
  const int gap_opening = wf_aligner_forward->penalties.gap_opening1;
  // Check junction i2i
  const wf_offset_t ioffset_forward = iwf_forward->offsets[k_forward];
  const wf_offset_t ioffset_reverse = iwf_reverse->offsets[k_reverse];
  const int ih_forward = WAVEFRONT_H(k_forward,ioffset_forward);
  const int ih_reverse = WAVEFRONT_H(k_reverse,ioffset_reverse);
  if (ih_forward + ih_reverse >= text_length &&
      score_forward + score_reverse - gap_opening < junction->score) {
    junction->score_forward = score_forward;
    junction->score_reverse = score_reverse;
    junction->score = score_forward + score_reverse - gap_opening;
    junction->k_forward = k_forward;
    junction->k_reverse = k_reverse;
    junction->offset_forward = ih_forward;
    junction->offset_reverse = ih_reverse;
    junction->component = affine_matrix_I;
  }
}
void wavefront_bialign_junction_d2d(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    wavefront_t* const dwf_forward,
    wavefront_t* const dwf_reverse,
    const int k_forward,
    const int k_reverse,
    wf_bialign_junction_t* const junction) {
  // Parameters
  const int text_length = wf_aligner_forward->text_length;
  const int gap_opening = wf_aligner_forward->penalties.gap_opening1;
  // Check junction d2d
  const wf_offset_t doffset_forward = dwf_forward->offsets[k_forward];
  const wf_offset_t doffset_reverse = dwf_reverse->offsets[k_reverse];
  const int dh_forward = WAVEFRONT_H(k_forward,doffset_forward);
  const int dh_reverse = WAVEFRONT_H(k_reverse,doffset_reverse);
  if (dh_forward + dh_reverse >= text_length &&
      score_forward + score_reverse - gap_opening < junction->score) {
    junction->score_forward = score_forward;
    junction->score_reverse = score_reverse;
    junction->score = score_forward + score_reverse - gap_opening;
    junction->k_forward = k_forward;
    junction->k_reverse = k_reverse;
    junction->offset_forward = dh_forward;
    junction->offset_reverse = dh_reverse;
    junction->component = affine_matrix_D;
  }
}
void wavefront_bialign_junction_m2m(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    wavefront_t* const mwf_forward,
    wavefront_t* const mwf_reverse,
    const int k_forward,
    const int k_reverse,
    wf_bialign_junction_t* const junction) {
  // Parameters
  const int text_length = wf_aligner_forward->text_length;
  // Check junction m2m
  const wf_offset_t moffset_forward = mwf_forward->offsets[k_forward];
  const wf_offset_t moffset_reverse = mwf_reverse->offsets[k_reverse];
  const int mh_forward = WAVEFRONT_H(k_forward,moffset_forward);
  const int mh_reverse = WAVEFRONT_H(k_reverse,moffset_reverse);
  if (mh_forward + mh_reverse >= text_length &&
      score_forward + score_reverse < junction->score) {
    junction->score_forward = score_forward;
    junction->score_reverse = score_reverse;
    junction->score = score_forward + score_reverse;
    junction->k_forward = k_forward;
    junction->k_reverse = k_reverse;
    junction->offset_forward = moffset_forward;
    junction->offset_reverse = moffset_reverse;
    junction->component = affine_matrix_M;
  }
}
/*
 * Wavefront diagonal overlap
 */
void wavefront_bialign_forward_overlap_diagonal(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    const int k_forward,
    wf_bialign_junction_t* const best_junction) {
  // Parameters
  const int pattern_length = wf_aligner_forward->pattern_length;
  const int text_length = wf_aligner_forward->text_length;
  const bool memory_modular = wf_aligner_forward->wf_components.memory_modular;
  const int max_score_scope = wf_aligner_forward->wf_components.max_score_scope;
  // Fetch wavefronts forward
  const int score_forward_mod = (memory_modular) ? (score_forward % max_score_scope) : score_forward;
  wavefront_t* const mwf_forward = wf_aligner_forward->wf_components.mwavefronts[score_forward_mod];
  wavefront_t* const dwf_forward = wf_aligner_forward->wf_components.d1wavefronts[score_forward_mod];
  wavefront_t* const iwf_forward = wf_aligner_forward->wf_components.i1wavefronts[score_forward_mod];
  // Traverse all reverse scores
  const int k_reverse = WAVEFRONT_K_INVERSE(k_forward,pattern_length,text_length);
  int i;
  for (i=0;i<max_score_scope;++i) {
    // Compute score
    const int score_rev = score_reverse - i;
    if (score_rev < 0) break;
    const int score_reverse_mod = (memory_modular) ? (score_rev % max_score_scope) : score_rev;
    // Fetch wavefronts reverse
    wavefront_t* const mwf_reverse = wf_aligner_reverse->wf_components.mwavefronts[score_reverse_mod];
    wavefront_t* const dwf_reverse = wf_aligner_reverse->wf_components.d1wavefronts[score_reverse_mod];
    wavefront_t* const iwf_reverse = wf_aligner_reverse->wf_components.i1wavefronts[score_reverse_mod];
    // Check junction m2m
    if (mwf_reverse != NULL && mwf_reverse->lo <= k_reverse && k_reverse <= mwf_reverse->hi) {
      wavefront_bialign_junction_m2m(
          wf_aligner_forward,wf_aligner_reverse,
          score_forward,score_rev,
          mwf_forward,mwf_reverse,
          k_forward,k_reverse,best_junction);
    }
    // Check junction d2d
    if (dwf_forward != NULL && dwf_reverse != NULL && dwf_reverse->lo <= k_reverse && k_reverse <= dwf_reverse->hi) {
      wavefront_bialign_junction_d2d(
          wf_aligner_forward,wf_aligner_reverse,
          score_forward,score_rev,
          dwf_forward,dwf_reverse,
          k_forward,k_reverse,best_junction);
    }
    // Check junction i2i
    if (iwf_forward != NULL && iwf_reverse != NULL && iwf_reverse->lo <= k_reverse && k_reverse <= iwf_reverse->hi) {
      wavefront_bialign_junction_i2i(
          wf_aligner_forward,wf_aligner_reverse,
          score_forward,score_rev,
          iwf_forward,iwf_reverse,
          k_forward,k_reverse,best_junction);
    }
  }
}
void wavefront_bialign_reverse_overlap_diagonal(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    const int k_reverse,
    wf_bialign_junction_t* const best_junction) {
  // Parameters
  const int pattern_length = wf_aligner_forward->pattern_length;
  const int text_length = wf_aligner_forward->text_length;
  const bool memory_modular = wf_aligner_forward->wf_components.memory_modular;
  const int max_score_scope = wf_aligner_forward->wf_components.max_score_scope;
  // Fetch wavefronts reverse
  const int score_reverse_mod = (memory_modular) ? (score_reverse % max_score_scope) : score_reverse;
  wavefront_t* const mwf_reverse = wf_aligner_reverse->wf_components.mwavefronts[score_reverse_mod];
  wavefront_t* const dwf_reverse = wf_aligner_reverse->wf_components.d1wavefronts[score_reverse_mod];
  wavefront_t* const iwf_reverse = wf_aligner_reverse->wf_components.i1wavefronts[score_reverse_mod];
  // Traverse all reverse scores
  const int k_forward = WAVEFRONT_K_INVERSE(k_reverse,pattern_length,text_length);
  int i;
  for (i=0;i<max_score_scope;++i) {
    // Compute score
    const int score_for = score_forward - i;
    if (score_for < 0) break;
    const int score_for_mod = (memory_modular) ? (score_for % max_score_scope) : score_for;
    // Fetch wavefronts forward
    wavefront_t* const mwf_forward = wf_aligner_forward->wf_components.mwavefronts[score_for_mod];
    wavefront_t* const dwf_forward = wf_aligner_forward->wf_components.d1wavefronts[score_for_mod];
    wavefront_t* const iwf_forward = wf_aligner_forward->wf_components.i1wavefronts[score_for_mod];
    // Check junction m2m
    if (mwf_forward != NULL && mwf_forward->lo <= k_forward && k_forward <= mwf_forward->hi) {
      wavefront_bialign_junction_m2m(
          wf_aligner_forward,wf_aligner_reverse,
          score_for,score_reverse,
          mwf_forward,mwf_reverse,
          k_forward,k_reverse,best_junction);
    }
    // Check junction d2d
    if (dwf_reverse != NULL && dwf_forward != NULL && dwf_forward->lo <= k_forward && k_forward <= dwf_forward->hi) {
      wavefront_bialign_junction_d2d(
          wf_aligner_forward,wf_aligner_reverse,
          score_for,score_reverse,
          dwf_forward,dwf_reverse,
          k_forward,k_reverse,best_junction);
    }
    // Check junction i2i
    if (iwf_reverse != NULL && iwf_forward != NULL && iwf_forward->lo <= k_forward && k_forward <= iwf_forward->hi) {
      wavefront_bialign_junction_i2i(
          wf_aligner_forward,wf_aligner_reverse,
          score_for,score_reverse,
          iwf_forward,iwf_reverse,
          k_forward,k_reverse,best_junction);
    }
  }
}
/*
 * Wavefront Bidirectional (find overlaps)
 */
void wavefront_bialign_forward_overlap(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    wf_bialign_junction_t* const junction) {
  // Parameters
  const bool memory_modular = wf_aligner_forward->wf_components.memory_modular;
  const int max_score_scope = wf_aligner_forward->wf_components.max_score_scope;
  // DEBUG
  //wavefront_aligner_print(stderr,wf_aligner_forward,score_forward,score_forward,6,0); // DEBUG
  //wavefront_aligner_print(stderr,wf_aligner_reverse,score_reverse,score_reverse,6,0); // DEBUG
  // Fetch m-wavefronts forward
  const int score_forward_mod = (memory_modular) ? (score_forward % max_score_scope) : score_forward;
  wavefront_t* const mwf_forward = wf_aligner_forward->wf_components.mwavefronts[score_forward_mod];
  if (mwf_forward == NULL) return;
  // Traverse all diagonals and look for overlaps
  int k_forward;
  for (k_forward=mwf_forward->lo;k_forward<=mwf_forward->hi;k_forward++) {
    // Find overlaps on diagonal
    wf_bialign_junction_t diagonal_junction = { .score = INT_MAX };
    wavefront_bialign_forward_overlap_diagonal(
        wf_aligner_forward,wf_aligner_reverse,
        score_forward,score_reverse,k_forward,
        &diagonal_junction);
    if (diagonal_junction.score < junction->score) {
      *junction = diagonal_junction;
    }
  }
}
void wavefront_bialign_reverse_overlap(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const int score_forward,
    const int score_reverse,
    wf_bialign_junction_t* const junction) {
  // Parameters
  const bool memory_modular = wf_aligner_forward->wf_components.memory_modular;
  const int max_score_scope = wf_aligner_forward->wf_components.max_score_scope;
  // DEBUG
  //wavefront_aligner_print(stderr,wf_aligner_forward,score_forward,score_forward,6,0); // DEBUG
  //wavefront_aligner_print(stderr,wf_aligner_reverse,score_reverse,score_reverse,6,0); // DEBUG
  // Fetch m-wavefronts forward
  const int score_reverse_mod = (memory_modular) ? (score_reverse % max_score_scope) : score_reverse;
  wavefront_t* const mwf_reverse = wf_aligner_reverse->wf_components.mwavefronts[score_reverse_mod];
  if (mwf_reverse == NULL) return;
  // Traverse all diagonals and look for overlaps
  int k_reverse;
  for (k_reverse=mwf_reverse->lo;k_reverse<=mwf_reverse->hi;k_reverse++) {
    // Find overlaps on diagonal
    wf_bialign_junction_t diagonal_junction = { .score = INT_MAX };
    wavefront_bialign_reverse_overlap_diagonal(
        wf_aligner_forward,wf_aligner_reverse,
        score_forward,score_reverse,
        k_reverse,&diagonal_junction);
    if (diagonal_junction.score < junction->score) {
      *junction = diagonal_junction;
    }
  }
}
/*
 * Wavefront Bidirectional (find junction)
 */
void wavefront_bialign_junction_find(
    wavefront_aligner_t* const wf_aligner_forward,
    wavefront_aligner_t* const wf_aligner_reverse,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const distance_metric_t distance_metric,
    alignment_form_t* const form,
    const affine_matrix_type component_begin,
    const affine_matrix_type component_end,
    wf_bialign_junction_t* const junction) {
  // Configure bialign
  wavefront_bialign_junction_configure(
      wf_aligner_forward,wf_aligner_reverse,
      pattern,pattern_length,text,text_length,
      distance_metric,form,component_begin,component_end);
  // Compute wavefronts of increasing score until both wavefronts overlap
  const int max_antidiagonal = DPMATRIX_ANTIDIAGONAL(pattern_length,text_length) - 1;
  void (*wf_align_compute)(wavefront_aligner_t* const,const int) = wf_aligner_forward->align_status.wf_align_compute;
  junction->score = INT_MAX;
  int score_forward = 0, score_reverse = 0;
  int forward_max_ak = wavefront_extend_end2end_max(wf_aligner_forward,score_forward);
  int reverse_max_ak = wavefront_extend_end2end_max(wf_aligner_reverse,score_reverse);
  int max_ak;
  bool last_wf_forward;
  while (true) {
    // Check if they are close to collision
    if (forward_max_ak + reverse_max_ak >= max_antidiagonal) break;
    // Compute-next & extend wavefront forward
    ++score_forward;
    (*wf_align_compute)(wf_aligner_forward,score_forward);
    max_ak = wavefront_extend_end2end_max(wf_aligner_forward,score_forward);
    if (forward_max_ak < max_ak) forward_max_ak = max_ak;
    last_wf_forward = true;
    // Check if they are close to collision
    if (forward_max_ak + reverse_max_ak >= max_antidiagonal) break;
    // Compute-next & extend wavefront reverse
    ++score_reverse;
    (*wf_align_compute)(wf_aligner_reverse,score_reverse);
    max_ak = wavefront_extend_end2end_max(wf_aligner_reverse,score_reverse);
    if (reverse_max_ak < max_ak) reverse_max_ak = max_ak;
    last_wf_forward = false;
  }
  // Advance until overlap is found
  const int max_score_scope = wf_aligner_forward->wf_components.max_score_scope;
  const int gap_opening = wf_aligner_forward->penalties.gap_opening1;
  int max_edge_cost = wf_aligner_forward->penalties.mismatch;
  switch (distance_metric) {
    case gap_affine_2p:
      max_edge_cost = MAX(max_edge_cost, wf_aligner->penalties.gap_opening1 + wf_aligner->penalties.gap_extension1);
      max_edge_cost = MAX(max_edge_cost, wf_aligner->penalties.gap_opening2 + wf_aligner->penalties.gap_extension2);
      break;
    case gap_affine:
      max_edge_cost = MAX(max_edge_cost, wf_aligner->penalties.gap_opening1 + wf_aligner->penalties.gap_extension1);
      break;
    case gap_linear:
      max_edge_cost = MAX(max_edge_cost, wf_aligner->penalties.gap_opening1);
  }
  while (true) {
    if (last_wf_forward) {
      // Check overlapping wavefronts
      const int min_score_reverse = (score_reverse > max_score_scope-1) ? score_reverse - (max_score_scope-1) : 0;
      if (score_forward + min_score_reverse - gap_opening - max_edge_cost >= junction->score) break; // Done!
      wavefront_bialign_forward_overlap(wf_aligner_forward,wf_aligner_reverse,score_forward,score_reverse,junction);
      // Compute-next and extend reverse-wavefront
      ++score_reverse;
      (*wf_align_compute)(wf_aligner_reverse,score_reverse);
      wavefront_extend_end2end(wf_aligner_reverse,score_reverse);
    }
    // Check overlapping wavefronts
    const int min_score_forward = (score_forward > max_score_scope-1) ? score_forward - (max_score_scope-1) : 0;
    if (min_score_forward + score_reverse - gap_opening - max_edge_cost >= junction->score) break; // Done!
    wavefront_bialign_reverse_overlap(wf_aligner_forward,wf_aligner_reverse,score_forward,score_reverse,junction);
    // Compute-next and extend forward-wavefront
    ++score_forward;
    (*wf_align_compute)(wf_aligner_forward,score_forward);
    wavefront_extend_end2end(wf_aligner_forward,score_forward);
    // Enable always
    last_wf_forward = true;
  }
}
void wavefront_bialign(
    wavefront_aligner_t* const wf_aligner,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    alignment_form_t* const form,
    const affine_matrix_type component_begin,
    const affine_matrix_type component_end,
    const int score_remaining,
    cigar_t* const cigar,
    const int rec_level) {
  // Check trivial cases
  if (text_length == 0) {
    cigar_append_deletion(cigar,pattern_length);
    return;
  }
  if (pattern_length == 0) {
    cigar_append_insertion(cigar,text_length);
    return;
  }
  // Check score of the remaining alignment
  if (score_remaining <= 100) {
    // Align the remaining
    wf_aligner->bidirectional_alignment = false; // Hacky...
    wf_aligner->component_begin = component_begin;
    wf_aligner->component_end = component_end;
    wf_aligner->alignment_form = *form;
    wavefront_align(wf_aligner,
        pattern,pattern_length,
        text,text_length);
    wf_aligner->bidirectional_alignment = true; // Hacky...
    cigar_append(cigar,&wf_aligner->cigar);
    return;
  }
  // Find junction in the alignment
  wf_bialign_junction_t junction;
  wavefront_bialign_junction_find(
      wf_aligner->aligner_forward,wf_aligner->aligner_reverse,
      pattern,pattern_length,text,text_length,
      wf_aligner->penalties.distance_metric,
      form,component_begin,component_end,&junction);
  // Junction found. Align the parts
  const int junction_h = WAVEFRONT_H(junction.k_forward,junction.offset_forward);
  const int junction_v = WAVEFRONT_V(junction.k_forward,junction.offset_forward);
  // DEBUG
  if (wf_aligner->system.verbose == 1) {
    //fprintf(stderr,"[WFA::BiAlign] ");
    //int i; for (i=0;i<rec_level;++i) fprintf(stderr,"   ");
    //fprintf(stderr,"[%d] Junction at (h,v,score) = (%d,%d,%d)\n",
    //    rec_level,junction_h,junction_v,junction.score);
  }
  // Align half_0
  alignment_span_t span_0 =
      (form->pattern_begin_free > 0 || form->text_begin_free > 0) ? alignment_endsfree : alignment_end2end;
  alignment_form_t form_0 = {
      .span = span_0,
      .pattern_begin_free = form->pattern_begin_free,
      .pattern_end_free = 0,
      .text_begin_free = form->text_begin_free,
      .text_end_free = 0,
  };
  wavefront_bialign(
      wf_aligner,pattern,junction_v,text,junction_h,
      &form_0,component_begin,junction.component,
      junction.score_forward,cigar,rec_level+1);
  // Align half_1
  alignment_span_t span_1 =
      (form->pattern_end_free > 0 || form->text_end_free > 0) ? alignment_endsfree : alignment_end2end;
  alignment_form_t form_1 = {
      .span = span_1,
      .pattern_begin_free = 0,
      .pattern_end_free = form->pattern_end_free,
      .text_begin_free = 0,
      .text_end_free = form->text_end_free,
  };
  wavefront_bialign(wf_aligner,
      pattern+junction_v,pattern_length-junction_v,
      text+junction_h,text_length-junction_h,
      &form_1,junction.component,component_end,
      junction.score_reverse,cigar,rec_level+1);
  // Set score
  cigar->score = -junction.score;
}





