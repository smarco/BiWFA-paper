# wfalm
Refinements of the WFA alignment algorithm with better complexity

## Introduction

This repository contains implementations of the [WFA algorithm](https://academic.oup.com/bioinformatics/article/37/4/456/5904262) as well as a few variations on it.

- A low-memory variation that reduces space complexity from O(s^2) to O(s^3/2).
- A recursive variation that reduces space complexity from O(s^2) to O(s log s) at the cost of an O(log s)-factor run time penalty.
- A memory-adaptive variation that dynamically chooses between the standard, low-memory, and recursive variations as necessary to keep memory below a user-defined limit.
- A suffix tree variation that reduces the time complexity from O(sN) to O(s^2 + N). 

The three memory-reducing variations are eminently practical, but the suffix tree algorithm has large constants that make it take more time in the majority of common use cases.

## Installation

Because many users will not need the suffix tree algorithm, the library is implemented in two headers. The non-suffix tree algorithm is implemented in `wfa_lm.hpp`. For environments with modern x86 processors, this header is entirely self-contained and depends on nothing except the C++14 standard library. Simply copy `wfa_lm.hpp` into the include directory of your project. You may need to add `-msse4.1` to your compiler invocation.

The library is usable in environments with non-x86 processors (such as many Apple notebooks), but it depends on [SIMDE](https://github.com/simd-everywhere/simde) for portability. SIMDE is a header-only library as well. To use it, first check out the SIMDE submodule:

	git submodule update --init --recursive

Next, copy `wfa_lm.hpp` and `simde/simde/` into the include directory of your project.

The suffix tree algorithm is implemented in `wfa_lm_st.hpp`. It is also single-header, but it depends on [SDSL v3](https://github.com/xxsds/sdsl-lite) for its implementation of a suffix tree. To use it, check out the submodules as described above and copy the headers in `external/sdsl-lite/include` to your project's include directory.

## Citation

Eizenga, JM, and Paten, B (2022) Improving the time and space complexity of the WFA algorithm and generalizing its scoring. _bioRxiv_. DOI [10.1101/2022.01.12.476087](https://doi.org/10.1101/2022.01.12.476087)
