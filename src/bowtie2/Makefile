#
# Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
#
# This file is part of Bowtie 2.
#
# Bowtie 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Bowtie 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
#

#
# Makefile for bowtie, bowtie2-build, bowtie2-inspect
#


INC_EXT = -I third_party
LIBS_BT2 = -pthread
#LIBS_BT2 = -lz -lpthread
#CXX_STD = CXX11

MM_DEF = -DBOWTIE_MM


EXTRA_FLAGS = -DPOPCNT_CAPABILITY

SHARED_CPPS = ccnt_lut.cpp ref_read.cpp alphabet.cpp shmem.cpp \
              edit.cpp bt2_idx.cpp bt2_io.cpp bt2_util.cpp \
              reference.cpp ds.cpp multikey_qsort.cpp limit.cpp \
			  random_source.cpp tinythread.cpp



SEARCH_CPPS = qual.cpp pat.cpp sam.cpp \
              read_qseq.cpp aligner_seed_policy.cpp \
              aligner_seed.cpp \
			  aligner_seed2.cpp \
			  aligner_sw.cpp \
			  aligner_sw_driver.cpp aligner_cache.cpp \
			  aligner_result.cpp ref_coord.cpp mask.cpp \
			  pe.cpp aln_sink.cpp dp_framer.cpp \
			  scoring.cpp presets.cpp unique.cpp \
			  simple_func.cpp \
			  random_util.cpp \
			  aligner_bt.cpp sse_util.cpp \
			  aligner_swsse.cpp outq.cpp \
			  aligner_swsse_loc_i16.cpp \
			  aligner_swsse_ee_i16.cpp \
			  aligner_swsse_loc_u8.cpp \
			  aligner_swsse_ee_u8.cpp \
			  aligner_driver.cpp

SEARCH_CPPS_MAIN = $(SEARCH_CPPS) bowtie_main.cpp

DP_CPPS = qual.cpp aligner_sw.cpp aligner_result.cpp ref_coord.cpp mask.cpp \
          simple_func.cpp sse_util.cpp aligner_bt.cpp aligner_swsse.cpp \
		  aligner_swsse_loc_i16.cpp aligner_swsse_ee_i16.cpp \
		  aligner_swsse_loc_u8.cpp aligner_swsse_ee_u8.cpp scoring.cpp

BUILD_CPPS = diff_sample.cpp
BUILD_CPPS_MAIN = $(BUILD_CPPS) bowtie_build_main.cpp


VERSION = 2.3.2


SSE_FLAG=-msse2


RELEASE_FLAGS  = -O3 -m64 $(SSE_FLAG) -funroll-loops -g3
RELEASE_DEFS   = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS) $(EXTRA_FLAGS)\""
NOASSERT_FLAGS = -DNDEBUG
FILE_FLAGS     = -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE




DEFS=-fno-strict-aliasing \
     -DBOWTIE2_VERSION="\"2.3.2\"" \
     -DBUILD_HOST="\"Rbowtie2\"" \
     -DBUILD_TIME="\"2017\"" \
     -DCOMPILER_VERSION="\"C++11\"" \
     $(FILE_FLAGS) \
     $(MM_DEF) \


CPPFLAGS = $(RELEASE_DEFS) $(EXTRA_FLAGS) $(DEFS) -DBOWTIE2 \
               $(NOASSERT_FLAGS) $(INC_EXT) $(LIBS_BT2)  \
               $(RELEASE_FLAGS) -Wall


CXXFLAGS = $(RELEASE_FLAGS) -Wall



SRCS_SEARCH = bt2_search.cpp $(SEARCH_CPPS) $(SHARED_CPPS) bowtie_main.cpp

SRCS_BUILD = bt2_build.cpp $(SHARED_CPPS) $(HEADERS) bowtie_build_main.cpp diff_sample.cpp

OBJS_SEARCH = $(SRCS_SEARCH:.cpp=.o)

OBJS_BUILD = $(SRCS_BUILD:.cpp=.o)

.PHONY: all clean bowtie2-align-s bowtie2-build-s install

all: bowtie2-build-s bowtie2-align-s install



.cc.o:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< 

 

bowtie2-build-s: $(OBJS_BUILD)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^
	
bowtie2-align-s: $(OBJS_SEARCH)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ 

install:
	cp bowtie2-align-s ../../inst
	cp bowtie2-build-s ../../inst

clean:
	rm -f ../../inst/bowtie2-align-s
	rm -f ../../inst/bowtie2-build-s
	rm -f bowtie2-align-s
	rm -f bowtie2-build-s
	rm -f *.o
