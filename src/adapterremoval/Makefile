

PKG_CXXFLAGS = -std=c++11 -O3
PKG_CPPFLAGS = -DR_ENV_FLAG -pthread -lz -DAR_GZIP_SUPPORT

include $(R_HOME)/etc${R_ARCH}/Makeconf

SRCS = adapterset.cpp \
alignment.cpp \
argparse.cpp \
debug.cpp \
demultiplex.cpp \
fastq.cpp \
fastq_enc.cpp \
fastq_io.cpp \
linereader.cpp \
linereader_joined.cpp \
main.cpp \
main_adapter_id.cpp \
main_adapter_rm.cpp \
main_demultiplex.cpp \
scheduler.cpp \
strutils.cpp \
threads.cpp \
timer.cpp \
trimmed_reads.cpp \
userconfig.cpp \

ADRM_LIBRARY = adrm


.PHONY: all clean

all: lib$(ADRM_LIBRARY).a

OBJS = $(SRCS:.cpp=.o)

lib$(ADRM_LIBRARY).a:$(OBJS)
	$(AR) cq $@ $^

clean:
	rm -f *.a
	rm -f *.o
