bin_PROGRAMS = petagstats

samtools_INCLUDE = $(top_srcdir)/lib/samtools-0.1.20
samtools_LIB = $(top_srcdir)/lib/samtools-0.1.20

CXXFLAGS = -g -O3 -std=c++11
petagstats_CPPFLAGS = -Wall -w -I$(samtools_INCLUDE)
petagstats_LDFLAGS = -L$(samtools_LIB)
petagstats_LDADD = -lbam -lz -lpthread
petagstats_SOURCES = petagstats.cpp
