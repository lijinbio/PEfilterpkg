bin_PROGRAMS = pefilter

samtools_INCLUDE = $(top_srcdir)/lib/samtools-0.1.20
samtools_LIB = $(top_srcdir)/lib/samtools-0.1.20

CXXFLAGS = -g -O3
pefilter_CPPFLAGS = -Wall -w -I$(samtools_INCLUDE)
pefilter_LDFLAGS = -L$(samtools_LIB)
pefilter_LDADD = -lbam -lz -lpthread
pefilter_SOURCES = pefilter.cpp
