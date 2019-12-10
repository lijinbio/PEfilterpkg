bin_PROGRAMS = pefilter

samtools_INCLUDE = $(top_srcdir)/lib/samtools-0.1.20
samtools_LIB = $(top_srcdir)/lib/samtools-0.1.20

CXXFLAGS = -g -O3 -std=c++11 -static
pefilter_CPPFLAGS = -Wall -w -I$(samtools_INCLUDE)
pefilter_LDFLAGS = -L$(samtools_LIB)
pefilter_LDADD = -lbam -lz -lpthread -lboost_program_options
pefilter_SOURCES = pefilter.cpp
