bin_PROGRAMS = pefilterpico

samtools_INCLUDE = $(top_srcdir)/lib/samtools-0.1.20
samtools_LIB = $(top_srcdir)/lib/samtools-0.1.20

CXXFLAGS = -g -O3 -std=c++11 -static
pefilterpico_CPPFLAGS = -Wall -w -I$(samtools_INCLUDE)
pefilterpico_LDFLAGS = -L$(samtools_LIB)
pefilterpico_LDADD = -lbam -lz -lpthread
pefilterpico_SOURCES = pefilterpico.cpp
