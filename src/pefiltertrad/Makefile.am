bin_PROGRAMS = pefiltertrad

samtools_INCLUDE = $(top_srcdir)/lib/samtools-0.1.20
samtools_LIB = $(top_srcdir)/lib/samtools-0.1.20

CXXFLAGS = -g -O3 -std=c++11
pefiltertrad_CPPFLAGS = -Wall -w -I$(samtools_INCLUDE)
pefiltertrad_LDFLAGS = -L$(samtools_LIB)
pefiltertrad_LDADD = -lbam -lz -lpthread
pefiltertrad_SOURCES = pefiltertrad.cpp
