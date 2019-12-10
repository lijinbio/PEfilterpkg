bin_PROGRAMS = pefiltertag

samtools_INCLUDE = $(top_srcdir)/lib/samtools-0.1.20
samtools_LIB = $(top_srcdir)/lib/samtools-0.1.20

CXXFLAGS = -g -O3 -std=c++11 -static
pefiltertag_CPPFLAGS = -Wall -w -I$(samtools_INCLUDE)
pefiltertag_LDFLAGS = -L$(samtools_LIB)
pefiltertag_LDADD = -lbam -lz -lpthread -lboost_program_options
pefiltertag_SOURCES = pefiltertag.cpp
