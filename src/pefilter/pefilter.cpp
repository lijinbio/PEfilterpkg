#include <boost/program_options.hpp>
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <iterator>
#include <thread>
#include <cstdio>
#include "sam.h"

using namespace boost::program_options;
using namespace std;

class Opts {
	public:
		string infile;
		string outfile;
		bool pico;
		bool statsonly;
		int numthreads;
	public:
		Opts():
			infile("")
			, outfile("")
			, pico(false)
			, statsonly(false)
			, numthreads(1) { }
	public:
		void out() {
			cout << "infile: " << infile << endl;
			cout << "outfile: " << outfile << endl;
			cout << "pico: " << std::boolalpha << pico << endl;
			cout << "statsonly: " << std::boolalpha << statsonly << endl;
			cout << "numthreads: " << numthreads << endl;
		}
} opts;

int parse_options(int ac, const char ** av) {
	try
	{
		options_description desc{"Allowed options"};
		desc.add_options()
			("help,h", "Produce help message. Example command:\npefilter -i in.bam -o out.bam\npefilter -i in.bam -p -s")
			("infile,i", value<string>()->default_value(""), "Input BAM file. It should be indexed.")
			("outfile,o", value<string>()->default_value(""), "Output BAM file. To save the filtered BAM file.")
			("pico,p", "Pico library preparation protocol. Default: traditional protocol.")
			("statsonly,s", "Report PE tag statistics only but not generate filtered BAM file. The statitics will show in stdout.")
			("numthreads,t", value<int>()->default_value(1), "Number of threads. Ensure enough memory for many threads. Default: 1.")
			;

		variables_map vm;
		store(parse_command_line(ac, av, desc), vm);
		notify(vm);

		if (vm.count("help")) {
			cout << desc << endl;
			exit(1);
		}

		for(map<string, variable_value>::iterator it=vm.begin(); it!=vm.end(); ++it) {
			string k=it->first;
			if( k == "infile"){
				opts.infile=vm[k].as<string>();
			} else if( k == "outfile"){
				opts.outfile=vm[k].as<string>();
			} else if( k == "numthreads"){
				opts.numthreads=vm[k].as<int>();
			} else if( k == "pico"){
				opts.pico=true;
			} else if( k == "statsonly"){
				opts.statsonly=true;
			} else {
				cerr << "Error: invalid option " << k << endl;
				exit(1);
			}
		}
		if (opts.infile.empty()) {
			cerr << "Error: -i|--infile must be specified." << endl;
			cout << desc << endl;
			exit(1);
		}
		if (opts.outfile.empty() && !opts.statsonly) {
			cerr << "Error: -o|--outfile must be specified." << endl;
			cout << desc << endl;
			exit(1);
		}
		opts.out();
	} catch (const error &ex) {
		cerr << ex.what() << endl;
	}
	return 0;
}


// From samtools 0.1.19
// callback function for bam_fetch() that prints nonskipped records

// Global map for each chr, remember to clear for each chr
map< int, map< string, vector< string > > > read2tag; // chrid->qname->[tag1,tag2]
static int addtag(const bam1_t *b, void *data) {
	uint32_t flag=b->core.flag;
	// Skip the multiple mapping @ 20191125
	if (flag & 0x100) return 1;

	string qname=string((char*)bam1_qname(b));
	string zs=string((char *)bam_aux2Z(bam_aux_get(b, "ZS")));
	uint32_t tid=b->core.tid;
	map< string, vector< string > > &read2tagchr=read2tag[tid]; // qname->[tag1,tag2]
	map< string, vector< string > > :: iterator it=read2tagchr.find(qname);
	if (read2tagchr.end()!=it) {
		if (flag & 0x40) {
			it->second[0]=zs;
		} else if (flag & 0x80) {
			it->second[1]=zs;
		}
	} else {
		vector< string > tag(2, "N");
		if (flag & 0x40) {
			tag[0]=zs;
		} else if (flag & 0x80) {
			tag[1]=zs;
		}
		read2tagchr[qname]=tag;
	}
	return 0;
}

map< string, map< string, int > > tagstats; // chr->tag->number
void petagstatschrbatch(string bamfile, vector< string > chrs) {
	samfile_t *in=0;
	if ((in=samopen(bamfile.c_str(), "rb", 0))==0) {
		cerr << "Error: not found " << bamfile << endl;
		return;
	}

	map< string, int > chr2tid;
	for (int i=0; i<in->header->n_targets; i++) {
		chr2tid[in->header->target_name[i]]=i;
	}

	bam_index_t *idx=0;
	idx = bam_index_load(bamfile.c_str());
	if (idx==0) {
		cerr << "Error: not found index file of " << bamfile << endl;
		return;
	}
	for (string &chr : chrs) {
		cout << "Start chromosome " << chr << endl;
		int tid, beg, end, result;
		bam_parse_region(in->header, chr.c_str(), &tid, &beg, &end);
		if (tid<0) { 
			cerr << "Error: unknown reference name " << chr << endl;
			return;
		}
		result=bam_fetch(in->x.bam, idx, tid, beg, end, NULL, addtag);
		if (result<0) {
			cerr << "Error: failed to retrieve region " << chr << endl;
			return;
		}

		map< string, int > &tagstatschr=tagstats[chr];
		map< string, vector< string > > &read2tagchr=read2tag[chr2tid[chr]];
		for (map< string, vector< string > > :: iterator it=read2tagchr.begin(); read2tagchr.end()!=it; ++it) {
			string tag=it->second[0] + "," + it->second[1];
			map< string, int > :: iterator tit=tagstatschr.find(tag);
			if (tagstatschr.end()!=tit) {
				tit->second++;
			} else {
				tagstatschr[tag]=1;
			}
		}
		read2tagchr.clear();
		cout << "End chromosome " << chr << endl;
	}
	samclose(in);
	bam_index_destroy(idx);
}

int petagstats(string bamfile)
{
	samfile_t *in=0;
	if ((in=samopen(bamfile.c_str(), "rb", 0))==0) {
		cerr << "Error: not found " << bamfile << endl;
		return 1;
	}
	vector< string > chroms;
	for (int i=0; i<in->header->n_targets; i++) {
		chroms.push_back(in->header->target_name[i]);
	}
	samclose(in);

	vector< vector< string > > chrbatch;
	for (int i=0; i<chroms.size(); i++) {
		if (i>=opts.numthreads) {
			chrbatch[i%opts.numthreads].push_back(chroms[i]);
		} else {
			vector< string > chrs {chroms[i]};
			chrbatch.push_back(chrs);
		}
	}

	vector<thread> threads;
	for (vector< string > &chrs : chrbatch) {
		threads.push_back(thread(petagstatschrbatch, bamfile, chrs));
	}
	for (auto& th : threads) {
		th.join();
	}

	map< string, int > tagsresult;
	for (map< string, map< string, int > > :: iterator itchr=tagstats.begin(); tagstats.end()!=itchr; ++itchr) {
		map< string, int > & tagstatschr=itchr->second;
		for (map< string, int > :: iterator it=tagstatschr.begin(); tagstatschr.end()!=it; ++it) {
			tagsresult[it->first]+=it->second;
		}
	}
	for (map< string, int > :: iterator it=tagsresult.begin(); tagsresult.end()!=it; ++it) {
		cout << it->first << "\t" << it->second << endl;
	}
	return 0;
}

// Six true PE mappings in traditional library preparation:
set< string > validtags_trad {
	"++,+-", "-+,--"
		, "++,N", "N,+-"
		, "-+,N", "N,--"
};
static int filter_trad(const bam1_t *b, void *data) {
	string qname=string((char*)bam1_qname(b));
	uint32_t tid=b->core.tid;
	map< string, vector< string > > &read2tagchr=read2tag[tid];
	map< string, vector< string > > :: iterator it=read2tagchr.find(qname);
	// skip multiple mapping in both ends
	if (read2tagchr.end()!=it) {
		string tags=it->second[0]+","+it->second[1];
		set< string > :: iterator sit=validtags_trad.find(tags);
		if (validtags_trad.end()!=sit) {
			samwrite((samfile_t*)data, b);
		}
	}
	return 0;
}
// 12 true PE mappings in Pico library preparation:
set< string > validtags_pico {
	"++,+-", "+-,++", "-+,--", "--,-+"
		, "++,N", "N,++", "+-,N", "N,+-"
		, "-+,N", "N,-+", "--,N", "N,--"
};
static int filter_pico(const bam1_t *b, void *data) {
	string qname=string((char*)bam1_qname(b));
	uint32_t tid=b->core.tid;
	map< string, vector< string > > &read2tagchr=read2tag[tid];
	map< string, vector< string > > :: iterator it=read2tagchr.find(qname);
	// skip multiple mapping in both ends
	if (read2tagchr.end()!=it) {
		string tags=it->second[0]+","+it->second[1];
		set< string > :: iterator sit=validtags_pico.find(tags);
		if (validtags_pico.end()!=sit) {
			samwrite((samfile_t*)data, b);
		}
	}
	return 0;
}


void pefilterchrbatch(string bamfile, string outfile, vector< string > chrs) {
	samfile_t *in=0;
	if ((in=samopen(bamfile.c_str(), "rb", 0))==0) {
		cerr << "Error: not found " << bamfile << endl;
		return;
	}

	map< string, int > chr2tid;
	for (int i=0; i<in->header->n_targets; i++) {
		chr2tid[in->header->target_name[i]]=i;
	}

	bam_index_t *idx=0;
	idx = bam_index_load(bamfile.c_str());
	if (idx==0) {
		cerr << "Error: not found index file of " << bamfile << endl;
		return;
	}
	for (string &chr : chrs) {
		cout << "Start chromosome " << chr << endl;
		string chroutfile=outfile+"_"+chr+".bam";
		samfile_t *out=0;
		if ((out=samopen(chroutfile.c_str(), "wb", in->header))==0) {
			cerr << "Error: can not write " << chroutfile << endl;
			return;
		}
		int tid, beg, end, result;
		bam_parse_region(in->header, chr.c_str(), &tid, &beg, &end);
		if (tid<0) { 
			cerr << "Error: unknown reference name " << chr << endl;
			return;
		}
		// 1. First scan to construct the tag directionary
		result=bam_fetch(in->x.bam, idx, tid, beg, end, NULL, addtag);
		if (result<0) {
			cerr << "Error: failed to retrieve region " << chr << endl;
			return;
		}
		// 2. Second scan to filter false paired mapping
		if (opts.pico) {
			result=bam_fetch(in->x.bam, idx, tid, beg, end, out, filter_pico);
		} else {
			result=bam_fetch(in->x.bam, idx, tid, beg, end, out, filter_trad);
		}
		if (result<0) {
			cerr << "Error: failed to filter region " << chr << endl;
			return;
		}
		samclose(out);
		// 3. Record the tag statistics
		map< string, int > &tagstatschr=tagstats[chr];
		map< string, vector< string > > &read2tagchr=read2tag[chr2tid[chr]];
		for (map< string, vector< string > > :: iterator it=read2tagchr.begin(); read2tagchr.end()!=it; ++it) {
			string tag=it->second[0] + "," + it->second[1];
			map< string, int > :: iterator tit=tagstatschr.find(tag);
			if (tagstatschr.end()!=tit) {
				tit->second++;
			} else {
				tagstatschr[tag]=1;
			}
		}
		read2tagchr.clear();
		cout << "End chromosome " << chr << endl;
	}
	samclose(in);
	bam_index_destroy(idx);
}

int mergebam(vector< string > & files, string & outfile) {
	string cmd = "samtools merge "+outfile;
	for (string &infile: files) {
		cmd += " " + infile;
	}
	cout << cmd << endl;
	FILE *fp;
	char info[10240];
	fp = popen(cmd.c_str(), "r");
	if (fp==NULL) {
		fprintf(stderr, "popen error.\n");
		return EXIT_FAILURE;
	}
	while (fgets(info, 10240, fp) != NULL) {
		printf("%s", info);
	}
	pclose(fp);
	return 0;
}

int rmtmpfiles(vector< string > & files) {
	string cmd = "rm -f";
	for (string &infile: files) {
		cmd += " " + infile;
	}
	cout << cmd << endl;
	FILE *fp;
	char info[10240];
	fp = popen(cmd.c_str(), "r");
	if (fp==NULL) {
		fprintf(stderr, "popen error.\n");
		return EXIT_FAILURE;
	}
	while (fgets(info, 10240, fp) != NULL) {
		printf("%s", info);
	}
	pclose(fp);
	return 0;
}

int pefilter(string bamfile, string outfile)
{
	samfile_t *in=0;
	if ((in=samopen(bamfile.c_str(), "rb", 0))==0) {
		cerr << "Error: not found " << bamfile << endl;
		return 1;
	}
	vector< string> chroms;
	for (int i=0; i<in->header->n_targets; i++) {
		chroms.push_back(in->header->target_name[i]);
	}
	samclose(in);

	vector< vector< string > > chrbatch;
	for (int i=0; i<chroms.size(); i++) {
		if (i>=opts.numthreads) {
			chrbatch[i%opts.numthreads].push_back(chroms[i]);
		} else {
			vector< string > chrs {chroms[i]};
			chrbatch.push_back(chrs);
		}
	}

	vector<thread> threads;
	for (vector< string > &chrs : chrbatch) {
		threads.push_back(thread(pefilterchrbatch, bamfile, outfile, chrs));
	}
	for (auto& th : threads) {
		th.join();
	}

	vector< string > tmpfiles;
	for (string &chr: chroms) {
		tmpfiles.push_back(outfile+"_"+chr+".bam");
	}
	mergebam(tmpfiles, outfile);
	rmtmpfiles(tmpfiles);

	map< string, int > tagsresult;
	for (map< string, map< string, int > > :: iterator itchr=tagstats.begin(); tagstats.end()!=itchr; ++itchr) {
		map< string, int > & tagstatschr=itchr->second;
		for (map< string, int > :: iterator it=tagstatschr.begin(); tagstatschr.end()!=it; ++it) {
			tagsresult[it->first]+=it->second;
		}
	}
	for (map< string, int > :: iterator it=tagsresult.begin(); tagsresult.end()!=it; ++it) {
		cout << it->first << "\t" << it->second << endl;
	}
	return 0;
}

int main(int argc, const char ** argv)
{
	parse_options(argc, argv);
	if (opts.statsonly) {
		petagstats(opts.infile);
	} else {
		pefilter(opts.infile, opts.outfile);
	}
	return 0;
}
