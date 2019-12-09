#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include "sam.h"

using namespace std;

// From samtools 0.1.19
// callback function for bam_fetch() that prints nonskipped records

// Global map for each chr, remember to clear for each chr
map< string, vector< string > > read2tag;
static int addtag(const bam1_t *b, void *data) {
	string qname=string((char*)bam1_qname(b));
	string zs=string((char *)bam_aux2Z(bam_aux_get(b, "ZS")));
	uint32_t flag=b->core.flag;

	// Skip the multiple mapping @ 20191125
	if (flag & 0x100) return 1;

	map< string, vector< string > > :: iterator it=read2tag.find(qname);
	if (read2tag.end()!=it) {
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
		read2tag[qname]=tag;
	}
	return 0;
}

// Six true PE mappings:
//   (++,+-)
//   (-+,--)
//   (++,N)
//   (-+,N)
//   (N,+-)
//   (N,--)
set< string > validtags {"++,+-", "-+,--", "++,N", "-+,N", "N,+-", "N,--"};
static int filter(const bam1_t *b, void *data) {
	string qname=string((char*)bam1_qname(b));
	map< string, vector< string > > :: iterator it=read2tag.find(qname);
	// skip multiple mapping in both ends
	if (read2tag.end()!=it) {
		string tags=it->second[0]+","+it->second[1];
		set< string > :: iterator sit=validtags.find(tags);
		if (validtags.end()!=sit) {
			samwrite((samfile_t*)data, b);
		}
	}
	return 0;
}

int main(int argc, char ** argv)
{
	string bamfile=argv[1];
	string outfile=argv[2];

	samfile_t *in=0;
	if ((in=samopen(bamfile.c_str(), "rb", 0))==0) {
		cerr << "Error: not found " << bamfile << endl;
		return 1;
	}

	vector< string> chroms;
	for (int i=0; i<in->header->n_targets; i++) {
		chroms.push_back(in->header->target_name[i]);
	}

	samfile_t *out=0;
	if ((out=samopen(outfile.c_str(), "wb", in->header))==0) {
		cerr << "Error: can not write " << outfile << endl;
		return 1;
	}

	bam_index_t *idx=0;
	idx = bam_index_load(bamfile.c_str());
	if (idx==0) {
		cerr << "Error: not found index file of " << bamfile << endl;
		return 1;
	}

	map< string, int > tagstats;
	for (string &chr : chroms) {
		int tid, beg, end, result;
		bam_parse_region(in->header, chr.c_str(), &tid, &beg, &end); // parse a region in the format like `chr2:100-200'
		if (tid<0) { 
			cerr << "Error: unknown reference name " << chr << endl;
			continue;
		}
		// 1. First scan to construct the tag directionary
		result=bam_fetch(in->x.bam, idx, tid, beg, end, NULL, addtag);
		if (result<0) {
			cerr << "Error: failed to retrieve region " << bamfile << endl;
			return 1;
		}
		// 2. Second scan to filter false paired mapping
		result=bam_fetch(in->x.bam, idx, tid, beg, end, out, filter);
		if (result<0) {
			cerr << "Error: failed to retrieve region " << bamfile << endl;
			return 1;
		}
		// 3. Record the tag statistics
		for (map< string, vector< string > > :: iterator it=read2tag.begin(); read2tag.end()!=it; ++it) {
			string tag=it->second[0] + "," + it->second[1];
			map< string, int > :: iterator tit=tagstats.find(tag);
			if (tagstats.end()!=tit) {
				tit->second++;
			} else {
				tagstats[tag]=1;
			}
		}
		read2tag.clear();
	}
	samclose(in);
	samclose(out);
	for (map< string, int > :: iterator it=tagstats.begin(); tagstats.end()!=it; ++it) {
		cout << it->first << "\t" << it->second << endl;
	}
	return 0;
}
