#ifndef _READS_CORRECTION_AUX_H
#define _READS_CORRECTION_AUX_H

#include <vector>
#include <cstring>

#include "dw.h"
#include "../common/packed_db.h"
#include "options.h"

struct CnsTableItem
{//修改的地方
	char base;
	uint1 mat_cnt;
	uint1 ins_cnt;
	uint1 del_cnt;
	uint1 skip_cnt;
	
	CnsTableItem() : base('N'), mat_cnt(0), skip_cnt(0), ins_cnt(0), del_cnt(0) {}
};

struct CnsTableItemCleaner
{
	void operator()(CnsTableItem& item) 
	{
		item.base = 'N';
		item.mat_cnt = 0;
		item.ins_cnt = 0;
		item.del_cnt = 0;
		item.skip_cnt = 0;
	}
};

#define MAX_CNS_OVLPS 100

struct MappingRange
{
	int start, end;
	
	MappingRange(int s, int e) : start(s), end(e) {}
};

struct CnsAln
{
	int soff, send, aln_idx, aln_size;
	char qaln[MAX_SEQ_SIZE];
	char saln[MAX_SEQ_SIZE];
	
	bool retrieve_aln_subseqs(int sb, int se, std::string& qstr, std::string& tstr, int& sb_out)
	{//首先判断是否有错误，很简单的阈值判断
		//其次把复杂区域copy到qstr,tstr
		if (se <= soff || sb >= send || aln_idx >= aln_size - 1) return false;
		sb_out = std::max(soff, sb);
		qstr.clear();
		tstr.clear();
		while(soff < sb && aln_idx < aln_size - 1)
		{
			++aln_idx;
			if (saln[aln_idx] != GAP) ++soff;
		}
		if(qaln[aln_idx] == 'N'){
			return false;
		}
		qstr += qaln[aln_idx];
		tstr += saln[aln_idx];
		while (soff < se && aln_idx < aln_size - 1)
		{
			++aln_idx;
			if (saln[aln_idx] != GAP) ++soff;
			if(qaln[aln_idx] == 'N'){
				return false;
			}
			qstr += qaln[aln_idx];
			tstr += saln[aln_idx];
		}
		return true;
	}
};

class CnsAlns
{
public:
	CnsAlns()
	{
		//std::cout<<"before new a CnsAlns.:)))."<<"\n"<<std::endl;
		safe_malloc(cns_alns_, CnsAln, MAX_CNS_OVLPS);
		//std::cout<<"after new a CnsAlns.:)))."<<"\n"<<std::endl;
		clear();
		//std::cout<<"initiall finished."<<"\n"<<std::endl;
	}
	~CnsAlns()
	{
		//std::cout<<"before safe_free."<<"\n"<<std::endl;
		safe_free(cns_alns_);
		//std::cout<<"safe_free end."<<"\n"<<std::endl;
	}
	void clear() { num_alns_ = 0; }
	int num_alns() { return num_alns_; }
	CnsAln* begin() { return cns_alns_; }
	CnsAln* end() { return cns_alns_ + num_alns_; }
	void add_aln(const int soff, const int send, const std::string& qstr, const std::string& tstr)
	{//cns_vec.add_aln(m5soff(*m5), m5send(*m5), nqstr, ntstr);
		//std::cout<<"MAX_CNS_OVLPS="<<MAX_CNS_OVLPS<<"&&&&num_alns_="<<num_alns_<<"\n";
		//std::cout<<"enter add_aln."<<"\n";
		r_assert(qstr.size() == tstr.size());
		//std::cout<<"r_assert end."<<"\n";

		CnsAln& a = cns_alns_[num_alns_++];
		//std::cout<<"1 end."<<"\n";
		//std::cout<<"soff="<<soff<<"\n";
		//std::cout<<"num_alns_="<<num_alns_<<"\n"<<std::endl;
		//std::cout<<"a="<<&a<<"\n"<<std::endl;
		//std::cout<<"a="<<&a.soff<<"\n"<<std::endl;
		//std::cout<<"*cns_alns_="<<*cns_alns_<<"\n"<<std::endl;
		//std::cout<<"cns_alns_[num_alns_]="<<int(cns_alns_[num_alns_])<<"\n"<<std::endl;
		//std::cout<<"&cns_alns_[num_alns_]="<<&cns_alns_[num_alns_]<<"\n"<<std::endl;
		//std::cout<<"*cns_alns_="<<*cns_alns_<<"\n"<<std::endl;
		//std::cout<<"a.soff="<<a.soff<<std::endl;//这里也没有输出
		a.soff = soff;//断在这里了
		//std::cout<<"1_1 end."<<"\n";
		a.send = send;
		//std::cout<<"1_2 end."<<"\n";
		a.aln_idx = 0;
		//std::cout<<"1_3 end."<<"\n";
		a.aln_size = qstr.size();
		//std::cout<<"2 end."<<"\n";
		memcpy(a.qaln, qstr.data(), qstr.size());
		//std::cout<<"3 end."<<"\n";
		a.qaln[qstr.size()] = '\0';
		//std::cout<<"4 end."<<"\n";
		memcpy(a.saln, tstr.data(), tstr.size());
		//std::cout<<"5 end."<<"\n";
		a.saln[tstr.size()] = '\0';
		//std::cout<<"6 end."<<"\n";
	}
	void get_mapping_ranges(std::vector<MappingRange>& ranges)
	{
		ranges.clear();
		for (int i = 0; i < num_alns_; ++i) ranges.push_back(MappingRange(cns_alns_[i].soff, cns_alns_[i].send));
	}
	
private:
	CnsAln* cns_alns_;
	int     num_alns_;
};

#define MAX_CNS_RESULTS 10000

struct ConsensusThreadData
{
	ReadsCorrectionOptions rco;
	int thread_id;
	PackedDB* reads;
	ExtensionCandidate* candidates;
	int num_candidates;
	ns_banded_sw::DiffRunningData* drd_s;
	ns_banded_sw::DiffRunningData* drd_l;
	M5Record* m5;
	CnsAlns cns_alns;
	//加入一个东西
	CnsAlns new_cns_alns;
	std::vector<CnsResult> cns_results;
	std::vector<char> query;
	std::vector<char> target;
	std::string qaln;
	std::string saln;
	CnsTableItem* cns_table;
	uint1* id_list;
	std::ostream* out;
	pthread_mutex_t out_lock;
	
	ConsensusThreadData(ReadsCorrectionOptions* prco, int tid, PackedDB* r, ExtensionCandidate* ec, int nec, std::ostream* output)
	{//ppctd[tid] = new ConsensusThreadData(prco, tid, reads, ec_list + i, j - i, out);
		rco = (*prco);
		thread_id = tid;
		reads = r;
		candidates = ec;
		num_candidates = nec;
		drd_s = new ns_banded_sw::DiffRunningData(ns_banded_sw::get_sw_parameters_small());
		drd_l = new ns_banded_sw::DiffRunningData(ns_banded_sw::get_sw_parameters_large());
		m5 = NewM5Record(MAX_SEQ_SIZE);
		out = output;
		
		query.reserve(MAX_SEQ_SIZE);
		target.reserve(MAX_SEQ_SIZE);
		qaln.reserve(MAX_SEQ_SIZE);
		saln.reserve(MAX_SEQ_SIZE);
		safe_malloc(cns_table, CnsTableItem, MAX_SEQ_SIZE);
		safe_malloc(id_list, uint1, MAX_SEQ_SIZE);
		pthread_mutex_init(&out_lock, NULL);
	}
	
	~ConsensusThreadData()
	{
		delete drd_s;
		delete drd_l;
		m5 = DeleteM5Record(m5);
		safe_free(cns_table);
		safe_free(id_list);
	}
};

void normalize_gaps(const char* qstr, const char* tstr, const index_t aln_size, std::string& qnorm, std::string& tnorm, const bool push);

void slide_window(std::string& str1, std::string& str2, const int k);
void slide_window2(const std::string& str1, const std::string& str2, std::string& newstr1, const int k, const double identity_threshold);

void
build_cns_thrd_data_can(ExtensionCandidate* ec_list, 
						const int nec,
						const idx_t min_rid,
						const idx_t max_rid,
						ReadsCorrectionOptions* prco,
						PackedDB* reads,
						std::ostream* out,
					    ConsensusThreadData** ppctd);

#endif // _READS_CORRECTION_AUX_H
