#include "reads_correction_aux.h"

void normalize_gaps(const char* qstr, const char* tstr, const index_t aln_size, std::string& qnorm, std::string& tnorm, const bool push)
{//normalize_gaps(m5qaln(*m5), m5saln(*m5), strlen(m5qaln(*m5)), nqstr, ntstr, true);
    qnorm.clear();
    tnorm.clear();
    const char kGap = '-';

#ifndef NDEBUG
    int qcnt = 0, tcnt = 0;
    for (index_t i = 0; i < aln_size; ++i)
    {
        const char qc = qstr[i];
        const char tc = tstr[i];
        if (qc != kGap) ++qcnt;
        if (tc != kGap) ++tcnt;
    }
#endif

    // convert mismatches to indels
    for (index_t i = 0; i < aln_size; ++i)
    {
        const char qc = qstr[i];
        const char tc = tstr[i];
        if (qc != tc && qc != kGap && tc != kGap)
        { 
            qnorm += kGap; qnorm += qc; tnorm += tc; tnorm += kGap;
            std::cout<<"这里有mismatch!!!!";
        }
        else
        { qnorm += qc; tnorm += tc; }
    }

    // push gaps to the right, but not pass the end
    if (push)
    {
        index_t qlen = qnorm.size();
        index_t tlen = tnorm.size();
        for (index_t i = 0; i < qlen - 1; ++i)
        {
            // push target gaps
            if (tnorm[i] == kGap)
            {
                index_t j = i;
                while (1)
                {
                    const char c = tnorm[++j];
                    if (c != kGap || j > qlen - 1)
                    {
                        if (c == qnorm[i]) { tnorm[i] = c; tnorm[j] = kGap; }
                        break;
                    }
                }
            }
            // push query gaps
            if (qnorm[i] == kGap)
            {
                index_t j = i;
                while (1)
                {
                    const char c = qnorm[++j];
                    if (c != kGap || j > tlen - 1)
                    {
                        if (c == tnorm[i]) { qnorm[i] = c; qnorm[j] = kGap; }
                        break;
                    }
                }
            }
        }
    }
    r_assert(qnorm.size() == tnorm.size());

#ifndef NDEBUG
    int qcnt2 = 0, tcnt2 = 0;
    for (std::string::const_iterator citer = qnorm.begin(); citer != qnorm.end(); ++citer)
        if ((*citer) != kGap) ++qcnt2;
    for (std::string::const_iterator citer = tnorm.begin(); citer != tnorm.end(); ++citer)
        if ((*citer) != kGap) ++tcnt2;
    d_assert(qcnt == qcnt2);
    d_assert(tcnt == tcnt2);
#endif
}

void slide_window(std::string& str1, std::string& str2, const int k){//slide windows strategy
    if(str1.size() != str2.size()){
        return;
    }
    int len = str1.size();
    int same_base = 0;
    double identity = 0;
    std::string window1 = str1.substr(0, k);
    std::string window2 = str2.substr(0, k);
    for(int i = 0; i < k; ++i){
        if(window1.substr(i, 1) == window2.substr(i, 1)){
            ++same_base;
        }
    }
    identity = (double)same_base/k;
    for(int index = 1; index <= len - k; ++index){
        if(window1.substr(0,1) == window2.substr(0,1)){
            --same_base;
        }
        window1 = str1.substr(index, k);;
        window2 = str2.substr(index, k);
        if(window1.substr(k-1, 1) == window2.substr(k-1, 1)){
            ++same_base;
        }
        identity = (double)same_base/k;
    }
    return;
}

void 
slide_window2(const std::string& str1, const std::string& str2, std::string& newstr1, const int k, const double identity_threshold){//slide windows strategy
    if(str1.size() != str2.size()){
        return;
    }
    int len = str1.size();
    int same_base = 0;
    double identity = 0;
    
    newstr1 = str1;
    
    std::string window1 = str1.substr(0, k);
    std::string window2 = str2.substr(0, k);
    for(int i = 0; i < k; ++i){
        if(window1.substr(i, 1) == window2.substr(i, 1)){
            ++same_base;
        }
    }
    identity = (double)same_base/k;
    if(identity < identity_threshold){
        newstr1[0] = 'N';
    }
    for(int index = 1; index <= len - k; ++index){
        if(window1.substr(0,1) == window2.substr(0,1)){
            --same_base;
        }
        window1 = str1.substr(index, k);
        window2 = str2.substr(index, k);
        if(window1.substr(k-1, 1) == window2.substr(k-1, 1)){
            ++same_base;
        }
        identity = (double)same_base/k;
        if(identity < identity_threshold){
            newstr1[index] = 'N';
        }
    }
    int newk = k;
    for(int index = len - newk + 1; index < len && newk > 0; ++index){
        if(window1.substr(0,1) == window2.substr(0,1)){
            --same_base;
        }
        --newk;
        window1 = str1.substr(index, newk);
        window2 = str2.substr(index, newk);
        identity = (double)same_base/newk;
        if(identity < identity_threshold){
            newstr1[index] = 'N';
        }
    }
    return;
}








struct CmpExtensionCandidateBySid
{
	bool operator()(const ExtensionCandidate& a, const ExtensionCandidate& b)
	{
		return a.sid < b.sid;
	}
};

void
build_cns_thrd_data_can(ExtensionCandidate* ec_list, 
						const int nec,
						const idx_t min_rid,
						const idx_t max_rid,
						ReadsCorrectionOptions* prco,
						PackedDB* reads,
						std::ostream* out,
					    ConsensusThreadData** ppctd)
{//build_cns_thrd_data_can(ec_list, num_ec, min_read_id, max_read_id, &rco, &reads, &out, pctds);
	const index_t num_reads = max_rid - min_rid + 1;
	const int num_threads = prco->num_threads;
    const index_t num_reads_per_thread = (num_reads + num_threads - 1) / num_threads;
	std::sort(ec_list, ec_list + nec, CmpExtensionCandidateBySid());
	idx_t max_id = min_rid;
	idx_t i = 0, j;
	int tid = 0;
	while (i < nec)
	{
		max_id += num_reads_per_thread;
		j = i + 1;
		while (j < nec && ec_list[j].sid < max_id) ++j;
		ppctd[tid] = new ConsensusThreadData(prco, tid, reads, ec_list + i, j - i, out);
		++tid;
		i = j;
	}
}
