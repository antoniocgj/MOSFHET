#include "spqlios-fft.h"

int32_t rev(int32_t x, int32_t M) {
    int32_t reps = 0;
    for (int32_t j = M; j > 1; j /= 2) {
        reps = 2 * reps + (x % 2);
        x /= 2;
    }
    return reps;
}

FFT_Processor_Spqlios new_FFT_Processor_Spqlios(const int32_t N){
    FFT_Processor_Spqlios res;
    res = (FFT_Processor_Spqlios) safe_malloc(sizeof(*res));
    res->_2N = 2*N;
    res->N = N;
    res->Ns2 = N/2;
    res->tables_direct = new_fft_table(N);
    res->tables_reverse = new_ifft_table(N);
    res->real_inout_direct = fft_table_get_buffer(res->tables_direct);
    res->real_inout_rev = fft_table_get_buffer(res->tables_reverse);
    return res;
}

FFT_Processor_Spqlios copy_FFT_Processor_Spqlios(FFT_Processor_Spqlios proc){
    FFT_Processor_Spqlios res;
    res = (FFT_Processor_Spqlios) safe_malloc(sizeof(*res));
    memcpy(res, proc, sizeof(*proc));
    res->real_inout_direct = (double *) safe_aligned_malloc(proc->N*sizeof(double));
    res->real_inout_rev = (double *) safe_aligned_malloc(proc->N*sizeof(double));
    return res;
}

void execute_reverse_int(double *res, const int32_t *a, FFT_Processor_Spqlios proc) {
    //for (int32_t i=0; i<N; i++) real_inout_rev[i]=(double)a[i];
    {
        double *dst = proc->real_inout_rev;
        const int32_t *ait = a;
        const int32_t *aend = a + proc->N;
        __asm__ __volatile__ (
        "0:\n"
                "vmovupd (%1),%%xmm0\n"
                "vcvtdq2pd %%xmm0,%%ymm1\n"
                "vmovapd %%ymm1,(%0)\n"
                "addq $16,%1\n"
                "addq $32,%0\n"
                "cmpq %2,%1\n"
                "jb 0b\n"
        : "=r"(dst), "=r"(ait), "=r"(aend)
        : "0"(dst), "1"(ait), "2"(aend)
        : "%xmm0", "%ymm1", "memory"
        );
    }
    ifft(proc->tables_reverse, proc->real_inout_rev);
    //for (int32_t i=0; i<N; i++) res[i]=real_inout_rev[i];
    {
        double *dst = res;
        double *sit = proc->real_inout_rev;
        double *send = proc->real_inout_rev + proc->N;
        __asm__ __volatile__ (
        "1:\n"
                "vmovapd (%1),%%ymm0\n"
                "vmovupd %%ymm0,(%0)\n"
                "addq $32,%1\n"
                "addq $32,%0\n"
                "cmpq %2,%1\n"
                "jb 1b\n"
                "vzeroall\n"
        : "=r"(dst), "=r"(sit), "=r"(send)
        : "0"(dst), "1"(sit), "2"(send)
        : "%ymm0", "memory"
        );
    }
}

void execute_reverse_torus32(double *res, const uint32_t *a, FFT_Processor_Spqlios proc) {
    int32_t *aa = (int32_t *) a;
    execute_reverse_int(res, aa, proc);
}

void execute_reverse_torus64(double* res, const uint64_t* a, FFT_Processor_Spqlios proc) {
    #ifdef AVX512_OPT
    __m512d * ri512 = (__m512d *) proc->real_inout_rev;
    __m512i * aa = (__m512i *) a;
    for (size_t i = 0; i < proc->N/8; i++) ri512[i] = _mm512_cvtepi64_pd (aa[i]);
    #else
    int64_t *aa = (int64_t *)a;
    for (int i=0; i<proc->N; i++) proc->real_inout_rev[i]=(double)aa[i];
    #endif
    ifft(proc->tables_reverse,proc->real_inout_rev);
    for (int i=0; i<proc->N; i++) res[i]=proc->real_inout_rev[i];
}

void execute_direct_torus32(uint32_t *res, const double *a, FFT_Processor_Spqlios proc) {
    //TODO: parallelization
    double _2sN = ((double) 2) / ((double) proc->N);
    //for (int32_t i=0; i<N; i++) real_inout_direct[i]=a[i]*_2sn;
    {
        double *dst = proc->real_inout_direct;
        const double *sit = a;
        const double *send = a + proc->N;
        //double __2sN = 2./N;
        const double *bla = &_2sN;
        __asm__ __volatile__ (
        "vbroadcastsd (%3),%%ymm2\n"
                "1:\n"
                "vmovupd (%1),%%ymm0\n"
                "vmulpd	%%ymm2,%%ymm0,%%ymm0\n"
                "vmovapd %%ymm0,(%0)\n"
                "addq $32,%1\n"
                "addq $32,%0\n"
                "cmpq %2,%1\n"
                "jb 1b\n"
        : "=r"(dst), "=r"(sit), "=r"(send), "=r"(bla)
        : "0"(dst), "1"(sit), "2"(send), "3"(bla)
        : "%ymm0", "%ymm2", "memory"
        );
    }
    fft(proc->tables_direct, proc->real_inout_direct);
    for (int32_t i = 0; i < proc->N; i++) res[i] = (uint32_t)((int64_t) proc->real_inout_direct[i]);
}

void execute_direct_torus64(uint64_t* res, const double* a, FFT_Processor_Spqlios proc) {
    double _2sN = ((double) 2) / ((double) proc->N);
    //static const double _2p64 = pow(2.,64);
    //for (int i=0; i<N; i++) real_inout_direct[i]=a[i]*_2sn;
    {
    double* dst = proc->real_inout_direct;
	const double* sit = a;
	const double* send = a+proc->N;
	//double __2sN = 2./N;
	const double* bla = &_2sN;
	__asm__ __volatile__ (
		"vbroadcastsd (%3),%%ymm2\n"
		"1:\n"
		"vmovupd (%1),%%ymm0\n"
		"vmulpd	%%ymm2,%%ymm0,%%ymm0\n"
		"vmovapd %%ymm0,(%0)\n"
		"addq $32,%1\n"
		"addq $32,%0\n"
		"cmpq %2,%1\n"
		"jb 1b\n"
		: "=r"(dst),"=r"(sit),"=r"(send),"=r"(bla)
		: "0"(dst),"1"(sit),"2"(send),"3"(bla)
		: "%ymm0","%ymm2","memory"
		);
    }
    fft(proc->tables_direct,proc->real_inout_direct); 
    // mod 2^64
    #ifdef AVX512_OPT
    __m512d * ri512 = (__m512d *) proc->real_inout_direct;
    __m512i * res512 = (__m512i *) res;
    const __m512d modc = {64, 64, 64, 64, 64, 64, 64, 64};
    for (size_t i = 0; i < proc->N/8; i++) {
        const __m512d _1 = _mm512_scalef_pd (ri512[i], -modc);
        const __m512d _2 = _mm512_reduce_pd (_1, 0);
        const __m512d _3 = _mm512_scalef_pd (_2, modc);
        res512[i] = _mm512_cvtpd_epi64 (_3);
    }
    #else
    const uint64_t* const vals = (const uint64_t*) proc->real_inout_direct;
    static const uint64_t valmask0 = 0x000FFFFFFFFFFFFFul;
    static const uint64_t valmask1 = 0x0010000000000000ul;
    static const uint16_t expmask0 = 0x07FFu;
    for (int i=0; i<proc->N; i++) {
        uint64_t val = (vals[i]&valmask0)|valmask1; //mantissa on 53 bits
        uint16_t expo = (vals[i]>>52)&expmask0; //exponent 11 bits
        // 1023 -> 52th pos -> 0th pos
        // 1075 -> 52th pos -> 52th pos
        int16_t trans = expo-1075;
        uint64_t val2 = trans>0?(val<<trans):(val>>-trans);
        res[i]=(vals[i]>>63)?-val2:val2;
    }
    #endif
}