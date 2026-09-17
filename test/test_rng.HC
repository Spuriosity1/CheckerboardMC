#include "../random.HC"


U8 Main(I32 argc, U8 **argv)
{
    if (argc < 3){
        "Usage: %s <num_repetitions> <seed>", argv[0];
        return 1;
    }
    auto i = 1;
    U64 count=Atoi(argv[1]);
    U64 seed  = Atoi(argv[2]);

    auto rng1 = NewSplitmix64(seed);
    auto rng2 = NewXoshiro128(seed);

    U32 N_CUMULANTS = 10;

    F64 *acc_sm64 = MAlloc(N_CUMULANTS*sizeof(F64));
    F64 *acc_xoshiro = MAlloc(N_CUMULANTS*sizeof(F64));
    
    for (U64 i=0; i<N_CUMULANTS; i++) {
        acc_sm64[i]=0.;
    }


    F64 tmp1;
    F64 tmp2;
    F64 r1, r2;

    for (I64 j=0; j<count; j++){
        r1 = to_01(NextSM64(&rng1));
        r2 = to_01(NextXS128(&rng1));

        tmp1=1.0;
        tmp2=1.0;
        for (I64 p=0; p<N_CUMULANTS; p++){
            tmp1 *= r1;
            acc_sm64[p] += tmp1;

            tmp2 *= r2;
            acc_xoshiro[p] += tmp2;
        }
        
    }


    "Cumulant\tanalytic\tsplitmix\terr_splitmix\txs128\terr_xs128\n";

    for (I64 p=0; p<N_CUMULANTS; p++){
        acc_sm64[p] /= count;
        acc_xoshiro[p] /= count;

        F64 analytic = 1.0 / (p+2) (F64);
        
        F64 err_sm64 = acc_sm64[p] - analytic;
        F64 err_xoshiro = acc_xoshiro[p] - analytic;

        "x^%d\t%f\t%f\t%e\t%f\t%e\n",p+1,analytic, acc_sm64[p], err_sm64, acc_xoshiro[p], err_xoshiro;
    }

    Free(acc_sm64);
    Free(acc_xoshiro);
    return 0;
}
