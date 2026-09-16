#include "random.HC"



U0 Main(I32 argc, U8 **argv)
{
    auto i = 1;
    U64 count=Atoi(argv[1]);

    auto rng = NewSplitmix64(55);

    F64 acc=0;
    F64 acc2=0;

    for (I64 j=0; j<count; j++){
        F64 r = to_01(NextSM64(&rng));
        acc += r;
        acc2 += r*r;
    }

    auto mu = acc/count;
    auto var = (acc2 / count - mu*mu);

    auto err_mu=mu-0.5;
    auto err_var=var- (1.0/12);

    "sum = %f, error = %e\n",mu,err_mu;
    "var = %f, error = %e\n",var,err_var;

    return;
}
