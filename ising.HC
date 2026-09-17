#include "random.HC"

F64 Atod(U8* nptr) {
    U8* end;
    return strtod(nptr, &end); // todo error handling...
}

class Settings {
    F64 T_hot;
    F64 T_cold;
    F64 J;
    U64 n_step;
    U64 n_sweep;
    U64 n_sample;
    splitmix64 rng;
};

class Spin {
    I8 S; // Ising spin state
    Spin** neighbours;
    I32 n_neighbours;
};


U32 lindex_at(I32 i, I32 j, U64 L) {
    return L*((i+L)%L) + (j+L)%L;
}

U0 assign_neighbours(Spin* spins, U64 L) {
    for (I32 i=0; i<L; i++){
    for (I32 j=0; j<L; j++){
        auto nn = 6;
        spins[L*i+j].n_neighbours=nn;
        spins[L*i+j].neighbours = MAlloc(sizeof(Spin**)*nn);

        // The horizontals
        spins[L*i+j].neighbours[0] = &(spins[lindex_at(i+1, j, L)]);
        spins[L*i+j].neighbours[1] = &(spins[lindex_at(i, j+1, L)]);

        spins[L*i+j].neighbours[3] = &(spins[lindex_at(i-1, j, L)]);
        spins[L*i+j].neighbours[4] = &(spins[lindex_at(i, j-1, L)]);

        // the diagonals (NB only sane if L is even!)
        I32 sl = ((i+j)%2)*2 -1; 
        spins[L*i+j].neighbours[2] = &(spins[lindex_at(i+1, j+sl, L)]);
        spins[L*i+j].neighbours[5] = &(spins[lindex_at(i-1, j-sl, L)]);
    }
    }
}

U0 free_neighbours(Spin* spins, U64 L){
    for (U64 i=0; i<L*L; i++){
        Free(spins[i].neighbours);
    }
}

F64 local_field(Spin* s){
    I64 h=0;
    for (U64 i=0; i<s->n_neighbours; i++){
        h += s->neighbours[i]->S;
    }
    return h (F64);
}

F64 energy(Spin* spins, U64 L, Settings* set){
    F64 acc=0;
    I64 tmp;

    for (I32 i=0; i<L; i++){
    for (I32 j=0; j<L; j++){
        auto s=spins[L*i+j];
        tmp=0;
        // include exactly half the bonds (note convention chosen)
        for (U32 i=0; i<s.n_neighbours/2; i++)
            tmp += s.neighbours[i]->S;

        tmp *= s.S;
        acc += tmp (F64);
    }}
    return acc * set->J;

}

U0 MetropolisSweep(Spin* spins, U64 L, Settings* set, F64 beta){
    for (I32 i=0; i<L*L; i++){
        auto s = &spins[i];

        F64 dE = -2*local_field(s) * s->S (F64) * set->J ;
        F64 r = to_01(NextSM64(&set->rng));

        if ( r< exp( -dE * beta )) {
            s->S *= -1;
        }   
    }
}


U8 Setup(I32 argc, U8 **argv, Settings* set, U64* L)
{
    auto help=StrNew("Call with ising <L> <Thot> <Tcold> <seed> <n_step=100> <n_sweep=100> <n_sample=30> <J=1>\n");
    U8 status=1;

    if (argc < 5) goto error;

    *L = 2*Atoi(argv[1]);
    set->T_hot = Atod(argv[2]);
    set->T_cold = Atod(argv[3]);
    set->rng = NewSplitmix64(Atoi(argv[4]));

    set->n_step=100;
    set->n_sweep=100;
    set->n_sample=30;
    set->J=1;

    switch (argc) {
        case 9: set->J = Atod(argv[8]);
        case 8: set->n_sample = Atoi(argv[7]);
        case 7: set->n_sweep = Atoi(argv[6]);
        case 6: set->n_step = Atoi(argv[5]);
        case 5:
            break;
        default:
            goto error;
    }
    status=0;
    goto done;
    
error:
    "%s",help;
done:
    Free(help);
    return status;
}

// does all initialisation work, leaving the system ready for 
// application of a MetropolisSweep
U0 Init(Spin* spins, U64 L)
{
    
    // initialisation: FM state
    for (I32 i=0; i<L*L; i++) {
        spins[i].S = 1;
    }
    // assign neighbours
    assign_neighbours(spins, L);
}


U8 Main(I32 argc, U8 **argv)
{
    Settings set;
    U64 L;
    if (Setup(argc, argv, &set, &L)) return 1;

    Spin* spins = MAlloc(L * L * sizeof(Spin));

    Init(spins, L);


    F64* E = MAlloc(set.n_step * sizeof(F64));
    F64* E2 = MAlloc(set.n_step * sizeof(F64));
    F64* T = MAlloc(set.n_step * sizeof(F64));


    F64 beta = 1./set.T_hot;
    F64 factor = exp((log(set.T_hot)-log(set.T_cold)) / (set.n_step-1));

    // the sweep loop

    F64 E_acc, E2_acc, tmp;


    for (U64 step=0; step < set.n_step; step++) {
        // burn for sweep_no steps
        for (U64 sweep_no=0; sweep_no < set.n_sweep; sweep_no++) {
            MetropolisSweep(spins, L, &set, beta);
        }

        T[step]=1./beta;

        E[step]=0;
        E2[step]=0;

        for (U64 sweep_no=0; sweep_no < set.n_sample; sweep_no++){
            MetropolisSweep(spins, L, &set, beta);
            tmp = energy(spins, L, &set);
            E[step] += tmp;
            E2[step] += tmp*tmp;
        }

        E[step]/= set.n_sample;
        E2[step]/= set.n_sample;

        "%e\t%e\t%e\n",1./beta,E[step], beta*beta*(E2[step]-E[step]*E[step]) / (L*L);

        beta *= factor;
    }

    free_neighbours(spins, L);
    Free(spins);

    return 0;
}
