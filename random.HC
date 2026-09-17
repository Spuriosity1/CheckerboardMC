
F64 to_01(U64 i){
    U64 mask1 = 0x3FF0000000000000;
    U64 mask2 = 0x3FFFFFFFFFFFFFFF;
    U64 to_12 = (i | mask1) & mask2;
    F64 d;
    MemCpy(&d, &to_12, 8);
    return d - 1;
}


class splitmix64 {
    U64 x;
};

class Xoshiro128 {
    U64 state[2];
};

splitmix64 NewSplitmix64(U64 seed) {
    splitmix64 sm;
    sm.x=seed;
    return sm;
}

U64 NextSM64(splitmix64* sm64) {
    auto z = (sm64->x += 0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EB;
    return z ^ (z >> 31);
}

Xoshiro128 NewXoshiro128(U64 seed) {
    auto sm = NewSplitmix64(seed);
    Xoshiro128 res;

    res.state[0] = NextSM64(&sm);
    res.state[1] = NextSM64(&sm);

    return res;
}

asm {
_ROTL::
    MOVQ     RAX,RDI
    MOVQ     RDX,RDI

    MOVQ     RCX,RSI

    SHLQ     RAX, CL
    NEGQ     RCX
    ADDQ     RCX, 64
    SHRQ     RDX, CL
    ORQ       RAX,RDX

    RET

_ROTR::
    MOVQ     RAX,RDI
    MOVQ     RDX,RDI

    MOVQ     RCX,RSI

    SHRQ     RAX, CL
    NEGQ     RCX
    ADDQ     RCX, 64
    SHLQ     RDX, CL
    ORQ       RAX,RDX

    RET
}


public _extern _ROTL U64 RotL(U64 x, U64 distance);
public _extern _ROTR U64 RotR(U64 x, U64 distance);


U64 NextXS128(Xoshiro128* x)
  {
    auto s0 = x->state[0];
    auto s1 = x->state[1];
    auto result = s0 + s1;

    s1 ^= s0;
    x->state[0] = RotL(s0, 24) ^ s1 ^ (s1 << 16);  // a, b
    x->state[1] = RotL(s1, 37);                    // c

    return result;
  }

