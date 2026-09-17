#include "../random.HC"


auto x = 0x0123456789abcdef;

for (U64 i=0; i<64; i++){
    auto xl = RotL(x, i);
    auto xr = RotR(x, i);
    "%2d|\t%016lx\t%016lx\t%016lx\n",i,x,xl,xr;
}
