#include <cstdio>
#include "firmware/simple_fullpfalgo.h"
#include "utils/random_inputs.h"
#include "utils/DiscretePFInputs_IO.h"
#include "utils/pattern_serializer.h"
#include "utils/test_utils.h"

#define NTEST 10


int main() {
    // input TP objects
    TkObj track[NTRACK]; 
    MuObj mu[NMU];

    // output PF objects
    PFChargedObj outch[NTRACK], outch_ref[NTRACK];
    PFChargedObj outmupf[NMU], outmupf_ref[NMU];

    CTP7PatternSerializer serInPatterns4( "ctp7_input_tkmu_patterns_nomux.txt",CTP7_NCHANN_IN,true); //"ctp7_input_patterns_nomux.txt",CTP7_NCHANN_IN, true);
    CTP7PatternSerializer serOutPatterns4("ctp7_output_patterns_nomux.txt",CTP7_NCHANN_OUT, false);

    HumanReadablePatternSerializer serHR("human_readable_patterns.txt");
    HumanReadablePatternSerializer debugHR("-"); // this will print on stdout, we'll use it for errors

    // -----------------------------------------
    // run multiple tests
    for (int test = 1; test <= NTEST; ++test) {

        // initialize TP objects
        for (int i = 0; i < NTRACK; ++i) {
            track[i].hwPt    = 0;
            track[i].hwPtErr = 0;
            track[i].hwEta   = 0;
            track[i].hwPhi   = 0;
            track[i].hwZ0    = 0;
        }
        for (int i = 0; i < NMU; ++i) {
            mu[i].hwPt    = 0;
            mu[i].hwPtErr = 0;
            mu[i].hwEta   = 0;
            mu[i].hwPhi   = 0;
        }

        MP7DataWord data_in[CTP7_NCHANN_IN], data_out[CTP7_NCHANN_OUT];
        // initialize
        for (int i=0; i<CTP7_NCHANN_IN; ++i)  data_in[i]  = 0;
        for (int i=0; i<CTP7_NCHANN_OUT; ++i) data_out[i] = 0;

        mp7wrapped_pack_in(track, mu, data_in);
        MP7_TOP_FUNC(data_in, data_out);
        mp7wrapped_unpack_out(data_out, outmupf);

        MP7_REF_FUNC(track, mu, outmupf_ref);

        // write out patterns for CTP7 board hardware or simulator test
        serInPatterns4(data_in,CTP7_NCHANN_IN);
        serOutPatterns4(data_out,CTP7_NCHANN_OUT);       

        if (!CTP7_VALIDATE) continue;

        // -----------------------------------------
        // validation against the reference algorithm
        int errors = 0; 
        int ntot   = 0;
        int nmu    = 0;

        for (int i = 0; i < NMU; ++i) {
            if (!pf_equals(outmupf_ref[i], outmupf[i], "PF Muon", i)) errors++;
            if (outmupf_ref[i].hwPt > 0) { ntot++; nmu++; }
        }        

        if (errors != 0) {
            printf("Error in computing test %d (%d)\n", test, errors);
            printf("Inputs: \n"); debugHR.dump_inputs(track, mu);
            printf("Reference output: \n"); debugHR.dump_outputs(outch_ref,outmupf_ref);
            printf("Current output: \n"); debugHR.dump_outputs(outch,outmupf);
            return 1;
        }
        else {
            printf("Passed test %d (%d)\n", test, nmu);
        }
    }

    return 0;
}
