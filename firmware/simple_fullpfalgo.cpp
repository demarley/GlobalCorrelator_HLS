#include "simple_fullpfalgo.h"
#include "mp7pf_encoding.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif



feta_t deltaR(MuObj mu, TkObj tk){
    /* Calculate eta and phi and deltaR match */
    feta_t tkEta       = tk.hwPropEta*INVETA_CONVERSION;
    fphiglobal_t tkPhi = normalizePhi(tk.hwPropPhi*INVPHI_CONVERSION);

    feta_m muEta = (1- 2*std::bitset<9>(mu.hwEta)[8]) * from_twos_complement<9>(mu.hwEta) * MUONETA_CONVERSION;
    fphi_m muPhi = normalizePhi(mu.hwPhi * MUONPHI_CONVERSION);

    feta_t dr = deltaR_hw( muEta,muPhi,tkEta,tkPhi );
    return dr;
}


int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}
template<int NB>
ap_uint<NB> dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<NB> max) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    int dr2 = deta*deta + dphi*dphi;
    return (dr2 < int(max) ? ap_uint<NB>(dr2) : max);
}
template<int NB, typename PTS_t>
ap_uint<NB> dr2_dpt_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, pt_t pt1, pt_t pt2, PTS_t ptscale, ap_uint<NB> dr2max, ap_uint<NB> max) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    int dr2 = deta*deta + dphi*dphi;
    pt_t dpt = pt1 - pt2;
    if (dpt < 0) dpt = 0;
    ap_int<26> dpt2 = (dpt > 5792) ? ap_int<26>((1<<25)-1) : ap_int<26>(dpt*dpt);
    int dq = dr2 + (dpt2*ptscale >> 8);
    return ((dr2 < int(dr2max)) && (dq < int(max))) ? ap_uint<NB>(dq) : max;
}

template<typename T, int NIn, int NOut>
void ptsort_hwopt(T in[NIn], T out[NOut]) {
    T tmp[NOut];
    #pragma HLS ARRAY_PARTITION variable=tmp complete

    for (int iout = 0; iout < NOut; ++iout) {
        #pragma HLS unroll
        tmp[iout].hwPt = 0;
    }

    for (int it = 0; it < NIn; ++it) {
        for (int iout = NOut-1; iout >= 0; --iout) {
            if (tmp[iout].hwPt <= in[it].hwPt) {
                if (iout == 0 || tmp[iout-1].hwPt > in[it].hwPt) {
                    tmp[iout] = in[it];
                } else {
                    tmp[iout] = tmp[iout-1];
                }
            }
        }
    }
    for (int iout = 0; iout < NOut; ++iout) {
        out[iout] = tmp[iout];
    }

}

template<int DR2MAX>
void init_dr2max_times_pterr2_inv(int vals[512]) {
    for (int i = 0; i < 512; ++i) {
    	int tmp = (DR2MAX<<8)/(i?i*i:1), int18_max = (1<<17)-1;
        vals[i] = (tmp > int18_max ? int18_max : tmp);
    }
}


//-------------------------------------------------------
// TK-MU Algos
//-------------------------------------------------------
void spfph_trkprop(TkObj track[NTRACK]){
    /* Propagate tracks to 2nd muon station */
    #pragma HLS ARRAY_PARTITION variable=track complete
    for (int tk=0; tk<NTRACK; ++tk){
        etaphiglobal_t etaphi;
        prop_hw(track[tk], etaphi);
        track[tk].hwPropEta = etaphi.first;
        track[tk].hwPropPhi = etaphi.second;
    }
    return;
}

void spfph_mu2trk_drvals(MuObj mu[NMU], TkObj track[NTRACK], feta_t mu_track_drval[NMU][NTRACK]) {
    const ap_uint<12> DR2MAX = PFALGO3_DR2MAX_TK_MU;
    for (int im = 0; im < NMU; ++im) {
        for (int it = 0; it < NTRACK; ++it) {
            feta_t dr = deltaR(mu[im],track[it]);
            if (dr<DR2MAX)
                mu_track_drval[im][it] = (dr > 0 ? dr : feta_t(-dr));
            else
                mu_track_drval[im][it] = DR2MAX;
        } // end loop over tracks
    } // end loop over muons

    return;
}

void spfph_mu2trk_linkstep(MuObj mu[NMU], feta_t mu_track_drval[NMU][NTRACK], ap_uint<NMU> mu_track_link_bit[NTRACK]) {
    /* Link tracks to muons */
    for (int im = 0; im < NMU; ++im) {
        for (int it = 0; it < NTRACK; ++it) {
            feta_t mydr = mu_track_drval[im][it];
            bool link = (mydr < 0.2);

            // confirm that this track has the smallest DeltaR value to link with muon
            for (int j=0; j < NTRACK; ++j) {
                if (it <= j) link = link && (mu_track_drval[im][j] >= mydr);
                else         link = link && (mu_track_drval[im][j] >  mydr);
            }   
            mu_track_link_bit[it][im] = link;
        } // end loop over tracks
    } // end loop over muons

    return;
}

void spfph_mutrk_link(MuObj mu[NMU], TkObj track[NTRACK], ap_uint<NMU> mu_track_link_bit[NTRACK]) {
    /* Function to calculate Delta(pT) between tracks and muons then link them */    
    #pragma HLS ARRAY_PARTITION variable=mu complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu_track_link_bit complete dim=0

    feta_t drvals[NMU][NTRACK];
    #pragma HLS ARRAY_PARTITION variable=drvals complete dim=0

    etaphiglobal_t etaphi;

    spfph_trkprop(track);
    spfph_mu2trk_drvals(mu, track, drvals);
    spfph_mu2trk_linkstep(mu, drvals, mu_track_link_bit);

    return;
}

void spfph_mualgo(MuObj mu[NMU], TkObj track[NTRACK], ap_uint<NMU> mu_track_link_bit[NTRACK], PFChargedObj pfmuout[NMU], bool isMu[NTRACK]) {
    #pragma HLS ARRAY_PARTITION variable=isMu complete

    for (int im = 0; im < NMU; ++im) {
        bool good = false;
        int ibest = -1;
        for (int it = 0; it < NTRACK; ++it) {
            if (mu_track_link_bit[it][im]){ good = true; ibest = it; }
        }
        if (mu[im].hwPt > 0 && good && ibest != -1) {
            pfmuout[im].hwPt  = track[ibest].hwPt;
            pfmuout[im].hwEta = track[ibest].hwEta;
            pfmuout[im].hwPhi = track[ibest].hwPhi;
            pfmuout[im].hwZ0  = track[ibest].hwZ0;
            isMu[ibest] = 1;
        } else {
            pfmuout[im].hwPt  = 0;
            pfmuout[im].hwEta = 0;
            pfmuout[im].hwPhi = 0;
            pfmuout[im].hwZ0  = 0;
        }
    }
}

//-------------------------------------------------------
// PF Algos
//-------------------------------------------------------

void pfalgo3_full(TkObj track[NTRACK], MuObj mu[NMU], PFChargedObj outmu[NMU]) {
    /* Build PF objects */
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete    
    #pragma HLS ARRAY_PARTITION variable=outmu complete

    #pragma HLS pipeline II=HLS_pipeline_II

    // ---------------------------------------------------------------
    // TK-MU Linking
    ap_uint<NMU> mu_track_link_bit[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=mu_track_link_bit complete

    bool isMu[NTRACK];
    for (int it=0; it<NTRACK; ++it) { isMu[it]=0; }  // initialize matches to 0

    spfph_mutrk_link(mu, track, mu_track_link_bit);
    spfph_mualgo(mu, track, mu_track_link_bit, outmu, isMu);

    return;
}


void mp7wrapped_pack_in(TkObj track[NTRACK], MuObj mu[NMU], MP7DataWord data[MP7_NCHANN]) {
    /* MP7 pack inputs */
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete

    // pack inputs
    assert(2*NEMCALO + 2*NTRACK + 2*NCALO + 2*NMU <= MP7_NCHANN);
    #define TKOFFS 0
    #define MUOFFS 3*NTRACK+TKOFFS
    mp7_pack<NTRACK,TKOFFS>(track, data);
    mp7_pack<NMU,MUOFFS>(mu, data);

    return;
}

void mp7wrapped_unpack_in(MP7DataWord data[MP7_NCHANN], TkObj track[NTRACK], MuObj mu[NMU]) {
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete

    // unpack inputs
    assert(2*NEMCALO + 2*NTRACK + 2*NCALO + 2*NMU <= MP7_NCHANN);
    #define TKOFFS 0
    #define MUOFFS 3*NTRACK+TKOFFS

    mp7_unpack<NTRACK,TKOFFS>(data, track);
    mp7_unpack<NMU,MUOFFS>(data, mu);
}

void mp7wrapped_pack_out( PFChargedObj pfmu[NMU], MP7DataWord data[MP7_NCHANN]) {
    /* MP7 pack outputs */
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete

    // pack outputs
    assert(2*NTRACK + 2*NPHOTON + 2*NSELCALO + 2*NMU <= MP7_NCHANN);
    #define PFMUOFFS 0

    for (unsigned int i = 0; i < NMU; ++i) {
        data[2*i+0+PFMUOFFS] = ( pfmu[i].hwPt, pfmu[i].hwEta, pfmu[i].hwPhi );
    }

    return;
}

void mp7wrapped_unpack_out( MP7DataWord data[MP7_NCHANN], PFChargedObj pfmu[NMU]) {
    /* MP7 unpack outputs */
    #pragma HLS ARRAY_PARTITION variable=data complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete

    // unpack outputs
    assert(2*NTRACK + 2*NPHOTON + 2*NSELCALO + 2*NMU <= MP7_NCHANN);
    #define PFMUOFFS 0

    for (unsigned int i = 0; i < NMU; ++i) {
        pfmu[i].hwPt  = data[2*i+0+PFMUOFFS](28,19);
        pfmu[i].hwEta = data[2*i+1+PFMUOFFS](19,10);
        pfmu[i].hwPhi = data[2*i+1+PFMUOFFS](10,0);
    }

    return;
}


void mp7wrapped_pfalgo3_full(MP7DataWord input[MP7_NCHANN], MP7DataWord output[MP7_NCHANN]) {
    /* Top-level function */
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=output complete
    #pragma HLS INTERFACE ap_none port=output

    #pragma HLS pipeline II=HLS_pipeline_II

    TkObj track[NTRACK];
    MuObj mu[NMU];
    #pragma HLS ARRAY_PARTITION variable=track complete
    #pragma HLS ARRAY_PARTITION variable=mu complete

    PFChargedObj pfmu[NMU];
    #pragma HLS ARRAY_PARTITION variable=pfmu complete

    mp7wrapped_unpack_in(input, track, mu);
    pfalgo3_full(track, mu, pfmu);
    mp7wrapped_pack_out(pfmu, output);

    return;
}

// THE END //
