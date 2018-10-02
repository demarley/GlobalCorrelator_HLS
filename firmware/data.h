#ifndef SIMPLE_PFLOW_DATA_H
#define SIMPLE_PFLOW_DATA_H

#include "ap_int.h"
#include "ap_fixed.h"

#define DEBUG 0
#define NSECTORS 27
#define CALC_PT_IN_PROPAGATOR true

// size of the LUTs
#define PT_TABLE_SIZE 16384  // 14 unsigned bits
#define RINV_TABLE_SIZE 16384 // 32768 // 16384  // 14 unsigned bits
#define ETA_TABLE_SIZE 8192  // 13 unsigned bits
#define Z0_TABLE_SIZE 1024   // 10 unsigned bits

// range for LUTs
#define PT_RANGE 175
#define INV_PT_RANGE 1/PT_RANGE
#define RINV_RANGE 57E-4  // (0.3*3.8*0.01)/2
#define INV_RINV_RANGE 175
#define SINHETA_RANGE 6
#define INV_SINHETA_RANGE 1/SINHETA_RANGE
#define ETA_RANGE 3
#define INV_ETA_RANGE 1/ETA_RANGE
#define COSH_RANGE 3
#define INV_COSH_RANGE 1/COSH_RANGE
#define Z0_RANGE 15
#define INV_Z0_RANGE 1/Z0_RANGE

// Conversions between binary and floating point (using example file to derive)
#define RINV_CONVERSION 792057              //1314229             // 1/0.000000760902077
#define INVRINV_CONVERSION 1262535E-12   //0.000001262538462  //original: 760902077E-15   
#define PT_CONVERSION 877192982456E-10  // 1/(0.01*0.3*3.8); 87719298E-6
#define INVPT_CONVERSION 114E-4
#define ETA_CONVERSION 512
#define INVETA_CONVERSION 19531261E-10
#define PHI_CONVERSION 211216
#define INVPHI_CONVERSION 47345E-10
#define Z_CONVERSION 17
#define INVZ_CONVERSION 5859375E-8

// Muon conversions
#define MUONPT_CONVERSION 0.5
#define MUONETA_CONVERSION 0.010875
#define MUONPHI_CONVERSION 0.010908

#define M_PI_144 0.0218              // phi propagation constant

const int N_BITS_TRACK_INVPT = 15;
const int N_BITS_TRACK_PT = 14;
const int N_BITS_TRACK_ETA = 14;
const int N_BITS_TRACK_PHI = 19;
const int N_BITS_TRACK_PHIGLOBAL = 23;
const int N_BITS_TRACK_CHI2 = 10;
const int N_BITS_TRACK_Z0 = 11;
const int N_BITS_TRACK_SECTOR = 5;

typedef ap_int<N_BITS_TRACK_INVPT> invpt_t;  // inverse pt [1% at 100 GeV]
typedef ap_uint<N_BITS_TRACK_PT> pt_t;       // convert from RINV
typedef ap_int<N_BITS_TRACK_ETA> eta_t;      // eta [sinh(eta) measure to 0.005]
typedef ap_int<N_BITS_TRACK_PHI> phi_t;      // phi (50 micro-rad)
typedef ap_int<N_BITS_TRACK_PHIGLOBAL> phiglobal_t;    // phi (50 micro-rad)
typedef ap_int<N_BITS_TRACK_CHI2> chisq_t;             // chi^2 (0 - 100; 0.1 steps)
typedef ap_uint<1> q_t;                                // charge
typedef ap_int<N_BITS_TRACK_Z0> z0_t;                  // z0  (1 mm over +/-14.9 cm)
typedef ap_int<3> bx_t;
typedef ap_uint<N_BITS_TRACK_SECTOR> sector_t;

typedef ap_int<10>  etaphi_t;
typedef ap_int<5>  vtx_t;
typedef ap_uint<3>  particleid_t;

typedef ap_fixed<24,2> finvpt_t;  // inverse pt [1% at 100 GeV]
typedef ap_fixed<14,9> fpt_t;     // 1/Rinv
typedef ap_fixed<14,4> feta_t;    // eta [sinh(eta) measure to 0.005]
typedef ap_fixed<19,3> fphi_t;    // phi (50 micro-rad)
typedef ap_fixed<23,3> fphiglobal_t;    // global phi (50 micro-rad)
typedef ap_fixed<10,7> fchisq_t;  // chi^2 (0 - 100; 0.1 steps) 
typedef ap_fixed<11,5> fz0_t;     // z0  (1 mm over +/-14.9 cm) 
typedef ap_fixed<10,4> feta_m;    // eta [sinh(eta) measure to 0.005]
typedef ap_fixed<9,3> fphi_m;     // phi (50 micro-rad)
// muon data
typedef ap_uint<9> pt_m;
typedef ap_int<9> eta_m; // muon eta goes from -2.4 to 2.4 
typedef ap_uint<10> phi_m; // muon phi goes from 0 to 2pi
typedef ap_uint<4> quality_m;
// propetaphi
typedef std::pair<eta_t, phiglobal_t> etaphiglobal_t;
typedef std::pair<eta_m, phi_m> etaphiglobal_m;

		
enum PID { PID_Charged=0, PID_Neutral=1, PID_Photon=2, PID_Electron=3, PID_Muon=4 };

// PF
#define NTRACK 11
#define NMU 5
#define NCALO 0
#define NEMCALO 0
#define NPHOTON NEMCALO
#define NSELCALO 0


static fphi_t phiOffSetValues[27] = {
   -0.0387851, // 1
   0.193925, // 2
   0.426636, // 3
   0.659347, // 4
   0.892057, // 5
   1.124768, // 6
   1.357478, // 7
   1.590189, // 8
   1.822899, // 9
   2.055610, // 10
   2.288321, // 11
   2.521031, // 12
   2.753742, // 13
   2.986452, // 14
   -3.064022, // 15
   -2.831312, // 16
   -2.598601, // 17
   -2.365891, // 18
   -2.133180, // 19
   -1.900470, // 20
   -1.667759, // 21
   -1.435049, // 22
   -1.202338, // 23
   -0.969627, // 24
   -0.736917, // 25
   -0.504206, // 26
   -0.271496  // 27
};



struct CaloObj {
	pt_t hwPt;
	etaphi_t hwEta, hwPhi; // relative to the region center, at calo
};
struct HadCaloObj : public CaloObj {
	pt_t hwEmPt;
   	bool hwIsEM;
};
inline void clear(HadCaloObj & c) {
    c.hwPt = 0; c.hwEta = 0; c.hwPhi = 0; c.hwEmPt = 0; c.hwIsEM = 0; 
}

struct EmCaloObj {
	pt_t hwPt, hwPtErr;
	etaphi_t hwEta, hwPhi; // relative to the region center, at calo
};
inline void clear(EmCaloObj & c) {
    c.hwPt = 0; c.hwPtErr = 0; c.hwEta = 0; c.hwPhi = 0; 
}
/*
struct TkObj {
	pt_t hwPt, hwPtErr;
	etaphi_t hwEta, hwPhi; // relative to the region center, at calo
	z0_t hwZ0;
	bool hwTightQuality;
};
inline void clear(TkObj & c) {
    c.hwPt = 0; c.hwPtErr = 0; c.hwEta = 0; c.hwPhi = 0; c.hwZ0 = 0; c.hwTightQuality = 0;
}

struct MuObj {
	pt_t hwPt, hwPtErr;
	etaphi_t hwEta, hwPhi; // relative to the region center, at vtx(?)
};
inline void clear(MuObj & c) {
    c.hwPt = 0; c.hwPtErr = 0; c.hwEta = 0; c.hwPhi = 0; 
}
struct PFChargedObj {
	pt_t hwPt;
	etaphi_t hwEta, hwPhi; // relative to the region center, at calo
	particleid_t hwId;
	z0_t hwZ0;
};
*/
struct PFNeutralObj {
	pt_t hwPt;
	etaphi_t hwEta, hwPhi; // relative to the region center, at calo
	particleid_t hwId;
  pt_t hwPtPuppi;
};
struct VtxObj {
	pt_t  hwSumPt;
	z0_t  hwZ0;
	vtx_t mult;
	particleid_t hwId;
};




struct TkObj {
    // Tracks
    invpt_t hwRinv;
    pt_t hwPt;
    pt_t hwPtErr;
    eta_t hwSinhEta;
    eta_t hwEta;
    phi_t hwPhi;
    phiglobal_t hwPhiGlobal;
    z0_t hwZ0;
    sector_t hwSector;
    q_t hwQ;
    chisq_t hwX2;
    q_t hwValid;
    bx_t hwBX;
    eta_t hwPropEta;
    phiglobal_t hwPropPhi;

    // constructor
    TkObj() : 
      hwRinv(0),
      hwPt(0),
      hwSinhEta(0),  
      hwEta(0),
      hwPhi(0),
      hwPhiGlobal(0),
      hwZ0(0),
      hwSector(0),
      hwQ(0),
      hwX2(0),
      hwValid(0),
      hwBX(0){}
};
inline void clear(TkObj & c) {
      c.hwRinv = 0;
      c.hwPt = 0;
      c.hwSinhEta = 0;
      c.hwEta = 0;
      c.hwPhi = 0;
      c.hwPhiGlobal = 0;
      c.hwZ0 = 0;
      c.hwSector = 0;
      c.hwQ = 0;
      c.hwX2 = 0;
      c.hwValid = 0;
      c.hwBX = 0;
}


struct MuObj {
    // Muons
    pt_m hwPt;
    pt_m hwPtErr;
    eta_m hwEta;
    phi_m hwPhi;
    q_t hwQ;
    q_t hwValid;
    bx_t hwBX;
    z0_t hwZ0;
    quality_m hwQuality;

    // constructor
    MuObj() : 
      hwPt(0),
      hwEta(0),
      hwPhi(0),
      hwQ(0),
      hwValid(0),
      hwBX(0),
      hwZ0(0),
      hwQuality(0){}
};
inline void clear(MuObj & c) {
      c.hwPt  = 0;
      c.hwEta = 0;
      c.hwPhi = 0;
      c.hwZ0  = 0;
      c.hwQ   = 0;
      c.hwBX  = 0;
      c.hwValid = 0;
      c.hwQuality = 0;
}


// need a 64-bit trackmuon at least
struct PFChargedObj {
    pt_t hwPt;
    pt_t hwPtErr;
    eta_m hwEta;
    phi_m hwPhi;
    q_t hwQ;
    z0_t hwZ0;
    q_t hwValid;
    bx_t hwBX;
    quality_m hwQuality;

    // constructor
    PFChargedObj() : 
      hwPt(0),
      hwEta(0),
      hwPhi(0),
      hwQ(0),
      hwValid(0),
      hwZ0(0),
      hwBX(0),
      hwQuality(0){}
};
inline void clear(PFChargedObj & c) {
      c.hwPt  = 0;
      c.hwEta = 0;
      c.hwPhi = 0;
      c.hwZ0  = 0;
      c.hwQ   = 0;
      c.hwBX  = 0;
      c.hwQuality = 0;
      c.hwValid = 0;
}







#define MP7_NCHANN 144
#define CTP7_NCHANN_IN 67
#define CTP7_NCHANN_OUT 48
typedef ap_uint<32> MP7DataWord;

#endif
