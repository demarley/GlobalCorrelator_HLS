// Encoding MP7

template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack(TkObj track[N], MP7DataWord data[]) {
    #pragma HLS inline
    ap_int<13> empty_phi = 0;
    ap_int<3>  empty_rinv_eta = 0;
    ap_int<16> empty_z0_sector = 0;
    for (unsigned int i = 0; i < N; ++i) {
        data[3*i+0+OFFS] = ( empty_phi,track[i].hwPhi );
        data[3*i+1+OFFS] = ( empty_rinv_eta,track[i].hwRinv, track[i].hwSinhEta);
        data[3*i+2+OFFS] = ( empty_z0_sector,track[i].hwZ0, track[i].hwSector );
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack(MuObj mu[N], MP7DataWord data[]) {
    #pragma HLS inline
    ap_int<4> empty = 0;
    for (unsigned int i = 0; i < N; ++i) {
        data[i+OFFS] = ( empty,mu[i].hwPt, mu[i].hwEta, mu[i].hwPhi );
    }
}


template<unsigned int N, unsigned int OFFS> 
inline void mp7_unpack(MP7DataWord data[], TkObj track[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        track[i].hwPhi     = data[3*i+0+OFFS](18, 0);
        track[i].hwRinv    = data[3*i+1+OFFS](29,14);
        track[i].hwSinhEta = data[3*i+1+OFFS](14, 0);
        track[i].hwZ0      = data[3*i+2+OFFS](16, 5);
        track[i].hwSector  = data[3*i+2+OFFS]( 5, 0);
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void mp7_unpack(MP7DataWord data[], MuObj mu[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        mu[i].hwPt    = data[i+OFFS](28,19);
        mu[i].hwEta   = data[i+OFFS](19,10);
        mu[i].hwPhi   = data[i+OFFS](10, 0);
    }
}


// THE END
