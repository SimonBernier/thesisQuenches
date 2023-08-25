//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"

using namespace itensor;

//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);
//calculate spin-spin correlator
std::tuple<double, double, double> spinspin(int,int,MPS,SiteSet);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    int N = 16;

    if(argc > 1)
        N = std::stoi(argv[1]);

    // We will write into a file with the time-evolved energy density at all times.
    char schar1[64], schar2[64];
    int n1 = std::sprintf(schar1,"N_%d_TFI1dCritical_En.dat",N);
    int n2 = std::sprintf(schar2,"N_%d_TFI1dCritical_SSC.dat",N);
    std::string s1(schar1), s2(schar2);
    std::ofstream enerfile, sscfile;

    enerfile.open(s1); // opens the file
    if( !enerfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    enerfile << "energy" << " " << "var" << " " << "SvN" << " " << "bondDim" << " " << "localEnergy" << " " << std::endl;

    sscfile.open(s2); // opens the file
    if( !sscfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    sscfile << "szsz" << " " << "sxsx" << " " << "sz" << " " << std::endl;

    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity=",true});

    // Create the Target Hamiltonian and find the Ground State Energy Density
    auto ampo = AutoMPO(sites);

    for (int b = 1; b < N; b++){
        ampo += -4.0,"Sx", b, "Sx", b+1;
    }
    for (int b = 1; b <= N; b++){
        ampo += -2.0,"Sz",b;
    }
    auto H = toMPO(ampo);
    //PrintData(H);

    //sweeps
    auto sweeps = Sweeps(10); //number of sweeps is 10
    sweeps.maxdim() = 10,20,40,80,100; //gradually increase states kept
    sweeps.cutoff() = 1E-10; //desired truncation error

    // Create the Local Energy Density Tensors
    std::vector<double> localEnergy(N-1);
    std::vector<ITensor> LED(N-1);
    for (int b = 1; b < N; b++){
        LED[b-1] =  -4.0*sites.op("Sx",b)*sites.op("Sx",b+1);
        LED[b-1] += -2.0*sites.op("Sz",b)*sites.op("Id",b+1);
    }

    // Create the SzSz and S+S- correlation vector
    std::vector<double> szszcorr(N), sxsxcorr(N), szexp(N);

    // make initial state with the paramagnetic parity
    auto state = InitState(sites);

    for(int i = 1; i <= N; i++){
        state.set(i,"Up");
    }
    auto initState = MPS(state);
    PrintData(totalQN(initState));

    // Find Initial Ground State
    auto [energy,psi] = dmrg(H,initState,sweeps,{"Silent=",true});
    auto var = inner(H, psi, H, psi) - energy*energy;
    //calculate entanglement
    auto SvN = vonNeumannS(psi, N/2);
    //calculate local energy <psi|H(x)|psi>
    for (int b = 1; b < N; b++){
        psi.position(b);
        auto ket = psi(b)*psi(b+1);
        auto hi = LED[b-1];
        hi += 1.0*sites.op("Sz",b)*sites.op("Id",b+1);
        hi += -1.0*sites.op("Id",b)*sites.op("Sz",b+1);
        localEnergy[b-1] = elt( dag(prime(ket,"Site")) * hi * ket);
    }
    enerfile << energy << " " << var << " " << SvN << " ";
    auto bonds = linkInds(psi); //get bond dimensions
    for (int j=1; j<N; j++){
        enerfile << dim(bonds[j-1]) << " ";
    }
    for (int j = 1; j < N; j++){
        enerfile << localEnergy[j-1] << " ";
    }
    enerfile << std::endl;

    //calculate spin-spin correlation
    for (int b = 1; b <= N; b++){
        auto [szsz,sxsx,sz] = spinspin(N/2+1,b,psi,sites);
        szszcorr[b-1] = szsz;
        sxsxcorr[b-1] = sxsx;
        szexp[b-1] = sz;
    }
    //store variables to spin spin correlation file
    for (int j = 0; j < N; j++){
        sscfile << szszcorr[j] << " ";
    }
    for (int j = 0; j < N; j++){
        sscfile << sxsxcorr[j] << " ";
    }
    for (int j = 0; j < N; j++){
        sscfile << szexp[j] << " ";
    }
    sscfile << std::endl;

    printfln("N = %d, energy = %0.3f, SvN = %0.3f, maxDim = %d", N, energy, SvN, maxLinkDim(psi));

    std::cout<< std::endl << " END PROGRAM. ";

    std::printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}

Real vonNeumannS(MPS psi, int b){
    Real SvN = 0.;

    //calculate entanglement
    psi.position(b);
    auto l = leftLinkIndex(psi,b);
    auto s = siteIndex(psi,b);
    auto [U,S,V] = svd(psi(b),{l,s});
    auto u = commonIndex(U,S);

    //Apply von Neumann formula
    //to the squares of the singular values
    for(auto n : range1(dim(u))){
        auto Sn = elt(S,n,n);
        auto p = sqr(Sn);
        if(p > 1E-12) SvN += -p*log(p);
    }
    return SvN;

}//vonNeumannS

//calculate spin-spin correlator
std::tuple<double, double, double> spinspin(int center, int b, MPS psi, SiteSet sites){

    double corrZ, corrX, expZ;

    psi.position(b);
    expZ = eltC(dag(prime(psi(b),"Site")) * sites.op("Sz",b) * psi(b)).real();
    if(b>center){ //bring site b next to the center from right
        for(int n=b-1; n>center; n--){
            auto g = BondGate(sites,n,n+1);
            auto AA = psi(n)*psi(n+1)*g.gate(); //contract over bond n
            AA.replaceTags("Site,1","Site,0");
            psi.svdBond(g.i1(), AA, Fromright); //svd from the right
            psi.position(g.i1()); //move orthogonality center to the left
        }
        auto ket = psi(center)*psi(center+1);
        auto SzSz = sites.op("Sz",center)*sites.op("Sz",center+1);
        auto SpSm = sites.op("Sx",center+1)*sites.op("Sx",center);
        corrZ = eltC( dag(prime(ket,"Site")) * SzSz * ket).real();
        corrX = eltC( dag(prime(ket,"Site")) * SpSm * ket).real();
    }
    else if(b<center){ //bring site b next to the center from left
        for(int n=b; n<center-1; n++){
          auto g = BondGate(sites,n,n+1);
          auto AA = psi(n)*psi(n+1)*g.gate(); //contract over bond n
          AA.replaceTags("Site,1","Site,0");
          psi.svdBond(g.i1(), AA, Fromleft); //svd from the right
          psi.position(g.i2()); //move orthogonality center to the right 
        }
        auto ket = psi(center-1)*psi(center);
        auto SzSz = sites.op("Sz",center-1)*sites.op("Sz",center);
        auto SpSm = sites.op("Sx",center-1)*sites.op("Sx",center);
        corrZ = eltC( dag(prime(ket,"Site")) * SzSz * ket).real();
        corrX = eltC( dag(prime(ket,"Site")) * SpSm * ket).real();
    }
    else{
        corrZ = 0.25; corrX = 0.25;
    }

    return {corrZ, corrX, expZ};

}//spinspin
