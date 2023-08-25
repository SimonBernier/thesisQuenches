//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"
#include "tdvp.h"

using namespace itensor;

//magnetic field vector
std::vector<double> hvector(int, double, double, double, double, double, double);
//makes gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int, std::vector<double>, double, SiteSet, std::vector<ITensor>);
//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);
//calculate spin-spin correlator
std::tuple<double, double, double> spinspin(int,int,MPS,SiteSet);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    if(argc < 2){ 
        printfln("Usage: %s input_file",argv[0]); 
        return 0; 
        }
    auto input = InputGroup(argv[1],"input");

    auto N = input.getInt("N");
    auto T0= input.getReal("T0", 2.);
    auto h = input.getReal("h", 3);
    auto tau = input.getReal("tau", 0.4);
    auto truncE = input.getReal("truncE", 1E-10);
    auto maxB = input.getInt("maxDim", 512);
    auto tanhshift = input.getReal("tanhshift", 2.0);
    auto dt = input.getReal("dt", 0.1);
    auto v = input.getReal("v", 1.5707963267948966); // pi/2
    
    printfln("N = %d, T0 = %0.1f, h = %0.1f, tau = %0.2f, cutoff = %0.1e, max bond dim = %d", 
                                                            N, T0, h, tau, truncE, maxB);

    // We will write into a file with the time-evolved energy density at all times.
    char schar1[128];
    int n1 = std::sprintf(schar1,"N_%d_T0_%0.1f_h_%0.1f_tau_%0.2f_cutoff_%0.0e_maxDim_%d_heisHyper.dat"
                                ,N,T0,h,tau,truncE,maxB);
  
    std::string s1(schar1);
    std::ofstream datafile;
    datafile.open(s1); // opens the file
    if( !datafile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    datafile << "time" << " " << "en0" << " " << "en-en0" << " " << "SvN" << " " << "localEn0" << " " 
            << "localEn-localEn0" << " " << "szsz" << " " << "spsm" << " " << "sz" << " " << std::endl;
    
    auto sites = SpinHalf(N);
    auto state = InitState(sites);
    for(int i = 1; i <= N; i++){
        if(i%2 == 0)
            state.set(i,"Dn");
        else
            state.set(i,"Up");
    }
    PrintData(totalQN(MPS(state)));
    
    // Create the Target Hamiltonian and find the Ground State Energy Density
    auto ampo = AutoMPO(sites);
    
    for (int b = 1; b < N; b++){
        ampo += 0.5,"S+", b, "S-", b+1;
        ampo += 0.5,"S-", b, "S+", b+1;
        ampo += 1.0,"Sz", b, "Sz", b+1;
    }
    
    //sweeps
    auto sweeps = Sweeps(5); //number of sweeps is 5
    sweeps.maxdim() = 20,50,100,200,maxB; //gradually increase states kept
    sweeps.cutoff() = truncE; //desired truncation error

    // Create the Local Energy Density Tensors
    std::vector<ITensor> PM(N-1), MP(N-1), ZZ(N-1);
    for (int b = 1; b < N; b++){
        PM[b-1] = 0.5*sites.op("S+",b)*sites.op("S-",b+1);
        MP[b-1] = 0.5*sites.op("S-",b)*sites.op("S+",b+1);
        ZZ[b-1] = 1.0*sites.op("Sz",b)*sites.op("Sz",b+1);
    }

    // Create the SzSz and S+S- correlation vector
    std::vector<double> szszcorr(N), spsmcorr(N), expSz(N);

    //magnetic field vector
    std::vector<double> hvals = hvector(N, 0.0, h, v, T0, tau, tanhshift);
    for(int b=1; b<=N; b++){
        ampo += hvals[b-1],"Sz",b;
    }
    auto H = toMPO(ampo);
    
    // Find Initial Ground State
    auto [en0,psi0] = dmrg(H,MPS(state),sweeps,{"Silent=",true});
    auto [en,psi] = dmrg(H,psi0,sweeps,{"Silent=",true});

    //calculate entanglement
    auto SvN = vonNeumannS(psi, N/2);
    //calculate local energy <psi|H(x)|psi>
    std::vector<double> localEn0(N-1), localEn(N-1);
    for (int b = 1; b < N; b++){
        psi0.position(b);
        auto ket = psi0(b)*psi0(b+1);
        localEn0[b-1] = elt( dag(prime(ket,"Site")) * PM[b-1] * ket);
        localEn0[b-1] += elt( dag(prime(ket,"Site")) * MP[b-1] * ket);
        localEn0[b-1] += elt( dag(prime(ket,"Site")) * ZZ[b-1] * ket);
        if(b==1){
            localEn0[b-1] += 1.0*hvals[b-1]*elt( dag(prime(psi0(b),"Site")) * sites.op("Sz",b) * psi0(b));
        }
        else{
            localEn0[b-1] += 0.5*hvals[b-1]*elt( dag(prime(psi0(b),"Site")) * sites.op("Sz",b) * psi0(b));
        }
        psi0.position(b+1);
        if(b == N-1){
            localEn0[b-1] += 1.0*hvals[b-1]*elt( dag(prime(psi0(b+1),"Site")) * sites.op("Sz",b+1) * psi0(b+1));
        }
        else{
            localEn0[b-1] += 0.5*hvals[b-1]*elt( dag(prime(psi0(b+1),"Site")) * sites.op("Sz",b+1) * psi0(b+1));
        }
        localEn[b-1] = localEn0[b-1]; //at t=0;
    }
    //calculate spin-spin correlation
    for (int b = 1; b <= N; b++){
        auto [szsz,spsm,sz] = spinspin(N/2+1,b,psi,sites);
        szszcorr[b-1] = szsz;
        spsmcorr[b-1] = spsm;
        expSz[b-1] = sz;
    }

    //store variables to file
    datafile << 0.0 << " " << en0 << " " << en-en0 << " " << SvN << " ";
    for (int j = 0; j < N-1; j++){
        datafile << localEn0[j] << " ";
    }
    for (int j = 0; j < N-1; j++){
        datafile << localEn[j]-localEn0[j] << " ";
    }
    //store variables to spin spin correlation file
    for (int j = 0; j < N; j++){
        datafile << szszcorr[j] << " ";
    }
    for (int j = 0; j < N; j++){
        datafile << spsmcorr[j] << " ";
    }
    for (int j = 0; j < N; j++){
        datafile << expSz[j] << " ";
    }
    datafile << std::endl;

    // time evolution parameters. Get time accuracy of 1E-4
    Real delta1 =  0.414490771794376*dt;
    Real delta2 = -0.657963087177503*dt;
    Real tval = 0.0;
    // 0.1*N/c + sqrt( (N/2/c)^2 + T0^2 ) + 2*tau*shift
    double finalTime = 0.5*double(N) + sqrt( pow(0.5*double(N)/v,2.0) + T0*T0) - T0 + 2.0*tau*tanhshift; 
    int nt = int(finalTime/dt);

    // instantaneous GS search parameters
    sweeps = Sweeps(5); //number of sweeps is 5
    sweeps.maxdim() = maxB;
    sweeps.cutoff() = truncE;
    // 4th order TDVP parameters
    auto sweeps1 = Sweeps(2); //two forward time steps of delta1
    sweeps1.maxdim() = maxB;
    sweeps1.cutoff() = truncE;
    sweeps1.niter() = 10;
    auto sweeps2 = Sweeps(1); //one backward time step of delta2
    sweeps2.maxdim() = maxB;
    sweeps2.cutoff() = truncE;
    sweeps2.niter() = 10;
    
    printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, en, SvN, maxLinkDim(psi));

    ////////////////////
    // TIME EVOLUTION //
    ////////////////////
    for (int n = 1; n <= nt; ++n){
        tval += dt;

        //update magnetic field vector
        hvals = hvector(N, tval, h, v, T0, tau, tanhshift);
        
        ampo = AutoMPO(sites);
        for (int b = 1; b < N; b++){
            ampo += 0.5,"S+", b, "S-", b+1;
            ampo += 0.5,"S-", b, "S+", b+1;
            ampo += 1.0,"Sz", b, "Sz", b+1;
        }
        for(int b=1; b<=N; b++){
            ampo += hvals[b-1],"Sz",b;
        }
        H = toMPO(ampo);

	    // calculate instantaneous ground state
        // use psi as initial condition
        en0 = dmrg(psi0,H,sweeps,{"Silent=",true});

	    // time evolve
        tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"NumCenter",2});
        tdvp(psi, H, -Cplx_i*delta2, sweeps2, {"Silent",true,"NumCenter",2});
        en = tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"NumCenter",2});
        
        //calculate entanglement entropy
        SvN = vonNeumannS(psi, N/2);
        //calculate local energy <psi|Hf(x)|psi>
        for (int b = 1; b < N; b++){
            psi0.position(b);
            psi.position(b);
            auto ket0 = psi0(b)*psi0(b+1);
            auto ket = psi(b)*psi(b+1);
            localEn0[b-1] =  eltC( dag(prime(ket0,"Site")) * PM[b-1] * ket0).real();
            localEn0[b-1] += eltC( dag(prime(ket0,"Site")) * MP[b-1] * ket0).real();
            localEn0[b-1] += eltC( dag(prime(ket0,"Site")) * ZZ[b-1] * ket0).real();
            localEn[b-1] =  eltC( dag(prime(ket,"Site")) * PM[b-1] * ket).real();
            localEn[b-1] += eltC( dag(prime(ket,"Site")) * MP[b-1] * ket).real();
            localEn[b-1] += eltC( dag(prime(ket,"Site")) * ZZ[b-1] * ket).real();

            if(b==1){
                localEn0[b-1] += hvals[b-1]*eltC( dag(prime(psi0(b),"Site")) * sites.op("Sz",b) * psi0(b)).real();
                localEn[b-1] += hvals[b-1]*eltC( dag(prime(psi(b),"Site")) * sites.op("Sz",b) * psi(b)).real();
            }
            else{
                localEn0[b-1] += 0.5*hvals[b-1]*eltC( dag(prime(psi0(b),"Site")) * sites.op("Sz",b) * psi0(b)).real();
                localEn[b-1] += 0.5*hvals[b-1]*eltC( dag(prime(psi(b),"Site")) * sites.op("Sz",b) * psi(b)).real();
            }
            psi0.position(b+1);
            psi.position(b+1);
            if(b == N-1){
                localEn0[b-1] += 1.0*hvals[b-1]*eltC( dag(prime(psi0(b+1),"Site")) * sites.op("Sz",b+1) * psi0(b+1)).real();
                localEn[b-1] += 1.0*hvals[b-1]*eltC( dag(prime(psi(b+1),"Site")) * sites.op("Sz",b+1) * psi(b+1)).real();
            }
            else{
                localEn0[b-1] += 0.5*hvals[b-1]*eltC( dag(prime(psi0(b+1),"Site")) * sites.op("Sz",b+1) * psi0(b+1)).real();
                localEn[b-1] += 0.5*hvals[b-1]*eltC( dag(prime(psi(b+1),"Site")) * sites.op("Sz",b+1) * psi(b+1)).real();
            }
        }
        //calculate spin-spin correlation
        for (int b = 1; b <= N; b++){
            auto [szsz,spsm,sz] = spinspin(N/2+1,b,psi,sites);
            szszcorr[b-1] = szsz;
            spsmcorr[b-1] = spsm;
            expSz[b-1] = sz;
        }

        //store variables to file
        datafile << tval << " " << en0 << " " << en-en0 << " " << SvN << " ";
        for (int j = 0; j < N-1; j++){
            datafile << localEn0[j] << " ";
        }
        for (int j = 0; j < N-1; j++){
            datafile << localEn[j]-localEn0[j] << " ";
        }
        //store variables to spin spin correlation file
        for (int j = 0; j < N; j++){
            datafile << szszcorr[j] << " ";
        }
        for (int j = 0; j < N; j++){
            datafile << spsmcorr[j] << " ";
        }
        for (int j = 0; j < N; j++){
            datafile << expSz[j] << " ";
        }
        datafile << std::endl;

        printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, en, SvN, maxLinkDim(psi));
    }// for n
    
    datafile.close();

    print(" END PROGRAM. TIME TAKEN :");
    printfln("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}

std::vector<double> hvector(int N, double tval, double h, double v, double T0, double tau, double tanhshift)
    {
    std::vector<double> hvals(N);
    for (int b = 1; b <= N; b++){
        double f = sqrt( pow( (double(b-N/2)-0.5)/v , 2.0)  + T0*T0) - tval - T0;
        if (b%2 == 0){
            hvals[b-1] = +h*(0.5 + 0.5*tanh( f/tau + tanhshift));
        }
        else{
            hvals[b-1] = -h*(0.5 + 0.5*tanh( f/tau + tanhshift));
        }
    }
    return hvals;
}

//calculate entanglement
Real vonNeumannS(MPS psi, int b){
    Real SvN = 0.;

    //choose orthogonality center and perform svd
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
    
    double corrZ, corrPM, expZ;

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
        auto SpSm = sites.op("S+",center+1)*sites.op("S-",center);
        corrZ = eltC( dag(prime(ket,"Site")) * SzSz * ket).real();
        corrPM = eltC( dag(prime(ket,"Site")) * SpSm * ket).real();
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
        auto SpSm = sites.op("S+",center-1)*sites.op("S-",center);
        corrZ = eltC( dag(prime(ket,"Site")) * SzSz * ket).real();
        corrPM = eltC( dag(prime(ket,"Site")) * SpSm * ket).real();
    }
    else{
        corrZ = 0.25;
        auto ket = psi(center);
        auto SpSm = 0.5*sites.op("Id",center) + sites.op("Sz",center);
        corrPM = eltC( dag(prime(ket,"Site")) * SpSm * ket).real();
    }

    return {corrZ, corrPM, expZ};

}//Szsz
