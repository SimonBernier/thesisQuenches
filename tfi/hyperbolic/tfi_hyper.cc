#include "itensor/all.h"

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
    auto T0 = input.getReal("T0", 2.);
    auto h = input.getReal("h", 1.5);
    auto tau = input.getReal("tau", 2.0);
    auto truncE = input.getReal("truncE", 1E-10);
    auto maxB = input.getInt("maxDim", 512);
    auto tanhshift = input.getReal("tanhshift", 2.0);
    auto dt = input.getReal("dt", 0.1);
    auto v = input.getReal("v", 2.0);
    
    printfln("N = %d, T0 = %0.1f, h = %0.1f, tau = %0.2f, cutoff = %0.1e, max bond dim = %d", 
                                                            N, T0, h, tau, truncE, maxB);

    // We will write into a file with the time-evolved energy density at all times.
    char schar2[128], schar3[128];
    int n2 = std::sprintf(schar2,"N_%d_T0_%0.1f_h_%0.1f_tau_%0.2f_cutoff_%0.0e_maxDim_%d_TFI1dHyper_En.dat"
                                        ,N,T0,h,tau,truncE,maxB);
    int n3 = std::sprintf(schar3,"N_%d_T0_%0.1f_h_%0.1f_tau_%0.2f_cutoff_%0.0e_maxDim_%d_TFI1dHyper_SSC.dat"
                                        ,N,T0,h,tau,truncE,maxB);

    std::string s2(schar2), s3(schar3);
    std::ofstream enerfile, sscfile;
    enerfile.open(s2); // opens the file
    if( !enerfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    enerfile << "time" << " " << "energy" << " " << "SvN" << " " << "truncerr" << " " << "bondDim" << " " << "localEnergy" << " " << std::endl;
    
    sscfile.open(s3); // opens the file
    if( !sscfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    sscfile << "time" << " " << "szsz" << " " << "sxsx" << " " << "sz" << " " << std::endl;

    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity=",true});
    
    // Create the Target Hamiltonian and find the Ground State Energy Density
    auto ampo = AutoMPO(sites);
    
    for (int b = 1; b < N; b++){
        ampo += -4.0,"Sx", b, "Sx", b+1;
    }
    for (int b = 1; b <= N; b++){
        ampo += -2.0,"Sz",b;
    }
    auto Hfinal = toMPO(ampo);
    
    //sweeps
    auto sweeps = Sweeps(5); //number of sweeps is 5
    sweeps.maxdim() = 10,20,100,200,maxB; //gradually increase states kept
    sweeps.cutoff() = truncE; //desired truncation error

    // Create the Local Energy Density Tensors
    std::vector<double> localEnergy(N-1);
    std::vector<ITensor> LED(N-1);
    for (int b = 1; b < N; b++){
        LED[b-1] =  -4.0*sites.op("Sx",b)*sites.op("Sx",b+1);
        LED[b-1] += -2.0*sites.op("Sz",b)*sites.op("Id",b+1);
    }

    // Create the SzSz and S+S- correlation vector
    std::vector<double> szszcorr(N), sxsxcorr(N), szexp(N);

    //magnetic field vector
    std::vector<double> hvals = hvector(N, 0.0, h, v, T0, tau, tanhshift);
    for(int b=1; b<=N; b++){
        ampo += -2.0*hvals[b-1],"Sz",b;
    }
    
    // make initial state with the paramagnetic parity
    auto state = InitState(sites);

    for(int i = 1; i <= N; i++){
        state.set(i,"Up");
    }
    auto initState = MPS(state);
    PrintData(totalQN(initState));

    // Find Initial Ground State
    auto [energy,psi] = dmrg(toMPO(ampo),initState,sweeps,{"Silent=",true});
    energy = inner(psi, Hfinal, psi);
    //calculate entanglement
    auto SvN = vonNeumannS(psi, N/2);
    //calculate local energy <psi|Hf(x)|psi>
    for (int b = 1; b < N; b++){
        psi.position(b);
        auto ket = psi(b)*psi(b+1);
        auto hi = LED[b-1];
        hi += 1.0*sites.op("Sz",b)*sites.op("Id",b+1);
        hi += -1.0*sites.op("Id",b)*sites.op("Sz",b+1);
        localEnergy[b-1] = elt( dag(prime(ket,"Site")) * hi * ket);
    }
    enerfile << 0.0 << " " << energy << " " << SvN << " " << 0.0 << " ";
    IndexSet bonds = linkInds(psi); //get bond dimensions
    for (int j=0; j < N-1; j++){
        enerfile << dim(bonds[j]) << " ";
    }
    for (int j=0; j < N-1; j++){
        enerfile << localEnergy[j] << " ";
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
    sscfile << 0.0 << " "; 
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

    Real delta1 =  0.414490771794376*dt;
    Real delta2 = -0.657963087177503*dt;
    Real tval = 0.0;
    // 0.5*N/c + sqrt( (N/2/c)^2 + T0^2 ) + 2*tau*shift
    double finalTime = 0.5*double(N)/v + sqrt( pow(0.5*double(N)/v,2.0) + T0*T0) - T0 + 2.0*tau*tanhshift; 
    int nt = int(finalTime/dt);
    auto args = Args("Cutoff=",truncE,"MaxDim=",maxB);
    
    printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, energy, SvN, maxLinkDim(psi));

    auto tot_truncerr=0.;

    ////////////////////
    // TIME EVOLUTION //
    ////////////////////
    for (int n = 1; n <= nt; ++n){
        tval += dt;

        //update magnetic field vector
        hvals = hvector(N, tval, h, v, T0, tau, tanhshift);
        
        // TEBD time update
        std::vector<BondGate> gates;
        auto gatesdelta1 = makeGates(N, hvals, delta1, sites, LED);
        auto gatesdelta2 = makeGates(N, hvals, delta2, sites, LED);
        gates = gatesdelta1;
        gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
        gates.insert(std::end(gates), std::begin(gatesdelta2), std::end(gatesdelta2));
        gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
        gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));

        // apply Trotter gates
        auto truncerr = gateTEvol(gates,dt,dt,psi,{args,"Verbose=",false});
        tot_truncerr += truncerr;
        psi.orthogonalize(args); //orthogonalize to minimize bond dimensions

        // calculate energy <psi|Hf|psi>
        auto en = innerC(psi, Hfinal, psi).real();
        //calculate entanglement entropy
        SvN = vonNeumannS(psi, N/2);
        //calculate local energy <psi|Hf(x)|psi>
        for (int b = 1; b < N; b++){
            psi.position(b);
            auto ket = psi(b)*psi(b+1);
            auto hi = LED[b-1];
            hi += 1.0*sites.op("Sz",b)*sites.op("Id",b+1);
            hi += -1.0*sites.op("Id",b)*sites.op("Sz",b+1);
            localEnergy[b-1] = eltC( dag(prime(ket,"Site")) * hi * ket).real();
        }

        enerfile << tval << " " << en << " " << SvN << " " << tot_truncerr << " ";
        IndexSet bonds = linkInds(psi); //get bond dimensions
        for (int j = 0; j < N-1; j++){
            enerfile << dim(bonds[j]) << " ";
        }
        for (int j = 0; j < N-1; j++){
            enerfile << localEnergy[j] << " ";
        }
        enerfile << std::endl;

        if (n % 10 == 0){
            //calculate spin-spin correlation
            for (int b = 1; b <= N; b++){
                auto [szsz,sxsx,sz] = spinspin(N/2+1,b,psi,sites);
                szszcorr[b-1] = szsz;
                sxsxcorr[b-1] = sxsx;
                szexp[b-1] = sz;
            }
            //store variables to spin spin correlation file
            sscfile << tval << " "; 
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
        } //if

        printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f \n\t maxDim = %d, total truncation error = %0.2e", 
                    tval, en, SvN, maxLinkDim(psi), tot_truncerr);
    }// for n
    
    enerfile.close();

    std::cout<< std::endl << " END PROGRAM. ";
    
    std::printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}

// return vector of h values
std::vector<double> hvector(int N, double tval, double h, double v, double T0, double tau, double tanhshift)
    {
    std::vector<double> hvals(N);
    for (int b = 1; b <= N; b++){
        double f = sqrt( pow( (double(b-N/2)-0.5)/v , 2.0)  + T0*T0) - tval - T0;
        hvals[b-1] = h*(0.5 + 0.5*tanh( f/tau + tanhshift));
    }
    return hvals;
}

// second order Trotter breakup of time step dt
// returns a vector of gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int L, std::vector<double> h, double dt, SiteSet sites, std::vector<ITensor> LED)
    {
    std::vector<BondGate> gates; 
    //Create the gates exp(-i*tstep/2*hterm)
    for(int i=1; i<=L; i++){
        if(i<L){
            auto hterm = LED[i-1];
            hterm += -2.0*h[i-1]*op(sites,"Sz",i)*op(sites,"Id",i+1);
            auto g = BondGate(sites,i,i+1,BondGate::tReal,dt/2.,hterm);
            gates.push_back(g);
        }
        else{
            auto hterm = -2.0*(1.+h[i-1])*op(sites,"Id",i-1)*op(sites,"Sz",i);
            auto g = BondGate(sites,i-1,i,BondGate::tReal,dt/2.,hterm);
            gates.push_back(g);
        }
    } // for i

  //Create the gates exp(-i*tstep/2*hterm) in reverse order 
  for(int i=L; i>=1; i--){
    if(i<L){
        auto hterm = LED[i-1];
        hterm += -2.0*h[i-1]*op(sites,"Sz",i)*op(sites,"Id",i+1);
        auto g = BondGate(sites,i,i+1,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);
    }
    else{
        auto hterm = -2.0*(1.+h[i-1])*op(sites,"Id",i-1)*op(sites,"Sz",i);
        auto g = BondGate(sites,i-1,i,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);
    }
  }// for i

  return gates;
  
}// makeGates

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