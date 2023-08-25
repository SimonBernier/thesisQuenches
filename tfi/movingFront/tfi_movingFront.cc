//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"

using namespace itensor;

//magnetic field vector
std::vector<double> hvector(int, double, double, double, double, double);
//makes gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int, std::vector<double>, double, SiteSet, std::vector<ITensor>);
//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);
//calculate spin-spin correlator
double spinspin(int,int,MPS,SiteSet);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    if(argc < 2){ 
        printfln("Usage: %s input_file",argv[0]); 
        return 0; 
        }
    auto input = InputGroup(argv[1],"input");

    auto N = input.getInt("N");
    auto v = input.getReal("v", 2.);
    auto h = input.getReal("h", 1.5);
    auto tau = input.getReal("tau", 2.0);
    auto truncE = input.getReal("truncE", 1E-10);
    auto maxB = input.getInt("maxDim", 512);
    auto tanhshift = input.getReal("tanhshift", 2.0);
    auto dt = input.getReal("dt", 0.1);

    printfln("N = %d, v = %0.1f, h = %0.1f, tau = %0.2f, cutoff = %0.1e, max bond dim = %d", 
                                                            N, v, h, tau, truncE, maxB);

    // We will write into a file with the time-evolved energy density at all times.
    char schar2[128];
    int n2 = std::sprintf(schar2,"N_%d_v_%0.2f_h_%0.1f_tau_%0.2f_cutoff_%0.0e_maxDim_%d_tfi.dat"
                                        ,N,v,h,tau,truncE,maxB);

    std::string s2(schar2);
    std::ofstream datafile;
    datafile.open(s2); // opens the file
    if( !datafile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    datafile << "time" << " " << "en0(t)" << " " << "en(t)" << " " << "SvN" << " " << "truncerr" << " " 
                << "localEn0" << " " << "localEn" << " " << "sxsx0" << " " << "sxsx" << std::endl;

    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity=",true});

    // vectors to store energy and spin correlation data
    std::vector<double> localEn0(N-1), localEn(N-1);
    std::vector<double> sxsx0(N), sxsx(N);

    // Create the Local Energy Density Tensors
    std::vector<ITensor> LED(N-1);
    for (int b = 1; b < N; b++){
        LED[b-1] =  -4.0*sites.op("Sx",b)*sites.op("Sx",b+1);
    }    

    //magnetic field vector
    std::vector<double> hvals = hvector(N, 0.0, h, v, tau, tanhshift);

    // make Hamiltonian MPO
    auto ampo = AutoMPO(sites);
    for (int b = 1; b < N; b++){
        ampo += -4.0,"Sx", b, "Sx", b+1;
    }
    for (int b = 1; b <= N; b++){
        ampo += -2.0*(1. + hvals[b-1]),"Sz",b;
    }
    auto H = toMPO(ampo);

    // make initial state with the paramagnetic parity
    auto state = InitState(sites);
    for(int i = 1; i <= N; i++){
        state.set(i,"Up");
    }
    auto initState = MPS(state);
    PrintData(totalQN(initState));

    // sweeps
    auto sweeps = Sweeps(5); //number of sweeps is 5
    sweeps.maxdim() = 10,20,100,maxB; //gradually increase states kept
    sweeps.cutoff() = truncE; //desired truncation error

    // Find initial ground state
    auto [en,psi] = dmrg(H,initState,sweeps,{"Silent=",true});
    //calculate entanglement
    auto SvN = vonNeumannS(psi, N/2);
    //calculate local energy <psi|Hf(x)|psi>
    for (int b = 1; b < N; b++){
        psi.position(b);
        auto ket = psi(b)*psi(b+1);
        auto hi = LED[b-1];
        if(b==1){
            hi += -2.0*(1.+hvals[b-1])*sites.op("Sz",b)*sites.op("Id",b+1);
            hi += -(1.+hvals[b])  *sites.op("Id",b)*sites.op("Sz",b+1);
        }
        else if(b==N-1){
            hi += -(1.+hvals[b-1])*sites.op("Sz",b)*sites.op("Id",b+1);
            hi += -2.0*(1.+hvals[b])  *sites.op("Id",b)*sites.op("Sz",b+1);
        }
        else{
            hi += -(1.+hvals[b-1])*sites.op("Sz",b)*sites.op("Id",b+1);
            hi += -(1.+hvals[b])  *sites.op("Id",b)*sites.op("Sz",b+1);
        }
        localEn[b-1] = elt( dag(prime(ket,"Site")) * hi * ket);
    }
    //calculate spin-spin correlation
    for (int b = 1; b <= N; b++){
        sxsx[b-1] = spinspin(N/2+1,b,psi,sites);
    }

    // store data to file
    datafile << 0.0 << " " << en << " " << en << " " << SvN << " " << 0.0 << " ";
    for (int j=0; j < N-1; j++){
        datafile << localEn[j] << " ";
    }
    for (int j=0; j < N-1; j++){
        datafile << localEn[j] << " ";
    }
    for (int j = 0; j < N; j++){
        datafile << sxsx[j] << " ";
    }
    for (int j = 0; j < N; j++){
        datafile << sxsx[j] << " ";
    }
    datafile << std::endl;

    Real delta1 =  0.414490771794376*dt;
    Real delta2 = -0.657963087177503*dt;
    Real tval = 0.0;
    double finalTime = double(N) + 2.0*tau*tanhshift; // 0.5*N/c + 2*tau*tshift
    int nt = int(finalTime/dt);
    auto args = Args("Cutoff=",truncE,"MaxDim=",maxB);

    printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, en, SvN, maxLinkDim(psi));

    auto tot_truncerr=0.;
    sweeps = Sweeps(5); //number of sweeps is 5
    sweeps.maxdim() = maxB; //gradually increase states kept
    sweeps.cutoff() = truncE; //desired truncation error

    ////////////////////
    // TIME EVOLUTION //
    ////////////////////
    for (int n = 1; n <= nt; ++n){
        tval += dt;

        //update magnetic field vector
        hvals = hvector(N, tval, h, v, tau, tanhshift);
        
        // make Hamiltonian MPO
        auto ampo = AutoMPO(sites);
        for (int b = 1; b < N; b++){
            ampo += -4.0,"Sx", b, "Sx", b+1;
        }
        for (int b = 1; b <= N; b++){
            ampo += -2.0*(1. + hvals[b-1]),"Sz",b;
        }
        auto H = toMPO(ampo);

        // calculate instantaneous ground state
        // use psi as initial condition
        auto [en0,psi0] = dmrg(H,psi,sweeps,{"Silent=",true});
        
        // make TEBD gates
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

        // calculate energy <psi(t)|H(t)|psi(t)>
        auto en = innerC(psi, H, psi).real();
        //calculate entanglement entropy
        SvN = vonNeumannS(psi, N/2);
        //calculate local energy <psi0|h(t)|psi0> and <psi(t)|h(t)|psi(t)>
        for (int b = 1; b < N; b++){
            psi0.position(b);
            psi.position(b);
            auto ket0 = psi0(b)*psi0(b+1);
            auto ket = psi(b)*psi(b+1);
            auto hi = LED[b-1];
            if(b==1){
                hi += -2.0*(1.+hvals[b-1])*sites.op("Sz",b)*sites.op("Id",b+1);
                hi += -(1.+hvals[b])  *sites.op("Id",b)*sites.op("Sz",b+1);
            }
            else if(b==N-1){
                hi += -(1.+hvals[b-1])*sites.op("Sz",b)*sites.op("Id",b+1);
                hi += -2.0*(1.+hvals[b])  *sites.op("Id",b)*sites.op("Sz",b+1);
            }
            else{
                hi += -(1.+hvals[b-1])*sites.op("Sz",b)*sites.op("Id",b+1);
                hi += -(1.+hvals[b])  *sites.op("Id",b)*sites.op("Sz",b+1);
            }
            localEn0[b-1] = eltC( dag(prime(ket0,"Site")) * hi * ket0).real();
            localEn[b-1] = eltC( dag(prime(ket,"Site")) * hi * ket).real();
        }
        //calculate spin correlations <psi0|XX|psi0> and <psi(t)|XX|psi(t)>
        for (int b = 1; b <= N; b++){
            sxsx0[b-1] = spinspin(N/2+1,b,psi0,sites);
            sxsx[b-1] = spinspin(N/2+1,b,psi,sites);
        }

        datafile << tval << " " << en0 << " " << en << " " << SvN << " " << tot_truncerr << " ";
        for (int j = 0; j < N-1; j++){
            datafile << localEn0[j] << " ";
        }
        for (int j = 0; j < N-1; j++){
            datafile << localEn[j] << " ";
        }
        for (int j = 0; j < N; j++){
            datafile << sxsx0[j] << " ";
        }
        for (int j = 0; j < N; j++){
            datafile << sxsx[j] << " ";
        }
        datafile << std::endl;

        printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f \n\t maxDim = %d, total truncation error = %0.2e", 
                    tval, en, SvN, maxLinkDim(psi), tot_truncerr);

    }// for n
    
    datafile.close();

    print(" END PROGRAM. ");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}

std::vector<double> hvector(int N, double tval, double h, double v, double tau, double tanhshift)
    {
    std::vector<double> hvals(N);
    for (int b = 1; b <= N; b++){
        double f = abs(double(b-N/2) - 0.5)/v - tval;
        hvals[b-1] = h * (0.5 + 0.5*tanh(f/tau + tanhshift ));
    }
    return hvals;
}

// second order Trotter breakup of time step dt
// returns a vector of gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int N, std::vector<double> h, double dt, SiteSet sites, std::vector<ITensor> LED)
    {
    std::vector<BondGate> gates; 
    //Create the gates exp(-i*tstep/2*hterm)
    for(int i=1; i<=N; i++){
        if(i<N){
            auto hterm = LED[i-1];
            hterm += -2.0*(1+h[i-1])*op(sites,"Sz",i)*op(sites,"Id",i+1);
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
  for(int i=N; i>=1; i--){
    if(i<N){
        auto hterm = LED[i-1];
        hterm += -2.0*(1+h[i-1])*op(sites,"Sz",i)*op(sites,"Id",i+1);
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
double spinspin(int center, int b, MPS psi, SiteSet sites){

    double corrX;

    psi.position(b);
    if(b>center){ //bring site b next to the center from right
        for(int n=b-1; n>center; n--){
            auto g = BondGate(sites,n,n+1);
            auto AA = psi(n)*psi(n+1)*g.gate(); //contract over bond n
            AA.replaceTags("Site,1","Site,0");
            psi.svdBond(g.i1(), AA, Fromright); //svd from the right
            psi.position(g.i1()); //move orthogonality center to the left
        }
        auto ket = psi(center)*psi(center+1);
        auto SxSx = 4.0*sites.op("Sx",center+1)*sites.op("Sx",center);
        corrX = eltC( dag(prime(ket,"Site")) * SxSx * ket).real();
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
        auto SxSx = 4.0*sites.op("Sx",center-1)*sites.op("Sx",center);
        corrX = eltC( dag(prime(ket,"Site")) * SxSx * ket).real();
    }
    else{
        corrX = 1.;
    }

    return corrX;

}//spinspin
