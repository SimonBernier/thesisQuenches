#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    int N = 16;
    if(argc > 1)
        N = std::stoi(argv[1]);

    printfln("N = %d", N);

    char schar1[128];
    int n1 = std::sprintf(schar1,"N_%d_tfi_gap.dat", N);

    std::string s1(schar1);
    std::ofstream enerfile;
    enerfile.open(s1); // opens the file
    if( !enerfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    enerfile << "h" << " " << "e0" << " " << "var0" << " " << "e1" << " "  << "var1" << " " << std::endl;
    
    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity=",true});

    auto state = InitState(sites);
    for(int i = 1; i <= N; i++){
        state.set(i,"Up");
    }
    auto initState = MPS(state);
    PrintData(totalQN(initState));

    std::vector<double> h(1001);
    for(int i = 0; i < int(h.size()); i++){
        h[i] = 0.5 + 1.0/(double(h.size())-1.)*double(i);
    }
    
    // calculate gap
    for(int i = 0; i < int(h.size()); i++){

        printf("h = %0.4f: ", h[i]);

        // Create the Heisenberg Hamiltonian
        auto ampo = AutoMPO(sites);
        for(int b=1; b<N; b++){
            ampo += -4.0, "Sx", b, "Sx", b+1;
        }
        for(int b=1; b<N; b++){
            ampo += -2.0*h[i], "Sz", b;
        }
        auto H = toMPO(ampo);

        //sweeps
        auto sweeps = Sweeps(10); //number of sweeps
        sweeps.maxdim() = 10,20,50,100,100,200; //gradually increase states kept
        sweeps.cutoff() = 1E-10; //desired truncation error
        sweeps.noise() = 1E-7,1E-8,1E-9,0;
        
        // Find Initial Ground State
        auto [en0,psi0] = dmrg(H,initState,sweeps,{"Silent=",true});
        auto var0 = inner(H, psi0, H, psi0) - en0*en0;
        printf("E0 = %0.8f, var0 = %0.8g, ", en0, var0);

        // find first excited state
        auto wfs = std::vector<MPS>(1);
        wfs.at(0) = psi0;
        auto [en1,psi1] = dmrg(H,wfs,initState,sweeps,{"Silent=",true,"Weight=",20.0});
        auto var1 = inner(H,psi1,H,psi1)-en1*en1;
        printfln("E1 = %0.8f, var1 = %0.8g", en1, var1);

        printfln("\t gap = %0.5g, overlap = %0.3g", en1-en0, inner(psi1,psi0));

        //store variables to energy file
        enerfile << h[i] << " " << en0 << " " << var0 << " " << en1 << " " << var1 << " " << std::endl;
    }

    enerfile.close();

    print(" END OF PROGRAM. ");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}