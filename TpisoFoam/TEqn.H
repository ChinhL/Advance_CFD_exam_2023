   volSymmTensorField tau = laminarTransport.nu()* twoSymm(fvc::grad(U));
   //laminarTransport.nu()*(fvc::grad(U) + fvc::grad(U)().T());
    
    // Momentum predictor

   

    tmp<fvScalarMatrix> tTEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      ==
        fvm::laplacian((laminarTransport.nu()/Pr),T) 
      + ( fvc::grad(U) && tau) / Cp
     
        
        
    );
    fvScalarMatrix& TEqn = tTEqn.ref();

    TEqn.relax();

    fvOptions.constrain(TEqn);
    
    TEqn.solve();

    
