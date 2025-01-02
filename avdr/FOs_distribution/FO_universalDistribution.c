distribution
{
  libs (utilityFunctionObjects);
  type coded;
  name computeDistribution;
    
  executeControl    timeStep;
  executeInterval   1000;
  writeControl      timeStep;
  writeInterval     1000;
  
  codeOptions
  #{
	  -I$(LIB_SRC)/finiteVolume/lnInclude \
	  -I$(LIB_SRC)/meshTools/lnInclude \
  #};
  codeInclude
  #{
    #include "wallDist.H"
  #};
  codeData
  #{
    const scalar ratioSpecHeats  = 1.4;
    const scalar isenTotPressure = 1.0e5;

    // kobra
    // - geometry
    const scalar gamma_c   = 38.35 * Foam::constant::mathematical::pi/180;
    const scalar gamma     =      0;

    //   - lengths 
    const scalar trueChord  = 0.100;
    const scalar pitch      = 0.075;

    const scalar axialChord = trueChord*sin(gamma_c);

    // - distribution
    const int    typeDf   =  1;
    const scalar a        = -0.25;
    const scalar b        =  1.25;
    const scalar c        =  1e-3; // minimum distance from the wall for non zero fluxes

    const scalar dist_min = a*axialChord;
    const scalar dist_max = b*axialChord;   
  #};
  codeWrite
  #{  
    // create field distribution, initialize as field of zeros
    volScalarField distribution
    (
      IOobject
    (
      "distribution",
      "0",
      mesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
      mesh(),
      dimensionedScalar("distribution", dimless, 0.0)
    );

    // get wall distance of the cells
    const volScalarField  wallDistance = wallDist(mesh()).y();

    // get centroids of the cells
    const volVectorField& cCenters = mesh().C();

    // transform [x,y] of centroids to [ax,tan] system 
    // (rotate by -gamma, but we want only the first component)
    const volScalarField axialPosition = cos(gamma) * cCenters.component(vector::X) + sin(gamma) * cCenters.component(vector::Y);

    // lambda function to asses if cell of given index is at required position
    auto desiredPosition = [&](int cellIndex) -> bool
    {       
      return axialPosition[cellIndex] >= dist_min && axialPosition[cellIndex] <= dist_max && wallDistance[cellIndex] > c ; 
    } ; 

    switch(typeDf) 
    {
      case 1: 
      {
        // distribution is binary 0/1, 
        // 1 where desiredPosition == True, else 0 
        forAll(cCenters, cellI)
        {
          if ( desiredPosition( cellI ) ) 
          {
             distribution[cellI] =  1.0;  
          }
         }
        break;
      }
      case 2: 
      {
        // distribution is binary 0/1, 
        // 1 where desiredPosition == True and (isentropically) subsonic, else 0
        const volScalarField& p = mesh().lookupObject<volScalarField>("p");
        const scalar p_crit = isenTotPressure *pow(2/(1+ratioSpecHeats), ratioSpecHeats/(ratioSpecHeats-1));
        forAll(cCenters, cellI)
        {
          if ( desiredPosition( cellI ) && p[cellI] > p_crit ) 
          {
             distribution[cellI] = 1.0;  
          }
        }
        break;
      }
      case 3:
      {
        // distribution is not binary
        // (p-p_in)/(p_out-p_in) where desiredPosition == True, else 0

        const volScalarField& p = mesh().lookupObject<volScalarField>("p");
        volScalarField::Boundary& pb = const_cast<volScalarField*>(&p)->boundaryFieldRef();
        label inletId = mesh().boundaryMesh().findPatchID("inlet");
        label outletId = mesh().boundaryMesh().findPatchID("outlet");
        
        const scalarField& magSfIn  = mesh().magSf().boundaryField()[ inletId];
        const scalarField& magSfOut = mesh().magSf().boundaryField()[outletId];

        const scalar p1 = gSum(pb[ inletId]*magSfIn )/gSum(magSfIn );
        const scalar p2 = gSum(pb[outletId]*magSfOut)/gSum(magSfOut);

        forAll(cCenters, cellI)
        {
          if ( desiredPosition( cellI ) ) 
          {
             distribution[cellI] = (p[cellI] - p1)/(p2-p1); 
          }
        }
        break;
      }
    }   

    distribution.write();
    Info << "Distribution was computed and written." << endl;
  #};
}