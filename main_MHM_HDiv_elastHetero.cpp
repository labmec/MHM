
#include "pzlog.h"
#include "pzgmesh.h"
#include "pzmanvector.h"
#include "TPZMHMeshControl.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZVTKGeoMesh.h"

#include "pzelasmat.h"
#include "TPZElasticity2DHybrid.h"
#include "TPZMixedElasticityMaterial.h"
//#include "pzmat1dlin.h"
#include "TPZNullMaterial.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "meshgen.h"

#include "TPZAnalyticSolution.h"

#ifndef USING_MKL
#include "pzskylstrmatrix.h"
#endif

#include "pzgengrid.h"
#include "pzcheckgeom.h"

/// Insert material objects for the MHM Mesh solution
void InsertMaterialObjects(TPZMHMeshControl &control);

/// Insert material objects for the MHM Mesh solution
void InsertMaterialObjects(TPZMHMixedMeshControl &control);

/// Insert material objects for the MHM Mesh solution
void InsertMaterialObjects(TPZCompMesh &control);

/// Compute an approximation using an H1 approximation
TPZAutoPointer<TPZCompMesh> ComputeH1Approximation(int nelx, int nely, int porder, std::string prefix);

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mainskeleton"));
#endif


const int matInterno = 1;
const int matCoarse = 2;
const int skeleton = 4;
const int secondskeleton = 3;
const int matpressure = 6;

const int dirichlet = 0;

int const bc1=-1;
int const bc2=-2;
int const bc3=-3;
int const bc4=-4;
int const bc5=-5;

TPZAnalyticSolution *example=0;
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    TExceptionManager except;
    
    //for(int k = 1; k<2; k++){
    //    for(int j=0; j<2; j++){
    
#ifdef _AUTODIFF
    if(0){
        example = new TElasticity2DAnalytic;
        {
            TElasticity2DAnalytic *elastc_examp = dynamic_cast<TElasticity2DAnalytic *>(example);
            if(!elastc_examp) DebugStop();
            elastc_examp->fProblemType = TElasticity2DAnalytic::Etest1;
            elastc_examp->gOscilatoryElasticity = 1;
            elastc_examp->fPlaneStress = 0;
        }
    }
    
#endif
    TRunConfig Configuration;
    
    //        // Hard coded setting for Figure 15.
    //        /// numhdiv - number of h-refinements
    //        int pOrder_skel = k;
    //        int ndiv_coarse = j;
    //        int nel_coarse = 2<<ndiv_coarse;
    //        //int j_int = 2 - j;//7-j
    //        int n_div_internal = 1;//7-j (7 - ndiv_coarse);
    //        Configuration.numHDivisions = n_div_internal;
    //        /// PolynomialOrder - p-order
    //        Configuration.pOrderInternal = pOrder_skel + 0;
    //        Configuration.pOrderSkeleton = pOrder_skel;
    //        Configuration.numDivSkeleton = 0;
    //        Configuration.nelxcoarse = nel_coarse;
    //        Configuration.nelycoarse = nel_coarse;
    //        Configuration.Hybridize = 0;
    //        Configuration.Condensed = 1;
    //        Configuration.LagrangeMult = 0;
    //        Configuration.n_threads = 8;
    //        Configuration.MHM_HDiv_Elast = true;
    
    //argv: mhmelast_exec.txt
    int ndiv_internal = 8;
    int ndiv_coarse = 3;
    int ndiv_skeleton = 3;
    int pOrder_skel = 1;
    int pOrder_internal = 1;
    if(argc == 6)
    {
        ndiv_coarse = atoi(argv[1]);
        int nel_coarse = 2<<ndiv_coarse;
        Configuration.nelxcoarse = nel_coarse;
        Configuration.nelycoarse = nel_coarse/2;
        Configuration.pOrderSkeleton = atoi(argv[2]);
        pOrder_skel = Configuration.pOrderSkeleton;
        Configuration.pOrderInternal = atoi(argv[3]);
        pOrder_internal = Configuration.pOrderInternal;
        Configuration.numHDivisions = atoi(argv[4]);
        ndiv_internal = Configuration.numHDivisions;
        Configuration.numDivSkeleton = atoi(argv[5]);
        ndiv_skeleton = Configuration.numDivSkeleton;
        Configuration.Hybridize = 0;
        Configuration.Condensed = 1;
        Configuration.LagrangeMult = 0;
        Configuration.n_threads = 8;
        Configuration.MHM_HDiv_Elast = true;
        
    }
    else if (argc == 3)
    {
        Configuration.nelxcoarse = atoi(argv[1]);
        Configuration.nelycoarse = Configuration.nelxcoarse;
        Configuration.pOrderSkeleton = atoi(argv[2]);
        Configuration.pOrderInternal = Configuration.pOrderSkeleton+1;
    }
    else
    {
        std::cout << "Executing using internal hard-code variables \n";
        // Hard coded setting for Figure 15.
        /// numhdiv - number of h-refinements
        int k = 1;
        int j = 3;
        int nel_coarse = 2<<ndiv_coarse;
        //int j_int = 2 - j;//7-j
        Configuration.numHDivisions = ndiv_internal;
        /// PolynomialOrder - p-order
        Configuration.pOrderInternal = pOrder_internal;
        Configuration.pOrderSkeleton = pOrder_skel;
        Configuration.numDivSkeleton = 3;
        Configuration.nelxcoarse = nel_coarse;
        Configuration.nelycoarse = nel_coarse/2;
        Configuration.Hybridize = 0;
        Configuration.Condensed = 1;
        Configuration.LagrangeMult = 0;
        Configuration.n_threads = 8;
        Configuration.MHM_HDiv_Elast = true;
    }
    
    
    if(0)
    {
        /// Compute an approximation using an H1 approximation
        int nx = 40;
        int porder = 2;
        std::string prefix = "H1Aprox";
        
        TPZAutoPointer<TPZCompMesh> cmeshH1 = ComputeH1Approximation(nx, nx, porder, prefix);
    }
    
    
    // to avoid singular internal matrices
    if (Configuration.numDivSkeleton == Configuration.numHDivisions && Configuration.pOrderInternal <= Configuration.pOrderSkeleton) {
        Configuration.pOrderInternal = Configuration.pOrderSkeleton+1;
    }
    
    
    //TPZGeoMesh *gmesh = 0;
    TPZStack<int64_t> coarseindices;
    
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.);
    int ndiv = Configuration.numHDivisions;
    //gmesh = MalhaGeomFredQuadrada(Configuration.nelxcoarse, Configuration.nelycoarse, x0, x1, coarseindices, ndiv);



#ifdef MACOSX
    std::string path1 = "../young_modulus.txt";
    std::string path2 = "../Poisson_coef.txt";
    std::string path3 = "../Poisson_coef.txt";
    path1 = "../Data_13_Set/VE.txt";
    path2 = "../Data_13_Set/VPoisso.txt";
    path3 = "../Data_13_Set/Vden.txt";

#else
    std::string path1 = "young_modulus.txt";
    std::string path2 = "Poisson_coef.txt";
    std::string path3 = "../Poisson_coef.txt";
    path1 = "Data_13_Set/VE.txt";
    path2 = "Data_13_Set/VPoisso.txt";
    path3 = "Data_13_Set/Vden.txt";
#endif

    int flag = ReadFromFile( ElastCoef, path1);
    if (!flag)
    {
        DebugStop();
    }
    ElastCoef *= 1.e-6;
    {
        std::ofstream out ("Elast.nb");
//        ElastCoef.Print("Elast",out,EMathematicaInput);
    }
    flag = ReadFromFile( PoissonCoef, path2);
    {
        std::ofstream out ("Poisson.nb");
//        PoissonCoef.Print("Poisson",out,EMathematicaInput);
    }
    if (!flag)
    {
        DebugStop();
    }
    flag = ReadFromFile( DensityCoef, path3);
    if (!flag)
    {
        DebugStop();
    }
    {
        std::ofstream out ("Density.nb");
//        DensityCoef.Print("Rho",out,EMathematicaInput);
    }


    TPZVec<int> nx(2,5);
    nx[0] =2;
    nx[1] =1;
    max_x = 10000.;
    max_y = 4500.;
    
    x1[0]=max_x;
    x1[1]=max_y;
    x1[2]=0.0;

    TPZGenGrid grid(nx,x0,x1);
    grid.SetRefpatternElements(true);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    //   grid.SetElementType(ETriangle);
    grid.Read(gmesh);
    grid.SetBC(gmesh,4,-1);
    grid.SetBC(gmesh,5,-2);
    grid.SetBC(gmesh,6,-3);
    grid.SetBC(gmesh,7,-4);

    TPZCheckGeom check(gmesh);
    int nrefine = ndiv_internal;
    check.UniformRefine(nrefine);
    
    {
        REAL Omega = max_x*max_y/2.;
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->HasSubElement() || gel->Dimension() != 2) continue;
            REAL area = gel->SideArea(gel->NSides()-1);
            REAL ratio = Omega/area;
            ratio = sqrt(ratio);
            std::cout << "the element side is " << ratio << " fraction of the total length\n";
            break;
        }
    }

    int nivel = ndiv_coarse;
    GetcoarseID(gmesh, nivel, coarseindices);

    TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
    TPZAutoPointer<TPZMHMeshControl> MHM;
    
    
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmesh);
        TPZMHMixedMeshControl *mhm = new TPZMHMixedMeshControl(gmeshauto);

        mhm->DefinePartitionbyCoarseIndices(coarseindices);
        MHM = mhm;
        TPZMHMixedMeshControl &meshcontrol = *mhm;
        
        meshcontrol.SetLagrangeAveragePressure(Configuration.LagrangeMult);
        
        InsertMaterialObjects(meshcontrol);
        meshcontrol.SetInternalPOrder(Configuration.pOrderInternal);
        meshcontrol.SetSkeletonPOrder(Configuration.pOrderSkeleton);

        meshcontrol.DivideSkeletonElements(Configuration.numDivSkeleton);
        if(Configuration.Hybridize)
        {
            meshcontrol.SetHybridize(true);
        }
        
        bool substructure = true;
        if (Configuration.Condensed == 0) {
            substructure = false;
        }
        meshcontrol.BuildComputationalMesh(substructure);
        if(1)
        {
            std::ofstream file("GMeshControl.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file);
        }
        
#ifdef PZDEBUG
        if(0)
        {
            std::ofstream out("MHMMeshControl.txt");
            meshcontrol.Print(out);
        }
#endif
        
        Configuration.fGlobalSystemSize = meshcontrol.fGlobalSystemSize;
        Configuration.fGlobalSystemWithLocalCondensationSize = meshcontrol.fGlobalSystemWithLocalCondensationSize;
        Configuration.fNumeq = meshcontrol.fNumeq;
        std::cout << "MHM Computational meshes created\n";
#ifdef PZDEBUG
        if(0)
        {
            std::ofstream gfile("geometry.txt");
            gmesh->Print(gfile);
            
            std::ofstream out_mhm("MHM_hdiv.txt");
            meshcontrol.CMesh()->Print(out_mhm);
            
        }
#endif
        std::cout << "Number of equations MHM equals " << MHM->CMesh()->NEquations() << std::endl;
        
    }
    std::string configuration;
    
    {
        std::stringstream sout;
        sout << "H" << Configuration.numHDivisions << "-P" << Configuration.pOrderInternal;
        configuration = sout.str();
    }
    
    // compute the MHM solution
        SolveProblem(MHM->CMesh(), MHM->GetMeshes(), example, "MHMElast_Hdiv", Configuration);
    
    //    }//fim_k
    //}//fim_j
    return 0;
}



void InsertMaterialObjects(TPZMHMixedMeshControl &control)
{
    TPZCompMesh &cmesh = control.CMesh();
    /// criar materiais
    int dim = cmesh.Dimension();
    TElasticity2DAnalytic *elastc_examp = dynamic_cast<TElasticity2DAnalytic *>(example);
    //if(!elastc_examp) DebugStop();
    control.SetProblemType(TPZMHMeshControl::EElasticity2D);
    
    STATE Young = 1000., nu = 0.3, fx = 0., fy = 0.;
    int planeStress = 0; //* @param plainstress = 1 \f$ indicates use of plain stress
    TPZMixedElasticityMaterial *material1 = new TPZMixedElasticityMaterial(matInterno,Young, nu, fx,fy,planeStress,dim);

    TPZAutoPointer<TPZFunction<STATE>> funcelast = new TPZDummyFunction<STATE> (ForcingFunction,1);
    material1->SetForcingFunction(funcelast);
    
    TPZAutoPointer<TPZFunction<STATE>> func = new TPZDummyFunction<STATE> (funcE2, 1);
    material1->SetElasticityFunction(func);


    control.fMaterialIds.insert(matInterno);
    if(example) material1->SetForcingFunction(example->ForcingFunction());
    
    if(elastc_examp)
    {
        material1->SetElasticityFunction(elastc_examp->ConstitutiveLawFunction());
        if(elastc_examp->fPlaneStress)
        {
            material1->SetPlaneStress();
        }
        else
        {
            material1->SetPlaneStrain();
        }
    }
    
    TPZMaterial * mat1(material1);
    
    TPZNullMaterial *materialCoarse = new TPZNullMaterial(matCoarse);
    materialCoarse->SetNStateVariables(2);
    materialCoarse->SetDimension(1);
    
    cmesh.InsertMaterialObject(materialCoarse);
    //    materialCoarse = new TPZMat1dLin(skeleton);
    //    materialCoarse->SetMaterial(xk, xc, xb, xf);
    //    cmesh.InsertMaterialObject(materialCoarse);
    //    materialCoarse = new TPZMat1dLin(secondskeleton);
    //    materialCoarse->SetMaterial(xk, xc, xb, xf);
    //    cmesh.InsertMaterialObject(materialCoarse);
    //    materialCoarse = new TPZMat1dLin(matpressure);
    //    materialCoarse->SetMaterial(xk, xc, xb, xf);
    //    cmesh.InsertMaterialObject(materialCoarse);
    
    
    
    //    REAL diff = -1.;
    //    REAL conv = 0.;
    //    TPZVec<REAL> convdir(3,0.);
    //    REAL flux = 8.;
    //
    //    material1->SetParameters(diff, conv, convdir);
    //    material1->SetInternalFlux( flux);
    
    cmesh.InsertMaterialObject(mat1);
    
    
    control.SwitchLagrangeMultiplierSign(false);
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    //BC -1
    /// mixed : sign 1x = 0 which means tau_xy = 0
    /// as there is no entry in y , u.1y =0 weakly
    val1(0,0) = 0;
    val2.Zero();
    val1(0,0) = 0;
    val1(1,1) = 0;
    val2(1,0) = 0.;
    TPZMaterial * BCondD1 = material1->CreateBC(mat1, bc1, 0, val1, val2);
    if(example) BCondD1->SetForcingFunction(example->Exact());
    cmesh.InsertMaterialObject(BCondD1);
    control.fMaterialBCIds.insert(bc1);
    //BC -2
    // Neumann condition
    // sigma_n = (10,0) which means sigx = 10, tau_xy = 0
    val1.Zero();
    val2.Zero();
    val1(1,1) = 0.;
    val2(0,0) = 0;
    TPZMaterial * BCondD2 = material1->CreateBC(mat1, bc2,1, val1, val2);
    if(example) BCondD2->SetForcingFunction(example->Exact());
    cmesh.InsertMaterialObject(BCondD2);
    control.fMaterialBCIds.insert(bc2);
    
    //BC -3
    // Neumann condition
    // applied force is (0,20)
    // sigma_n = (0,20)
    val1.Zero();
    val2.Zero();
    val1(0,0) = 0;
    val2(1,0) = 0.;
    TPZMaterial * BCondD3 = material1->CreateBC(mat1, bc3,1, val1, val2);
    if(example) BCondD3->SetForcingFunction(example->Exact());
    cmesh.InsertMaterialObject(BCondD3);
    control.fMaterialBCIds.insert(bc3);
    
    //BC -4
    // mixed condition
    // sigma_n .1_y = 0 -> tau_xy = 0
    // u_x = 0 weakly
    val1.Zero();
    val1(1,1) = 0;
    val2.Zero();
    TPZMaterial * BCondD4 = material1->CreateBC(mat1, bc4,1, val1, val2);
    if(example) BCondD4->SetForcingFunction(example->Exact());
    cmesh.InsertMaterialObject(BCondD4);
    control.fMaterialBCIds.insert(bc4);
    
    //BC -5: dirichlet nulo
    TPZMaterial * BCondD5 = material1->CreateBC(mat1, bc5,dirichlet, val1, val2);
    cmesh.InsertMaterialObject(BCondD5);
    control.fMaterialBCIds.insert(bc5);
    
}
void InsertMaterialObjects(TPZCompMesh &cmesh)
{
    /// criar materiais
    //    int dim = cmesh.Dimension();
    STATE Young = 1000., nu = 0.3, fx = 0., fy = 0.;
    TPZElasticityMaterial *material1 = new TPZElasticityMaterial(matInterno,Young,nu,fx,fy);
    material1->SetPlaneStrain();
    
    if(example) material1->SetForcingFunction(example->ForcingFunction());
    TElasticity2DAnalytic *example_elast = dynamic_cast<TElasticity2DAnalytic *>(example);
    if(example_elast) material1->SetElasticityFunction(example_elast->ConstitutiveLawFunction());
    
    TPZMaterial * mat1(material1);
    
    
    //    REAL diff = -1.;
    //	REAL conv = 0.;
    //	TPZVec<REAL> convdir(3,0.);
    //	REAL flux = 8.;
    //
    //	material1->SetParameters(diff, conv, convdir);
    //	material1->SetInternalFlux( flux);
    
    cmesh.InsertMaterialObject(mat1);
    
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    //BC -1
    val1(0,0) = 0.;
    val2.Zero();
    val1(0,0) = 0;
    val1(1,1) = 0;
    TPZMaterial * BCondD1 = material1->CreateBC(mat1, bc1,dirichlet, val1, val2);
    if(example) BCondD1->SetForcingFunction(example->Exact());
    cmesh.InsertMaterialObject(BCondD1);
    
    //BC -2
    val1.Zero();
    val2(0,0) = 10.;
    TPZMaterial * BCondD2 = material1->CreateBC(mat1, bc2,dirichlet, val1, val2);
    if(example) BCondD2->SetForcingFunction(example->Exact());
    cmesh.InsertMaterialObject(BCondD2);
    
    //BC -3
    val1.Zero();
    val2.Zero();
    TPZMaterial * BCondD3 = material1->CreateBC(mat1, bc3,dirichlet, val1, val2);
    if(example) BCondD3->SetForcingFunction(example->Exact());
    cmesh.InsertMaterialObject(BCondD3);
    
    //BC -4
    val1(0,0) = 0;
    val1(1,1) = 1.e9;
    val2(0,0) = -1.;
    TPZMaterial * BCondD4 = material1->CreateBC(mat1, bc4,dirichlet, val1, val2);
    if(example) BCondD4->SetForcingFunction(example->Exact());
    cmesh.InsertMaterialObject(BCondD4);
    
    //BC -5: dirichlet nulo
    TPZMaterial * BCondD5 = material1->CreateBC(mat1, bc5,dirichlet, val1, val2);
    cmesh.InsertMaterialObject(BCondD5);
}

/// Compute an approximation using an H1 approximation
TPZAutoPointer<TPZCompMesh> ComputeH1Approximation(int nelx, int nely, int porder, std::string prefix)
{
    TPZGeoMesh *gmesh = 0;
    TPZVec<int64_t> coarseindices;
    
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.);
    x1[2] = 0.;
    int ndiv = 0;
    gmesh = MalhaGeomFredQuadrada(nelx, nely, x0, x1, coarseindices, ndiv);
    
    gmesh->SetDimension(2);
    TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
    
    TPZAutoPointer<TPZCompMesh> cmeshauto = new TPZCompMesh(gmeshauto);
    cmeshauto->SetDimModel(2);
    InsertMaterialObjects(cmeshauto);
    
    cmeshauto->SetAllCreateFunctionsContinuous();
    
    cmeshauto->AutoBuild();
    
    
    //calculo solution
    bool shouldrenumber = true;
    TPZAnalysis an(cmeshauto,shouldrenumber);
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmeshauto.operator->());
    strmat.SetNumThreads(8);
    
#else
    TPZSkylineStructMatrix strmat(cmeshauto.operator->());
    strmat.SetNumThreads(0);
#endif
    
    
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();
    if(0)
    {
        std::string filename = prefix;
        filename += "_Global.nb";
        std::ofstream global(filename.c_str());
        TPZAutoPointer<TPZStructMatrix> strmat = an.StructMatrix();
        an.Solver().Matrix()->Print("Glob = ",global,EMathematicaInput);
        an.Rhs().Print("Rhs = ",global,EMathematicaInput);
    }
    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
    an.LoadSolution(); // compute internal dofs
                       //    an.Solution().Print("sol = ");
    
    
#ifdef PZDEBUG
    {
        std::ofstream out(prefix+"_MeshWithSol.txt");
        cmeshauto->Print(out);
    }
#endif
    
    std::string configuration;
    {
        std::stringstream sout;
        sout << nelx << "-" << nely;
        configuration = sout.str();
    }
    //    TPZBuildMultiphysicsMesh::TransferFromMeshes(cmeshes, an.Mesh());
    //    for (int i=0; i<cmeshes.size(); i++) {
    //        cmeshes[i]->Solution().Print("sol = ");
    //    }
    //    cmeshes[0]->Solution().Print("solq = ");
    //    cmeshes[1]->Solution().Print("solp = ");
    std::stringstream sout;
    sout << prefix << "Approx-";
    sout << configuration << ".vtk";
    std::string plotfile = sout.str();
    std::cout << "plotfile " << plotfile.c_str() << std::endl;
    TPZStack<std::string> scalnames,vecnames;
    TPZMaterial *mat = cmeshauto->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }
    if (mat->NStateVariables() == 2)
    {
        scalnames.Push("SigmaX");
        scalnames.Push("SigmaY");
        scalnames.Push("TauXY");
        vecnames.Push("Displacement");
    }
    else if(mat->NStateVariables() == 1)
    {
        scalnames.Push("Pressure");
        vecnames.Push("Flux");
        vecnames.Push("Derivative");
    }
    an.DefineGraphMesh(cmeshauto->Dimension(), scalnames, vecnames, plotfile);
    int resolution = 2;
    an.PostProcess(resolution,cmeshauto->Dimension());
    
#ifdef _AUTODIFF
    std::cout << "Computing errors\n";
    int64_t neq = cmeshauto->NEquations();
    an.SetExact(example->ExactSolution());
    TPZVec<REAL> errors(3,0.);
    an.PostProcessError(errors);
    std::cout << "Errors computed " << errors << std::endl;
    
    std::stringstream filename;
    filename << prefix << "Errors.txt";
    std::ofstream out (filename.str(),std::ios::app);
    out << "nelx " << nelx << " nely " << nely << " porder " << porder << " neq " << neq <<  " Energy " << errors[0] << " L2 " << errors[1] << " H1 " << errors[2] << std::endl;
#endif
    return cmeshauto;
    
}

