#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "pzpoisson3d.h"
//#include "pzhybridpoisson.h"
#include "pzpoisson3dreferred.h"
#include "mixedpoisson.h"
#include "pzelasmat.h"
#include "pzelasthybrid.h"
#include "pzmat1dlin.h"
#include "TPZVecL2.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZLagrangeMultiplier.h"


#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "TPZCompMeshTools.h"
#include "pzcondensedcompel.h"
#include "pzfunction.h"
#include "pzgraphmesh.h"
#include "pzfmatrix.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"
#include "pzvisualmatrix.h"
#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "pzcheckgeom.h"

#include "TPZMHMeshControl.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZMHMixedHybridMeshControl.h"

#include "meshgen.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>



using namespace std;

struct TRunConfig;

TPZGeoMesh *MalhaGeomFred(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, const std::string quad, const std::string triangle, TPZVec<int64_t> &coarseindices, int ndiv);

/// Create a Refinement Pattern that divides a quadrilateral by two triangles
TPZAutoPointer<TPZRefPattern> DivideQuadbyTriangles(const std::string refpatname);

/// Create a Refinement Patterns that divides a triangle into nine triangles
TPZAutoPointer<TPZRefPattern> DivideTriangleby9Triangles(const std::string refpatname);

/// Insert material objects for the MHM Mesh solution
void InsertMaterialObjects(TPZMHMeshControl &control);
/// Insert material objects for the MHM-H(div) solution
void InsertMaterialObjects(TPZMHMixedMeshControl &control);
/// Insert material objects for the MHM-H(div) solution
void InsertMaterialObjects(TPZMHMixedHybridMeshControl &control);

/// Compute the differences at the submesh level
void ComputeDifferencesBySubmesh(TRunConfig &config, TPZMHMeshControl &MHM, TPZMHMixedMeshControl &MHMixed, const std::string &filename);

/// Create a reference geometric mesh starting with nelx by nely domains
// called in the main program
TPZGeoMesh *CreateReferenceGMesh(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, int numref);

/// compute the reference solution and return created mesh
// call CreateReferenceCMesh
TPZCompMesh *ComputeReferenceSolution(TPZGeoMesh *gmesh, int porder, TPZVec<TPZCompMesh *> &meshvec);

/// create the computational mesh of the reference solution
// call CreateHDivMHMMesh and CreatePressureMHMMesh
TPZCompMesh *CreateReferenceCMesh(TPZGeoMesh *gmesh, TPZVec<TPZCompMesh *> &meshvec, int porder);

/// Create an HDiv mesh used as a reference mesh
TPZCompMesh * CreateHDivMHMMesh(TPZGeoMesh * gmesh, int porder);
/// Create a pressure mesh for the reference mesh
TPZCompMesh * CreatePressureMHMMesh(TPZGeoMesh * gmesh, int porder, int dimension);

/// Create a multiphysics mesh from an H(div) and Pressure mesh
/// Called by the method which creates the reference computacional mesh
TPZCompMesh * CreateHDivPressureMHMMesh(TPZVec<TPZCompMesh * > &cmesh);


/// Analise the regularity of the subdomain problems
void AnalyseRegularity(const TPZVec<int> &pos0,const TPZVec<int> &nelx, TPZVec<int> &nsub, TPZFMatrix<REAL> &lowestexp);

/// Print the elements with geometric information and connect values
void PrintElements(TPZCompMesh *cmesh, std::ostream &out);

/// copy the solution between one computation mesh to the other assuming the geometric elements match
void CopySolution(TPZCompMesh *from, TPZCompMesh *to);

/// unwrap de TPZCondensedCompel and TPZElementGroup elements
void UnwrapMesh(TPZCompMesh *cmesh);

/// function that returns the permeability for a given coordinate
void Permeability(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &diff);

/// function that randomly refines some elements
void RandomRefine(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices, int nref);

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mainskeleton"));
#endif


const int matInterno = 1;
const int matCoarse = 2;
const int skeleton = 4;
const int secondskeleton = 3;
const int matpressure = 6;

const int dirichlet = 0;
const int neumann = 1;

int const bc1=-1;
int const bc2=-2;
int const bc3=-3;
int const bc4=-4;
int const bc5=-5;

static void DirichletValidacao(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {   ///Jorge 2017 , TPZFMatrix<STATE> &gradres){
    result[0] = loc[0]+1;
}


/**
 Extract arguments for <inputs>: n_dx n_dy k_skeleton m_div k_subelement

 @param argc number or arguments
 @param argv vector of arguments
 @return std::vector containing integers
 */
std::vector <int64_t> ExtractArguments(int argc, char *argv[]);

TPZFMatrix<REAL> gPorous(500,100,0.);

TAnalyticSolution *example = 0;

int main(int argc, char *argv[])
{
    std::vector <int64_t> parameters = ExtractArguments(argc, argv);
    TExceptionManager except;
    
#ifdef _AUTODIFF
//    example = new TLaplaceExampleSmooth; //  Problem 1
    example = new TLaplaceExample1; // Problem 2
#endif
    
    TRunConfig Configuration;
    
    /// computation type :
    // (0) - compute reference mesh
    // (1) - compute MHM H1 mesh and compute MHM(div) mesh
    int ComputationType = 1;
    /// numhdiv - number of h-refinements
    Configuration.numHDivisions = 2;
    /// PolynomialOrder - p-order
    Configuration.pOrderInternal = 1;
    
    Configuration.Hybridize = 0;
    Configuration.Condensed = 1;
    Configuration.n_threads = 32;
    
    Configuration.pOrderSkeleton = 1;
    Configuration.numDivSkeleton = 0;
    
    TPZManVector<REAL,3> x0(2,0.),x1(2,0.);
    // for using the aligned mesh
    x0[0] = 0;
    if (!example)
    {
        int nelxref = 16;
        int nelyref = 4;
        Configuration.nelxcoarse = nelxref;
        Configuration.nelycoarse = nelyref;
    }
    else
    {
        int nelxref = 2;
        int nelyref = 2;
        Configuration.nelxcoarse = nelxref;
        Configuration.nelycoarse = nelyref;
    }

    int64_t n_par = parameters.size();
    if (n_par != 0 && n_par == 5)
    {
        std::cout << "****************************************" << std::endl;
        std::cout << " Executing using command line arguments " << std::endl;
        std::cout << "n_dx          = " << parameters[0] << std::endl;
        std::cout << "n_dy          = " << parameters[1] << std::endl;
        std::cout << "k_skeleton    = " << parameters[2] << std::endl;
        std::cout << "m_div         = " << parameters[3] << std::endl;
        std::cout << "k_subelement  = " << parameters[4] << std::endl;
        std::cout << "****************************************" << std::endl;
        // <inputs>: n_dx n_dy k_skeleton m_div k_subelement
        Configuration.nelxcoarse = parameters[0];
        Configuration.nelycoarse = parameters[1];
        Configuration.numHDivisions = parameters[3];
        Configuration.pOrderInternal = parameters[4];
        Configuration.numDivSkeleton = 0;
        Configuration.pOrderSkeleton = parameters[2];
        
    }else{
        std::cout << "Executing using internal hard-code variables \n";
    }
    
    // to avoid singular internal matrices
    if (Configuration.numHDivisions == 0 && Configuration.pOrderInternal <= Configuration.pOrderSkeleton) {
        Configuration.pOrderInternal = Configuration.pOrderSkeleton+1;
    }
    // Hereogeneous flow
    x1[0] = x0[0]+1.0;
    x1[1] = x0[1]+0.2;
    
    if(example)
    {
        x0.Fill(0.);
        x1.Fill(1.);
    }

    HDivPiola = 1;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    {
#ifdef MACOSX
        std::ifstream pores("../porous_scaled.txt");
#else
        std::ifstream pores("porous_scaled.txt");
#endif
        for (int j=0; j<100; j++) {
            for (int i=0; i<500; i++) {
                pores >> gPorous(i,j);
                if (!pores) {
                    DebugStop();
                }
            }
        }
    }
    // tototo
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
//    gRefDBase.InitializeUniformRefPattern(ECube);
    
    TPZManVector<TPZCompMesh *,2> ReferenceMeshVec(2,0);
    TPZGeoMesh *ReferenceGMesh = 0;
    TPZCompMesh *ReferenceCMesh = 0;

//    gRefDBase.InitializeRefPatterns();
    
    if(ComputationType == 0)
    {
        // generate the reference solution, save it on disk and exit
        int nelx = Configuration.nelxcoarse, nely = Configuration.nelycoarse;
        int numref = 1;
        TPZGeoMesh *gmesh = CreateReferenceGMesh(nelx, nely, x0, x1, numref);
        TPZManVector<TPZCompMesh *,2> meshvec(2);
        int porder = 1;
        TPZCompMesh *cmesh = ComputeReferenceSolution(gmesh,porder,meshvec);
        TPZPersistenceManager::OpenWrite("Ref.bin");
        gmesh->ResetReference();
        TPZPersistenceManager::WriteToFile(gmesh);
        TPZPersistenceManager::WriteToFile(meshvec[0]);
        TPZPersistenceManager::WriteToFile(meshvec[1]);
        TPZPersistenceManager::CloseWrite();
        
        //       cmesh->Write(meshfile, false);
        if(0)
        {
            ofstream out("gmesh1.txt");
            gmesh->Print(out);
        }
        if(0)
        {
            ofstream out1("cmeshwrite0.txt");
            meshvec[0]->Print(out1);
            ofstream out2("cmeshwrite1.txt");
            meshvec[1]->Print(out2);
        }
        UnwrapMesh(cmesh);

        if(0)
        {
            ofstream out1("mfmeshwrite.txt");
            cmesh->Print(out1);
        }

        delete cmesh;
        delete meshvec[0];
        delete meshvec[1];
        delete gmesh;
        exit(0);
    }
    if(0)
    {
        // read the reference solution from the file
        TPZPersistenceManager::OpenRead("Ref.bin");
        TPZGeoMesh *gmesh = dynamic_cast<TPZGeoMesh *>(TPZPersistenceManager::ReadFromFile());
        ReferenceGMesh = gmesh;
        if(0)
        {
            ofstream out("gmesh2.txt");
            gmesh->Print(out);
        }
        TPZManVector<TPZCompMesh *,2> meshvec(2);
        meshvec[0] = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::ReadFromFile());
        meshvec[1] = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::ReadFromFile());
        TPZPersistenceManager::CloseRead();
        ReferenceMeshVec = meshvec;
        if(0)
        {
            ofstream out1("cmeshread0.txt");
            meshvec[0]->Print(out1);
            ofstream out2("cmeshread1.txt");
            meshvec[1]->Print(out2);
        }
        TPZCompMesh *cmesh = CreateHDivPressureMHMMesh(meshvec);
        for (int i=0; i<10; i++) {
            TPZMaterial *mat = cmesh->FindMaterial(i);
            if (mat) {
                TPZMixedPoisson *mixed = dynamic_cast<TPZMixedPoisson *>(mat);
                if (!mixed) {
                    DebugStop();
                }
                TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Permeability);
                dummy->SetPolynomialOrder(0);
                TPZAutoPointer<TPZFunction<STATE> > func(dummy);
                mixed->SetPermeabilityFunction(func);
            }
        }

        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, cmesh);
        ReferenceCMesh = cmesh;
        if(0)
        {
            std::string plotfile("referencesolutionRead.vtk");
            TPZStack<std::string> scalnames,vecnames;
            scalnames.Push("Pressure");
            scalnames.Push("Permeability");
            vecnames.Push("Derivative");
            vecnames.Push("Flux");
            TPZAnalysis an(cmesh,false);
            if(0)
            {
                ofstream out1("mfmeshread.txt");
                cmesh->Print(out1);
            }
            an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotfile);
            int resolution = 0;
            an.PostProcess(resolution,cmesh->Dimension());
        }
    }
    TPZGeoMesh *gmesh = 0;
    TPZManVector<int64_t> coarseindices;
    if(0)
    {
        // original research paper - the mesh was not aligned with the heterogeneities
        std::string quad = "QuadByTriangles";
        std::string triangle = "TriangleBy9Triangles";
        TPZAutoPointer<TPZRefPattern> refpatquad = DivideQuadbyTriangles(quad);
        quad = refpatquad->Name();
        TPZAutoPointer<TPZRefPattern> refpattriangle = DivideTriangleby9Triangles(triangle);
        int nelx = 15;
        int nely = 5;
        Configuration.nelxcoarse = nelx;
        Configuration.nelycoarse = nely;
        int ndiv = Configuration.numHDivisions;
        gmesh = MalhaGeomFred(nelx, nely, x0, x1, quad, triangle, coarseindices, ndiv);
        {
            std::ofstream out("DiffResults.nb",std::ios::app);
            out << "(* Running triangular mesh with subdomains " << nelx << " " << nely << " *)\n";
        }
    }
    else if(!example)
    {
        // verifying differences between the MHM-original and MHM with mixed approximations
        int nelx = Configuration.nelxcoarse;
        int nely = Configuration.nelycoarse;
        {
            std::ofstream out("DiffResults.nb",std::ios::app);
            out << "(* Running quadrilateral mesh with numsubdomains " << nelx << ", " << nely << " *)\n";
        }
        /// Analise the regularity of the subdomain problems
        TPZManVector<int,3> nelvec(2),nsub(2);
        nelvec[0] = Configuration.nelxcoarse;
        nelvec[1] = Configuration.nelycoarse;
        nsub[0] = nelx;
        nsub[1] = nely;
        TPZFMatrix<REAL> lowestexp;
        TPZManVector<int,2> pos0(2,0);
        pos0[0] = 100;

        AnalyseRegularity(pos0, nelvec,  nsub,  lowestexp);

        VisualMatrixVTK(lowestexp, "regularity.vtk");
        {
            std::ofstream out("regularity.nb");
            lowestexp.Print("Regularity=",out,EMathematicaInput);
        }

        int ndiv = Configuration.numHDivisions;
        gmesh = MalhaGeomFredQuadrada(nelx, nely, x0, x1, coarseindices, ndiv);
//        RandomRefine(gmesh, coarseindices,1);
    }
    else
    {
        {
            std::ofstream out("DiffResults.nb",std::ios::app);
            out << "(* Running quadrilateral mesh with Config { ";
            Configuration.MathematicaInlinePrint(out);
            out << "} *)\n";
        }
        int ndiv = Configuration.numHDivisions;
        gmesh = MalhaGeomFredQuadrada(Configuration.nelxcoarse, Configuration.nelycoarse, x0, x1, coarseindices, ndiv);
//        RandomRefine(gmesh, coarseindices,1);
        
    }
    
    TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
    TPZAutoPointer<TPZMHMeshControl> MHM;
    TPZAutoPointer<TPZMHMixedMeshControl> MHMixed;
    
    std::stringstream MHMPref, MHMMixedPref;


    if(1)
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmesh);
        TPZMHMeshControl *mhm = new TPZMHMeshControl(gmeshauto);
        {
            std::set<int> matids;
            matids.insert(1);
            mhm->fMaterialIds = matids;
            matids.clear();
            matids.insert(-1);
            matids.insert(-2);
            matids.insert(-3);
            matids.insert(-4);
            mhm->fMaterialBCIds = matids;
        }
        mhm->DefinePartitionbyCoarseIndices(coarseindices);
        MHMPref << "MHM";
        MHM = mhm;
        TPZMHMeshControl &meshcontrol = *mhm;
        MHM->SwitchLagrangeMultiplierSign(true);

        if (Configuration.LagrangeMult) {
            meshcontrol.SetLagrangeAveragePressure(true);
        }
        
        InsertMaterialObjects(*mhm);

        meshcontrol.SetInternalPOrder(Configuration.pOrderInternal);
        meshcontrol.SetSkeletonPOrder(Configuration.pOrderSkeleton);
        
        meshcontrol.DivideSkeletonElements(Configuration.numDivSkeleton);
        if (Configuration.Hybridize)
        {
            meshcontrol.SetHybridize(true);
        }
        
        bool substructure = (bool) Configuration.Condensed;
        meshcontrol.BuildComputationalMesh(substructure);
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream file("GMeshControl.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file,true);
        }
#endif
#ifdef PZDEBUG
        if(0)
        {
            std::ofstream out("MHMMeshControl.txt");
            meshcontrol.Print(out);
        }
#endif

        std::cout << "MHM Computational meshes created\n";
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream gfile("geometry.txt");
            gmesh->Print(gfile);

            std::ofstream out_mhm("MHM_hybrid.txt");
            meshcontrol.CMesh()->Print(out_mhm);

        }
#endif
        std::cout << "Number of equations MHM equals " << MHM->CMesh()->NEquations() << std::endl;
    
    }
    
    if(1)
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmesh);
        TPZMHMixedMeshControl *mhm = new TPZMHMixedMeshControl(gmeshauto);
        // criam-se apenas elementos geometricos
        mhm->DefinePartitionbyCoarseIndices(coarseindices);
        MHMMixedPref << "MHMixed";
        MHMixed = mhm;
        TPZMHMixedMeshControl &meshcontrol = *mhm;
        {
            std::set<int> matids;
            matids.insert(1);
            mhm->fMaterialIds = matids;
            matids.clear();
            matids.insert(-1);
            matids.insert(-2);
            matids.insert(-3);
            matids.insert(-4);
            mhm->fMaterialBCIds = matids;
        }

        
        InsertMaterialObjects(*mhm);
        
        meshcontrol.SetInternalPOrder(Configuration.pOrderInternal);
        meshcontrol.SetSkeletonPOrder(Configuration.pOrderSkeleton);
        
        meshcontrol.DivideSkeletonElements(Configuration.numDivSkeleton);

        if (Configuration.Hybridize)
        {
            meshcontrol.SetHybridize(true);
        }
        
        bool substructure = (bool) Configuration.Condensed;
        meshcontrol.BuildComputationalMesh(substructure);
#ifdef PZDEBUG
        if(0)
        {
            std::ofstream file("GMeshControlHDiv.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file);
        }
#endif
#ifdef PZDEBUG
        if(0)
        {
            std::ofstream out("MixedMeshControlHDiv.txt");
            meshcontrol.Print(out);
        }
#endif
        
        std::cout << "MHM Hdiv Computational meshes created\n";
#ifdef PZDEBUG
        if(0)
        {
            std::ofstream gfile("geometryMHMHdiv.txt");
            gmeshauto->Print(gfile);
            std::ofstream out_mhm("MHM_hdiv.txt");
            meshcontrol.CMesh()->Print(out_mhm);
            
        }
#endif
        
        std::cout << "Number of equations MHMixed " << MHMixed->CMesh()->NEquations() << std::endl;

    }
    else
    {
        DebugStop();
    }
    std::string configuration;
    
    {
        std::stringstream sout;
        sout << "H" << Configuration.numHDivisions << "-P" << Configuration.pOrderInternal;
        configuration = sout.str();
    }
    if(Configuration.LagrangeMult)
    {
        MHMPref << "_Lagr";
        MHMMixedPref << "_Lagr";
    }
    if (Configuration.Hybridize) {
        MHMPref << "_Hybr";
        MHMMixedPref << "_Hybr";
    }
    
    // compute the MHM solution
    Configuration.fGlobalSystemWithLocalCondensationSize = MHM->fGlobalSystemWithLocalCondensationSize;
    Configuration.fGlobalSystemSize = MHM->fGlobalSystemSize;
    Configuration.fNumeq = MHM->fNumeq;
    std::cout<< "Begin: Solving and Error post-processing for MHM. " << std::endl;
    SolveProblem(MHM->CMesh(), MHM->GetMeshes(), example, MHMPref.str(), Configuration);
    std::cout<< "End: Solving and Error post-processing for MHM. " << std::endl;
    
    // compute the MHM H(div) solution
    Configuration.fGlobalSystemWithLocalCondensationSize = MHMixed->fGlobalSystemWithLocalCondensationSize;
    Configuration.fGlobalSystemSize = MHMixed->fGlobalSystemSize;
    Configuration.fNumeq = MHMixed->fNumeq;
    std::cout<< "Begin: Solving and Error post-processing for MHM-Hdiv. " << std::endl;
    SolveProblem(MHMixed->CMesh(), MHMixed->GetMeshes(), example, MHMMixedPref.str(), Configuration);
    std::cout<< "End: Solving and Error post-processing for MHM-Hdiv. " << std::endl;
    
//    CopySolution(MHMixed->CMesh().operator->(), MHM->CMesh().operator->());
    
//    if (Configuration.Condensed)
//    {
//        std::string filename = "MHMixed_" + configuration + ".txt";
//        std::ofstream out(filename);
//        PrintElements(MHMixed->CMesh().operator->(), out);
//    }
//    if(Configuration.Condensed)
//    {
//        std::string filename = "MHM_" + configuration + ".txt";
//        std::ofstream out(filename);
//        PrintElements(MHM->CMesh().operator->(), out);
//    }
    
//    ComputeDifferencesBySubmesh(Configuration, MHM, MHMixed, "DiffResults.nb");
//    if(0 && !example)
//    {
//        TPZManVector<STATE,10> square_errors(3,0.);
//        TPZCompMeshTools::ComputeDifferenceNorm(MHMixed->CMesh().operator->(), MHM->CMesh().operator->(), square_errors);
//        std::cout << "Difference between both formulations " << square_errors << std::endl;
//        {
//            std::ofstream out("DiffResults.nb",std::ios::app);
//            out << "(* domain size " << Configuration.nelxcoarse << " " << Configuration.nelycoarse << " num subdomains " << MHM->Coarse_to_Submesh().size() << " *)\n";
//            out << "AppendTo[results, {";
//            out << " ";
//            Configuration.MathematicaInlinePrint(out);
//            out << " ,";
//            out << " {";
//            out << square_errors;
//            out << " } }];\n";
//        }
//    }

    return 0;
}

std::vector <int64_t> ExtractArguments(int argc, char *argv[]){
    
    std::vector <int64_t> parameters;
    if (argc < 6) { // We expect 3 arguments: the program name, the source path and the destination path
        std::cerr << "Usage: " << argv[0] << " <inputs> " << std::endl;
        std::cerr << " <inputs>: n_dx n_dy k_skeleton m_div k_subelement " << std::endl;
        std::cerr << " n_dx number of x division on macro mesh " << std::endl;
        std::cerr << " n_dy number of y division on macro mesh " << std::endl;
        std::cerr << " k_skeleton polynomial order for skeleton variable " << std::endl;
        std::cerr << " m_div number of uniform refinements for subelement mesh " << std::endl;
        std::cerr << " k_subelements polynomial order for subelements variables " << std::endl;
        return parameters;
    }
    
    std::string destination;
    for (int i = 1; i < argc; ++i) { // Remember argv[0] is the path to the program, we want from argv[1] onwards
        parameters.push_back(atoi(argv[i])); // Add all but the last argument to the vector.
    }
    
    return parameters;
}




void InsertMaterialObjects(TPZMHMeshControl &control)
{
    TPZCompMesh &cmesh = control.CMesh();
	/// criar materiais
	int dim = cmesh.Dimension();
    TPZMatLaplacianHybrid *material1 = new TPZMatLaplacianHybrid(matInterno,dim);
    
    material1->SetParameters(1., 0.);
    if(!example)
    {
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Permeability);
        dummy->SetPolynomialOrder(0);
        TPZAutoPointer<TPZFunction<STATE> > func(dummy);
        material1->SetPermeabilityFunction(func);
    }
    else
    {
        material1->SetPermeabilityFunction(example->ConstitutiveLawFunction());
        material1->SetForcingFunction(example->ForcingFunction());
    }
    
    
    
	TPZMaterial * mat1(material1);
    
    TPZMat1dLin *materialCoarse = new TPZMat1dLin(matCoarse);
    TPZFNMatrix<1,STATE> xk(1,1,0.),xb(1,1,0.),xc(1,1,0.),xf(1,1,0.);
    materialCoarse->SetMaterial(xk, xc, xb, xf);    
    cmesh.InsertMaterialObject(materialCoarse);
    
    
    //    REAL diff = -1.;
    //	REAL conv = 0.;
    //	TPZVec<REAL> convdir(3,0.);
    //	REAL flux = 8.;
    //
    //	material1->SetParameters(diff, conv, convdir);
    //	material1->SetInternalFlux( flux);
    
	cmesh.InsertMaterialObject(mat1);
	
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,1.), val2(2,1,0.);
	
    //BC -1
    if (example) {
        TPZBndCond * BCondD1 = dynamic_cast<TPZBndCond *>( material1->CreateBC(mat1, bc1,dirichlet, val1, val2));
        BCondD1->SetType(dirichlet);
        BCondD1->TPZDiscontinuousGalerkin::SetForcingFunction(example->ValueFunction());
        cmesh.InsertMaterialObject(BCondD1);
    }else{
        TPZMaterial * BCondD1 = material1->CreateBC(mat1, bc1,neumann, val1, val2);
//        TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet1 = new TPZDummyFunction<STATE>(DirichletValidacao);
//        BCondD1->SetForcingFunction(bcmatDirichlet1);
        cmesh.InsertMaterialObject(BCondD1);
    }

    
    //BC -2
	TPZMaterial * BCondD2 = material1->CreateBC(mat1, bc2,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet2 = new TPZDummyFunction<STATE>(DirichletValidacao);
    BCondD2->SetForcingFunction(bcmatDirichlet2);
    if (example) {
        BCondD2->SetForcingFunction(example->ValueFunction());
    }
    cmesh.InsertMaterialObject(BCondD2);
    
    //BC -3
    if (example) {
        TPZBndCond* BCondD3 = material1->CreateBC(mat1, bc3,dirichlet, val1, val2);
        BCondD3->SetType(dirichlet);
        BCondD3->TPZDiscontinuousGalerkin::SetForcingFunction(example->ValueFunction());
        cmesh.InsertMaterialObject(BCondD3);
    }
    else{
        TPZMaterial * BCondD3 = material1->CreateBC(mat1, bc3,neumann, val1, val2);
//        TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet3 = new TPZDummyFunction<STATE>(DirichletValidacao);
//        BCondD3->SetForcingFunction(bcmatDirichlet3);
        cmesh.InsertMaterialObject(BCondD3);
    }
    
    //BC -4
	TPZMaterial * BCondD4 = material1->CreateBC(mat1, bc4,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > bcmatDirichlet4 = new TPZDummyFunction<STATE>(DirichletValidacao);
    BCondD4->SetForcingFunction(bcmatDirichlet4);
    if (example) {
        BCondD4->SetForcingFunction(example->ValueFunction());
    }
    cmesh.InsertMaterialObject(BCondD4);
    
    //BC -5: dirichlet nulo
    TPZMaterial * BCondD5 = material1->CreateBC(mat1, bc5,dirichlet, val1, val2);
    cmesh.InsertMaterialObject(BCondD5);
}

void InsertMaterialObjects(TPZMHMixedMeshControl &control)
{
    TPZCompMesh &cmesh = control.CMesh();

    TPZGeoMesh &gmesh = control.GMesh();
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(1,1,0.), val2Flux(1,1,0.), val2Pressure(1,1,10.);
    val2Pressure(0,0) = 1000.;

    int dim = gmesh.Dimension();
    cmesh.SetDimModel(dim);
    
    TPZCompMesh *MixedFluxPressureCmesh = &cmesh;
    
    // Material medio poroso
    TPZMixedPoisson * mat = new TPZMixedPoisson(1,dim);
    mat->SetSymmetric();
    mat->SetPermeability(1.);
    if(!example)
    {
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Permeability);
        dummy->SetPolynomialOrder(0);
        TPZAutoPointer<TPZFunction<STATE> > func(dummy);
        mat->SetPermeabilityFunction(func);
    } else
    {
        mat->SetPermeabilityFunction(example->ConstitutiveLawFunction());
        mat->SetForcingFunction(example->ForcingFunction());
    }
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat);
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, bc1, typeFlux, val1, val2Flux);
    TPZAutoPointer<TPZFunction<STATE> > force = new TPZDummyFunction<STATE>(DirichletValidacao);
//    bcN->SetForcingFunction(0, force);
    if (example) {
        bcN->SetType(typePressure);
        bcN->TPZMaterial::SetForcingFunction(example->ValueFunction());
    }
    
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    bcN = mat->CreateBC(mat, bc3, typeFlux, val1, val2Flux);
//    bcN->SetForcingFunction(0, force);
    if (example) {
        bcN->SetType(typePressure);
        bcN->TPZDiscontinuousGalerkin::SetForcingFunction(example->ValueFunction());
    }
    
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, bc2, typePressure, val1, val2Pressure);
    bcS->SetForcingFunction(0, force);
    if (example) {
        bcS->TPZMaterial::SetForcingFunction(example->ValueFunction());
    }

    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    bcS = mat->CreateBC(mat, bc4, typePressure, val1, val2Pressure);
    bcS->SetForcingFunction(0, force);
    if (example) {
        bcS->TPZMaterial::SetForcingFunction(example->ValueFunction());
    }
    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    
//    TPZMatLaplacian *matdim = new TPZMatLaplacian(1);
//    matdim->SetDimension(gmesh.Dimension());
//    control.PressureMesh()->InsertMaterialObject(matdim);

    
//    control.InsertPeriferalMaterialObjects();
    

}

/// Insert material objects for the MHM-H(div) solution
void InsertMaterialObjects(TPZMHMixedHybridMeshControl &control)
{
    TPZGeoMesh &gmesh = control.GMesh();
    
    /// Tem que ver isso melhor : nao eh qualquer malha MHM que quer este material
    if(0)
    {
        int meshdim = gmesh.Dimension();
        control.fFractureFlowDim1MatId.insert(10);
        // Material medio poroso
        TPZMixedPoisson * mat = new TPZMixedPoisson(10,meshdim-1);
        mat->SetSymmetric();
        mat->SetPermeability(1.e-3);
        TPZFNMatrix<9,REAL> K(3,3,0.),KInv(3,3,0.);
        K(0,0) = 1.;
        K(1,1) = 1.e-3;
        KInv(0,0) = 1.;
        KInv(1,1) = 1000.;
        mat->SetPermeabilityTensor(K, KInv);
        control.CMesh()->InsertMaterialObject(mat);
    }
    
    
    InsertMaterialObjects((TPZMHMixedMeshControl &) control);
    
}



TPZCompMesh * CreateHDivMHMMesh(TPZGeoMesh * gmesh, int porder)
{
    int meshdim = gmesh->Dimension();
    TPZCompMesh * cmeshHDiv = new TPZCompMesh(gmesh);
    cmeshHDiv->SetDimModel(meshdim);
    cmeshHDiv->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    cmeshHDiv->SetDefaultOrder(porder);
    TPZVecL2 *matl2 = new TPZVecL2(1);
    matl2->SetDimension(2);
    cmeshHDiv->InsertMaterialObject(matl2);
    for (int matid = 2; matid<10; matid++) {
        TPZVecL2 *matl2 = new TPZVecL2(matid);
        matl2->SetDimension(2);
        cmeshHDiv->InsertMaterialObject(matl2);
    }
    TPZFNMatrix<1,STATE> val1(1,1,0.),val2(1,1,0.);
    TPZBndCond *bc = matl2->CreateBC(matl2, -1, 0, val1, val2);
    cmeshHDiv->InsertMaterialObject(bc);
    bc = matl2->CreateBC(matl2, -2, 0, val1, val2);
    cmeshHDiv->InsertMaterialObject(bc);
    cmeshHDiv->AutoBuild();
    
#ifdef PZDEBUG
    {
        std::ofstream outmesh("BigHDivMesh.txt");
        cmeshHDiv->Print(outmesh);
    }
#endif
    return cmeshHDiv;
}


TPZCompMesh * CreatePressureMHMMesh(TPZGeoMesh * gmesh, int porder, int dimension)
{
    TPZCompMesh * cmeshPressure = new TPZCompMesh(gmesh);
    cmeshPressure->SetDimModel(dimension);
    cmeshPressure->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmeshPressure->ApproxSpace().CreateDisconnectedElements(true);
    cmeshPressure->SetDefaultOrder(porder);
    TPZMatLaplacian *matl2 = new TPZMatLaplacian(1);
    matl2->SetDimension(dimension);
    cmeshPressure->InsertMaterialObject(matl2);
    for (int matid = 2; matid<10; matid++) {
        TPZMatLaplacian *matl2 = new TPZMatLaplacian(matid);
        matl2->SetDimension(dimension);
        cmeshPressure->InsertMaterialObject(matl2);
    }

    cmeshPressure->AutoBuild();
    int64_t nc = cmeshPressure->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        cmeshPressure->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    return cmeshPressure;
}

TPZCompMesh * CreateHDivPressureMHMMesh(TPZVec<TPZCompMesh * > & cmeshes)
{
    TPZGeoMesh *gmesh = cmeshes[0]->Reference();
    if(!gmesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    int dim = gmesh->Dimension();
    
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(1,1,0.), val2Flux(1,1,0.), val2Pressure(1,1,10.);
    val2Pressure(0,0) = 1000.;
    
    // Malha computacional
    TPZCompMesh * MixedFluxPressureCmesh = new TPZCompMesh(gmesh);
    
    // Material medio poroso
    TPZMixedPoisson * mat = new TPZMixedPoisson(1,dim);
    mat->SetSymmetric();
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat);

    for (int matid = 2; matid<10; matid++)
    {
        // Material medio poroso
        TPZMixedPoisson * mat = new TPZMixedPoisson(matid,dim);
        mat->SetSymmetric();
        //    mat->SetForcingFunction(One);
        MixedFluxPressureCmesh->InsertMaterialObject(mat);
        

    }
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, -1, typeFlux, val1, val2Flux);
    TPZAutoPointer<TPZFunction<STATE> > force = new TPZDummyFunction<STATE>(DirichletValidacao);
    //    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, -2, typePressure, val1, val2Pressure);
    bcS->SetForcingFunction(0, force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    
    
    
    
    
    MixedFluxPressureCmesh->SetDimModel(dim);
    MixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElem();
    MixedFluxPressureCmesh->AutoBuild();
    
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    
    
    meshvector[0] = cmeshes[0];
    meshvector[1] = cmeshes[1];
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, MixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, MixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, MixedFluxPressureCmesh);
    
    return MixedFluxPressureCmesh;

}


TPZAutoPointer<TPZRefPattern> DivideQuadbyTriangles(const std::string refpatname)
{
    TPZGeoMesh gmesh;
    gmesh.NodeVec().Resize(5);
    REAL nodeco[][3] =
    {
        {-1,-1,0},
        {1,-1,0},
        {1,1,0},
        {-1,1,0},
        {0,0,0}
    };
    int64_t nodeindexes[][3] = {
        {0,1,4},
        {1,2,4},
        {2,3,4},
        {3,0,4}
    };
    for (int i=0; i<5; i++) {
        TPZManVector<REAL,3> coord(3);
        for (int c=0; c<3; c++) {
            coord[c] = nodeco[i][c];
        }
        gmesh.NodeVec()[i].Initialize(coord, gmesh);
    }
    TPZManVector<int64_t> corners(4);
    for (int64_t i=0; i<4; i++) {
        corners[i] = i;
    }
    int64_t elindex;
    gmesh.CreateGeoElement(EQuadrilateral, corners, 1, elindex);
    
    int64_t fatherindex = elindex;
    
    for (int is=0; is<4; is++)
    {
        for (int64_t i=0; i<3; i++) {
            corners[i] = nodeindexes[is][i];
        }
        gmesh.CreateGeoElement(ETriangle, corners, 1, elindex);
        gmesh.Element(elindex)->SetFather(fatherindex);
    }
    gmesh.BuildConnectivity();
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(gmesh);
    refpat->SetName(refpatname);
    TPZAutoPointer<TPZRefPattern> found = gRefDBase.FindRefPattern(refpat);
    if(!found)
    {
        gRefDBase.InsertRefPattern(refpat);
        refpat->InsertPermuted();
    }
    else
    {
        refpat = found;
    }
    return refpat;
}


TPZAutoPointer<TPZRefPattern> DivideTriangleby9Triangles(const std::string refpatname)
{
    TPZGeoMesh gmesh;
    gmesh.NodeVec().Resize(10);
    REAL nodeco[][3] =
    {
        {0,0,0}, //0
        {1,0,0}, //1
        {2,0,0},  //2
        {3,0,0},  //3
        {0,1,0},  //4
        {1,1,0},  //5
        {2,1,0},  //6
        {0,2,0},  //7
        {1,2,0},  //8
        {0,3,0} //9
    };
    int64_t nodeindexes[][3] = {
        {0,3,9},
        {0,1,4},
        {1,5,4},
        {1,2,5},
        {2,6,5},
        {2,3,6},
        {4,5,7},
        {5,8,7},
        {5,6,8},
        {7,8,9}
    };
    for (int i=0; i<10; i++) {
        TPZManVector<REAL,3> coord(3);
        for (int c=0; c<3; c++) {
            coord[c] = nodeco[i][c];
        }
        gmesh.NodeVec()[i].Initialize(coord, gmesh);
    }
    TPZManVector<int64_t> corners(3);
    for (int64_t i=0; i<3; i++) {
        corners[i] = nodeindexes[0][i];
    }
    int64_t elindex;
    gmesh.CreateGeoElement(ETriangle, corners, 1, elindex);
    
    int64_t fatherindex = elindex;
    
    for (int is=1; is<10; is++)
    {
        for (int64_t i=0; i<3; i++) {
            corners[i] = nodeindexes[is][i];
        }
        gmesh.CreateGeoElement(ETriangle, corners, 1, elindex);
        gmesh.Element(elindex)->SetFather(fatherindex);
    }
    gmesh.BuildConnectivity();
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(gmesh);
    refpat->SetName(refpatname);
    if(!gRefDBase.FindRefPattern(refpat))
    {
        gRefDBase.InsertRefPattern(refpat);
        refpat->InsertPermuted();
    }
    return refpat;

}

TPZGeoMesh *MalhaGeomFred(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, const std::string quad, const std::string triangle, TPZVec<int64_t> &coarseindices, int ndiv)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    int dimension = 2;
    gmesh->SetDimension(dimension);
    TPZManVector<int,2> nx(2,3);
    nx[0] = nelx;
    nx[1] = nely;
    TPZGenGrid gengrid(nx, x0, x1);
    gengrid.SetRefpatternElements(true);
    gengrid.Read(gmesh, 1);
    gengrid.SetBC(gmesh, 4, -1);
    gengrid.SetBC(gmesh, 5, -2);
    gengrid.SetBC(gmesh, 6, -1);
    gengrid.SetBC(gmesh, 7, -2);
    
    TPZAutoPointer<TPZRefPattern> refquad,reftriangle;
    refquad = gRefDBase.FindRefPattern(quad);
    reftriangle = gRefDBase.FindRefPattern(triangle);
    if (!refquad || ! reftriangle) {
        DebugStop();
    }
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->Type() == EQuadrilateral   ) {
            gel->SetRefPattern(refquad);
            TPZManVector<TPZGeoEl *,4> subs;
            gel->Divide(subs);
        }
    }
    nel = gmesh->NElements();

    coarseindices.resize(nel);
    int64_t elcount = 0;
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement() ||  gel->Dimension() != dimension) {
            continue;
        }
        coarseindices[elcount] = el;
        elcount++;
    }
    coarseindices.resize(elcount);
#ifdef PZDEBUG
    {
        std::ofstream file("GMeshFredCoarse.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }
#endif

    if(1)
    {
    
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel->HasSubElement() &&  gel->Type() == ETriangle) {
                gel->SetRefPattern(reftriangle);
                TPZManVector<TPZGeoEl *,12> subs;
                gel->Divide(subs);
            }
        }
        nel = gmesh->NElements();
        
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel->HasSubElement() &&  gel->Type() == EOned) {
                TPZAutoPointer<TPZRefPattern> refpat = TPZRefPatternTools::PerfectMatchRefPattern(gel);
                if (!refpat) {
                    DebugStop();
                }
                gel->SetRefPattern(refpat);
                TPZManVector<TPZGeoEl *,12> subs;
                gel->Divide(subs);
            }
        }
    }
    
#ifdef PZDEBUG
    {
        std::ofstream file("GMeshFredInit.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }
#endif

    TPZCheckGeom geom(gmesh);
    geom.UniformRefine(ndiv);
//    InsertInterfaceElements(gmesh,1,2);

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
#ifdef PZDEBUG
    {
        std::ofstream file("GMeshFred.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }
#endif
    
    return gmesh;
}


/// Create a reference geometric mesh starting with nelx by nely domains
TPZGeoMesh *CreateReferenceGMesh(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, int numref)
{
    TPZManVector<int,3> nx(2);
    nx[0] = nelx;
    nx[1] = nely;
    TPZGenGrid gengrid(nx,x0,x1);
    gengrid.SetRefpatternElements(true);
    TPZGeoMesh *result = new TPZGeoMesh;
    int matid = 1;
    gengrid.Read(result, matid);
    gengrid.SetBC(result, 4, -1);
    gengrid.SetBC(result, 5, -2);
    gengrid.SetBC(result, 6, -1);
    gengrid.SetBC(result, 7, -2);
    int matidpoint = 10;
    
    int64_t firstnode = nelx+2;
    int64_t numnodes = nelx-1;
    for (int64_t ynode = 1; ynode < nely; ynode++)
    {
        for (int64_t node = firstnode; node < firstnode+numnodes; node++) {
            TPZManVector<int64_t,2> nodeindices(1);
            nodeindices[0] = node;
            int64_t index;
            result->CreateGeoElement(EPoint, nodeindices, matidpoint, index);
        }
        firstnode += nelx+1;
    }
    result->BuildConnectivity();
    // refina a malha uma vez uniformemente
    int numuni = 1;
    for (int uni=0; uni<numuni; uni++)
    {
        int64_t nelem = result->NElements();
        for (int64_t el=0; el<nelem; el++) {
            TPZGeoEl *gel = result->Element(el);
            if (gel->Dimension() == 0) {
                continue;
            }
            TPZManVector<TPZGeoEl *,8> subs;
            gel->Divide(subs);
        }
    }
    // refina a malha na direcao dos elementos ponto
    std::set<int> matids;
    matids.insert(matidpoint);
    for (int cycle = 0; cycle < numref; cycle++) {
        int64_t nelem = result->NElements();
        for (int64_t el=0; el<nelem; el++) {
            TPZGeoEl *gel = result->Element(el);
            int targetmatid = cycle+2;
            TPZRefPatternTools::RefineDirectional(gel, matids, targetmatid);
        }
    }
    
    {
        std::ofstream out("../ReferenceGMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(result, out);
    }
    return result;
}

void Permeability(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &diff)
{
    int factor = 500;
    int64_t ix = x[0]*factor;
    int64_t iy = x[1]*factor;
    static int count = 0;
    if((fabs(ix-x[0]*factor) < 1.e-6 || fabs(ix-x[1]*factor) < 1.e-6) && count < 20)
    {
        count++;
#ifdef PZDEBUG
        std::cout << "probing for a permeability at the interface of two regions\n";
        std::cout << "x = " << x << std::endl;
#endif
    }
    if (IsZero(x[1]-0.2)) {
        iy = 99;
    }
    if (IsZero(x[0]-1.)) {
        ix = 499;
    }
    
//    std::cout << "ix = " << ix << std::endl;
//    std::cout << "iy = " << iy << std::endl;
    
    REAL valporous = gPorous(ix,iy);
    // totototototo
    //    valporous = 1.+0.3*sin(x[0]*50)*cos(x[1]*50.);
    for (int i=0; i<2; i++) {
        diff(i,i) = valporous;
        diff(i,1-i)=0.;
        diff(2+i,i) = 1./valporous;
        diff(2+i,1-i) = 0.;
    }
//    for (int i=0; i<2; i++) {
//        diff(i,i) = 1.;
//        diff(2+i,i) = 1.;
//    }
}

/// compute the reference solution and return created mesh
TPZCompMesh *ComputeReferenceSolution(TPZGeoMesh *gmesh, int porder, TPZVec<TPZCompMesh *> &meshvec)
{
    TPZCompMesh *CHDivPressureMesh = CreateReferenceCMesh(gmesh, meshvec, porder);
    //calculo solution
    TPZAnalysis an(CHDivPressureMesh);
    std::cout << "Assembling and Solving " << CHDivPressureMesh->NEquations() << " equations\n";
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(CHDivPressureMesh);
    strmat.SetNumThreads(0);
    an.SetStructuralMatrix(strmat);

#else
    TPZSkylineStructMatrix strmat(CHDivPressureMesh);
#endif
#ifndef PZDEBUG
    //    skyl.SetNumThreads(16);
#endif
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();
    if(0)
    {
        std::ofstream global("Global.nb");
        TPZAutoPointer<TPZStructMatrix> strmat = an.StructMatrix();
        an.Solver().Matrix()->Print("Glob = ",global,EMathematicaInput);
        an.Rhs().Print("Rhs = ",global,EMathematicaInput);
    }
    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
#ifdef PZDEBUG
    {
        std::ofstream out("MeshWithSol.txt");
        CHDivPressureMesh->Print(out);
    }
#endif
    an.LoadSolution(); // compute internal dofs
                       //    an.Solution().Print("sol = ");
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, an.Mesh());
    //    TPZBuildMultiphysicsMesh::TransferFromMeshes(cmeshes, an.Mesh());
    //    for (int i=0; i<cmeshes.size(); i++) {
    //        cmeshes[i]->Solution().Print("sol = ");
    //    }
    //    cmeshes[0]->Solution().Print("solq = ");
    //    cmeshes[1]->Solution().Print("solp = ");
    UnwrapMesh(an.Mesh());
    std::string plotfile("referencesolution.vtk");
    TPZStack<std::string> scalnames,vecnames;
    scalnames.Push("Pressure");
    scalnames.Push("Permeability");
    vecnames.Push("Derivative");
    vecnames.Push("Flux");
    an.DefineGraphMesh(CHDivPressureMesh->Dimension(), scalnames, vecnames, plotfile);
    int resolution = 0;
    an.PostProcess(resolution,CHDivPressureMesh->Dimension());

    return CHDivPressureMesh;
}

/// create the computational mesh of the reference solution
TPZCompMesh *CreateReferenceCMesh(TPZGeoMesh *gmesh, TPZVec<TPZCompMesh *> &meshvec, int porder)
{
    meshvec.resize(2);
    int dimension = 2;
    meshvec[0] = CreateHDivMHMMesh(gmesh,porder);
    meshvec[1] = CreatePressureMHMMesh(gmesh, porder, dimension);
    int hdivplusplus=1;
    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshvec[0], hdivplusplus);
    TPZCompMeshTools::SetPressureOrders(meshvec[0], meshvec[1]);
    TPZCompMesh *cmesh = CreateHDivPressureMHMMesh(meshvec);
    for (int i=0; i<10; i++) {
        TPZMaterial *mat = cmesh->FindMaterial(i);
        if (mat) {
            TPZMixedPoisson *mixed = dynamic_cast<TPZMixedPoisson *>(mat);
            if (!mixed) {
                DebugStop();
            }
            TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Permeability);
            dummy->SetPolynomialOrder(0);
            TPZAutoPointer<TPZFunction<STATE> > func(dummy);
            mixed->SetPermeabilityFunction(func);
        }
    }
    TPZCompMeshTools::GroupElements(cmesh);
    TPZCompMeshTools::CreatedCondensedElements(cmesh, true);
    return cmesh;
}

void UnwrapMesh(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    bool change = true;
    while(change)
    {
        change = false;
        for (int64_t el=0; el<nel; el++) {
            
            TPZCompEl *cel = cmesh->Element(el);
            TPZCondensedCompEl *condense = dynamic_cast<TPZCondensedCompEl *>(cel);
            if (condense) {
                condense->Unwrap();
                change = true;
            }
            cel = cmesh->Element(el);
            TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
            if (elgr) {
                elgr->Unwrap();
                change = true;
            }
        }
    }
}

REAL objectivefunc(REAL K1, REAL K2, REAL K3, REAL K4, REAL Lambda)
{
    if (abs(Lambda-1.) < 1.e-6)
    {
        REAL val = 8*K2*K3*(K1*K1*(K3*(-K3 + K4) + K2*(K3 + K4)) +
                             K2*K4*(K2*(K3 - K4) + K3*(K3 + K4)) +
                             K1*(K2*K2*(K3 + K4) + K3*K4*(K3 + K4) +
                                 K2*(K3*K3 + 6*K3*K4 + K4*K4)) +
                             (K1 + K2)*(K2 + K3)*(K1 + K4)*(K3 + K4)*cos(M_PI*Lambda));
        return val;
    }
    REAL val = (8*K2*K3*(K1*K1*(K3*(-K3 + K4) + K2*(K3 + K4)) +
              K2*K4*(K2*(K3 - K4) + K3*(K3 + K4)) +
              K1*(K2*K2*(K3 + K4) + K3*K4*(K3 + K4) +
                  K2*(K3*K3 + 6*K3*K4 + K4*K4)) +
              (K1 + K2)*(K2 + K3)*(K1 + K4)*(K3 + K4)*cos(M_PI*Lambda))*tan((M_PI*Lambda)/2.))/
    (K4*(K1*(K2*(K3 - K4) - K3*(K3 + K4)) + K2*(K2*(K3 - K4) + K3*(K3 + K4)) +
         (K1 + K2)*(K2 + K3)*(K3 + K4)*cos(M_PI*Lambda)));
    return val;
}

REAL Power(REAL val, int expon)
{
    if (expon != 2) {
        DebugStop();
    }
    return val*val;
}

REAL Sec(REAL val)
{
    return 1./cos(val);
}

REAL Dobjectivefunc(REAL K1, REAL K2, REAL K3, REAL K4, REAL Lambda)
{
    const double Pi = M_PI;
    REAL nom = 4*K2*K3*Pi*(2*Power(K1 + K2,2)*Power(K2 + K3,2)*(K1 + K4)*Power(K3 + K4,2)*cos(Pi*Lambda) +
                           2*(K1 + K2)*(K2 + K3)*(K3 + K4)*
                           (K1*K3*(K1*(K2 - 3*K3) + K2*(K2 + K3)) +
                            (Power(K1,2)*(K2 + K3) + K2*K3*(K2 + K3) +
                             K1*(Power(K2,2) + 10*K2*K3 + Power(K3,2)))*K4 +
                            (K2*(-3*K2 + K3) + K1*(K2 + K3))*Power(K4,2) -
                            2*K1*(K2 + K3)*K4*(K1 + K2 + K3 + K4)*cos(Pi*Lambda)) +
                           4*Power(K1*K3 - K2*K4,2)*(Power(K2,2)*K4 + K1*(Power(K3,2) + (K2 + K3)*K4))*
                           Power(Sec((Pi*Lambda)/2.),2));
//    4*K2*K3*M_PI*(
//    2*(K1 + K2)*(K1+K2)*(K2 + K3)*(K2+K3)*(K1 + K4)*(K3 + K4)*(K3+K4)*cos(M_PI*Lambda) +
//                            2*(K1 + K2)*(K2 + K3)*(K3 + K4)*(
//                                K1*K3*(K1*(K2 - 3*K3) + K2*(K2 + K3)) +
//                                    ((K1*K1)*(K2 + K3) + K2*K3*(K2 + K3) + K1*((K2*K2) + 10*K2*K3 + (K3*K3)))*K4 +
//                             (K2*(-3*K2 + K3) + K1*(K2 + K3))*(K4*K4) -
//                             2*K1*(K2 + K3)*K4*(K1 + K2 + K3 + K4)*cos(M_PI*Lambda)
//                                                             ) +
//                            4*(K1*K3 - K2*K4)*(K1*K3 - K2*K4)*((K2*K2)*K4 + K1*((K3*K3) + (K2 + K3)*K4))/((cos((M_PI*Lambda)/2.)*cos((M_PI*Lambda)/2.)))
//                             );
    REAL denom =
    (K4*
     (K1*(K2 - K3)*K3 - K1*(K2 + K3)*K4 + K2*(K2*(K3 - K4) + K3*(K3 + K4)) + (K1 + K2)*(K2 + K3)*(K3 + K4)*cos(M_PI*Lambda))
     *(K1*(K2 - K3)*K3 - K1*(K2 + K3)*K4 + K2*(K2*(K3 - K4) + K3*(K3 + K4)) + (K1 + K2)*(K2 + K3)*(K3 + K4)*cos(M_PI*Lambda))
     );
    REAL val = nom/denom;
    return val;
}

REAL StartGuess(REAL K1, REAL K2, REAL K3, REAL K4)
{
    REAL inc = 0.1;
    REAL start = 0.60;
    REAL val = objectivefunc(K1,K2,K3,K4,start);
    if(val > 0.) DebugStop();
    while(abs(inc) > 1.e-2)
    {
        REAL val2 = objectivefunc(K1,K2,K3,K4,start+inc);
        if (val2 < 0) {
            start += inc;
            if (start > 1.-1.e-6) {
                start -=inc;
                inc /=2.;
            }
        }
        else
        {
            inc /= 2.;
        }
    }
    return start;
}

REAL ConvergeNewton(REAL K1, REAL K2, REAL K3, REAL K4, REAL start)
{
    REAL val = objectivefunc(K1, K2, K3, K4, start);
    int nstep = 0;
    while(abs(val) > 1.e-6 && nstep < 15)
    {
        REAL deriv = Dobjectivefunc(K1, K2, K3, K4, start);
        start -= val/deriv;
        val = objectivefunc(K1, K2, K3, K4, start);
        nstep++;
    }
    if (nstep == 15) {
        std::cout << "val = " << val << " lambda = " << start << " nsteps " << nstep;
    }
    return start;
}

REAL SolutionRegularity(TPZFMatrix<REAL> &perms)
{
    REAL K1,K2,K3,K4;
    K1 = perms(0,0);
    K2 = perms(1,0);
    K3 = perms(1,1);
    K4 = perms(0,1);
    REAL lambda = StartGuess(K1, K2, K3, K4);
    lambda = ConvergeNewton(K1, K2, K3, K4, lambda);
    return lambda;
}

/// Analise the regularity of the subdomain problems
void AnalyseRegularity(const TPZVec<int> &pos0,const TPZVec<int> &nelx, TPZVec<int> &nsub, TPZFMatrix<REAL> &lowestexp)
{
#ifdef PZDEBUG
    if(nelx[0]%nsub[0] || nelx[1]%nsub[1])
    {
        DebugStop();
    }
#endif
    lowestexp.Redim(nelx[0]-1, nelx[1]-1);
    for (int ir=0; ir<nelx[0]-1; ir++) {
        for (int ic=0; ic<nelx[1]-1; ic++) {
            TPZFNMatrix<4,REAL> perm(2,2);
            perm(0,0) = gPorous(pos0[0]+ir,pos0[1]+ic);
            perm(1,0) = gPorous(pos0[0]+ir+1,pos0[1]+ic);
            perm(0,1) = gPorous(pos0[0]+ir,pos0[1]+ic+1);
            perm(1,1) = gPorous(pos0[0]+ir+1,pos0[1]+ic+1);
            lowestexp(ir,ic) = SolutionRegularity(perm);
            if (lowestexp(ir,ic) < 0 || lowestexp(ir,ic) > 1. || isnan(lowestexp(ir,ic))) {
                DebugStop();
            }
        }
    }
}

/// Print the elements with geometric information and connect values
void PrintElements(TPZCompMesh *cmesh, std::ostream &out)
{
    int64_t nelem = cmesh->NElements();
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        if (gel->Dimension() == cmesh->Dimension()) {
            DebugStop();
        }
        TPZManVector<REAL,3> co1(3),co2(3,-100.);
        gel->Node(0).GetCoordinates(co1);
        if(gel->NCornerNodes() > 1) gel->Node(1).GetCoordinates(co2);
        out << "gel index " << gel->Index() << " node loc " << co1 << " and " << co2 << std::endl;
        int nc = cel->NConnects();
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if (c.NShape()) {
                c.Print(*cmesh,out);
            }
        }
    }
}

/// copy the solution between one computation mesh to the other assuming the geometric elements match
void CopySolution(TPZCompMesh *from, TPZCompMesh *to)
{
    int64_t nelem = from->NElements();
    TPZGeoMesh *gto = to->Reference();
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *celfrom = from->Element(el);
        if(!celfrom) continue;
        TPZGeoEl *gelfrom = celfrom->Reference();
        if(!gelfrom) continue;
        if (gelfrom->Dimension() != 1) {
            DebugStop();
        }
        int64_t index = gelfrom->Index();
        TPZGeoEl *gelto = gto->Element(index);
        TPZCompEl *celto = gelto->Reference();
        int nc = celfrom->NConnects();
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &cfrom = celfrom->Connect(ic);
            TPZConnect &cto = celto->Connect(ic);
            int64_t seqfrom = cfrom.SequenceNumber();
            int64_t seqto = cto.SequenceNumber();
            int64_t size = from->Block().Size(seqfrom);
            for (int i=0; i<size; i++) {
                to->Block()(seqto,0,i,0) = from->Block()(seqfrom,0,i,0);
            }
        }
    }
    to->LoadSolution(to->Solution());
}

/// Compute the differences at the submesh level
void ComputeDifferencesBySubmesh(TRunConfig &config, TPZMHMeshControl &MHM, TPZMHMixedMeshControl &MHMixed, const std::string &filename)
{
    int64_t nsubmeshes = MHM.Coarse_to_Submesh().size();
    TPZManVector<STATE> difference(nsubmeshes,0.);
    TPZManVector<STATE,10> square_errors(3,0.);
    std::map<int64_t,int64_t>::iterator it;
    int64_t count = 0;
    for (it = MHM.Coarse_to_Submesh().begin(); it != MHM.Coarse_to_Submesh().end(); it++)
    {
        square_errors[0] = 0;
        square_errors[1] = 0;
        square_errors[2] = 0;
        int64_t skelindex = it->first;
        int64_t MHM_index = it->second;
        if (MHMixed.Coarse_to_Submesh().find(skelindex) == MHMixed.Coarse_to_Submesh().end()) {
            DebugStop();
        }
        int64_t MHMixed_index = MHMixed.Coarse_to_Submesh()[skelindex];
        if (MHM_index == -1 || MHMixed_index == -1) {
            continue;
        }
        TPZSubCompMesh *submhm = dynamic_cast<TPZSubCompMesh *>(MHM.CMesh()->Element(MHM_index));
        TPZSubCompMesh *submhmixed = dynamic_cast<TPZSubCompMesh *>(MHMixed.CMesh()->Element(MHMixed_index));
        TPZCompMeshTools::ComputeDifferenceNorm(submhmixed, submhm, square_errors);
        difference[count] = square_errors[0];
        //        count++;
        //        TPZCompMeshTools::ComputeDifferenceNorm(submhm, submhmixed, square_errors);
        //        difference[count] = square_errors[0];
        count++;
    }
    //    ComputationalMesh->Reference()->ResetReference();
    //    ComputationalMesh->LoadReferences();
    //    return 0;
    {
        std::ofstream out("DiffResults.nb",std::ios::app);
        out << "(* domain size " << config.nelxcoarse << " " << config.nelycoarse << " num subdomains " << count << " *)\n";
        out << "AppendTo[results, {";
        out << " "  << config.numHDivisions << " , " << config.pOrderInternal << " ,";
        out << " {";
        for (int64_t el=0; el<difference.size(); el++) {
            out << difference[el];
            if (el != difference.size()-1) {
                out << " , ";
            }
        }
        out << " } }];\n";
    }
    
}

/// function that randomly refines some elements
void RandomRefine(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices, int nref)
{
    int64_t nel = coarseindices.size();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(coarseindices[el]);
        while (gel->HasSubElement()) {
            int nsub = gel->NSubElements();
            int isub = rand()%nsub;
            gel = gel->SubElement(isub);
        }
        for (int iref = 0; iref<nref; iref++)
        {
            TPZManVector<TPZGeoEl *,10> gelsub;
            gel->Divide(gelsub);
            int nsub = gel->NSubElements();
            int isub = rand()%nsub;
            gel = gel->SubElement(isub);
        }
    }
}
