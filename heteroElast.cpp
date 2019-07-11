
#include "pzlog.h"
#include "pzgmesh.h"
#include "pzmanvector.h"
#include "TPZMHMeshControl.h"
#include "TPZVTKGeoMesh.h"

#include "pzelasmat.h"
#include "TPZElasticity2DHybrid.h"
#include "pzmat1dlin.h"
#include "TPZNullMaterial.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "TPZAnalyticSolution.h"
#include "meshgen.h"
#include "pzgengrid.h"
#include "pzcheckgeom.h"

#ifndef USING_MKL
#include "pzskylstrmatrix.h"
#endif

/// Insert material objects for the MHM Mesh solution
void InsertMaterialObjects(TPZCompMesh *cmesh);

//int ReadFromFile(TPZFMatrix<double> &mat, string path);
//int ReadFromFile2(TPZFMatrix<double> &mat, string path);

//std::pair<double,double> funcE(TPZFMatrix<double> &ElastCoef, TPZFMatrix<double> &PoissonCoef, double x,double y,double min_x, double max_x,double min_y, double max_y, int nx, int ny);



/*void funcE2(const TPZVec<REAL> &x, TPZVec<STATE> &func, TPZFMatrix<STATE> &deriv)
{
    std::pair<double,double> val;
    val = funcE(ElastCoef, PoissonCoef, x[0], x[1], min_x, max_x, min_y, max_y, nx, ny);
    func[0] = val.first;
    func[1] = val.second;
}*/


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
const int neumann = 1;

int const bc1=-1;
int const bc2=-2;
int const bc3=-3;
int const bc4=-4;
int const bc5=-5;


//TPZFMatrix<REAL> ElastCoef;
//REAL **ElastCoef;




TPZAnalyticSolution *example;
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif


    int flag = 0;
    string path1, path2;

#ifdef MACOSX
    path1 = "../young_modulus.txt";
    path2 = "../Poisson_coef.txt";
#else
    path1 = "young_modulus.txt";
    path2 = "Poisson_coef.txt";
#endif

    flag = ReadFromFile( ElastCoef, path1);
    if (!flag)
    {
        DebugStop();
    }
    flag = ReadFromFile( PoissonCoef, path2);
    if (!flag)
    {
        DebugStop();
    }




    max_x = 10000;
    max_y = 4500;
    min_x = 0.0;
    min_y = 0.0;


/*
    if (0) {
        std::pair<double, double> numE = funcE(ElastCoef, PoissonCoef, xx, yy, min_x, max_x, min_y, max_y, nx, ny);

        TPZManVector<double> x_coord(2), result(2);
        TPZFMatrix<double> grad;
        x_coord[0] = xx;
        x_coord[1] = yy;
        for (int i = 1; i < 5; i++) {
            for (int j = 1; j < 5; j++) {
                x_coord[0] = max_x * i / 4;
                x_coord[1] = max_y * j / 4;
                funcE2(x_coord, result, grad);
                std::cout << "x_coord = " << x_coord << " Val result = " << result << std::endl;
            }
        }
        funcE2(x_coord, result, grad);


        std::cout << "Val numE[0] = " << numE.first << "Val numE[1] = " << numE.second << " Val result = " << result
                  << std::endl;
    }
*/



    HDivPiola = 1;

    TPZVec<int> nx(2,5);
    nx[0] =2;
    nx[1] =1;
    TPZVec<REAL> x0(3,0.), x1(3,0.);
    x1[0]=max_x;
    x1[1]=max_y;
    x1[2]=0.0;

    TPZGenGrid grid(nx,x0,x1);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
 //   grid.SetElementType(ETriangle);
    grid.Read(gmesh);
    grid.SetBC(gmesh,4,-1);
    grid.SetBC(gmesh,5,-2);
    grid.SetBC(gmesh,6,-3);
    grid.SetBC(gmesh,7,-4);

    TPZCheckGeom check(gmesh);
    check.UniformRefine(7);

    {
        std::ofstream out("gmesh.colsvtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,out,true);

    }


    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    InsertMaterialObjects(cmesh);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(3);
    cmesh->AutoBuild();
    {
        std::ofstream out("cmesh.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh,out,true);

    }

    TPZAnalysis an(cmesh);

    TPZSymetricSpStructMatrix strmat(cmesh);
    strmat.SetNumThreads(4);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> solver;
    solver.SetDirect(ECholesky);
    an.SetSolver(solver);

    an.Run();

    TPZStack<std::string> scalnames,vecnames;
    scalnames.Push("J2");
    scalnames.Push("I1");
    scalnames.Push("SigmaX");
    scalnames.Push("SigmaY");
    scalnames.Push("TauXY");
    scalnames.Push("Young_Modulus");
    scalnames.Push("Poisson");
    vecnames.Push("displacement");
    an.DefineGraphMesh(2,scalnames,vecnames,"saida.vtk");
    an.PostProcess(0);

    return 0;
}





void InsertMaterialObjects(TPZCompMesh *cmesh)
{
    /// criar materiais
//    int dim = cmesh.Dimension();
    int matInterno = 1;
    STATE Young = 1000., nu = 0.3, fx = 0., fy = 0.;
    TPZElasticityMaterial *material1 = new TPZElasticityMaterial(matInterno,Young,nu,fx,fy,0);
    material1->SetPlaneStrain();

    TPZAutoPointer<TPZFunction<STATE>> func = new TPZDummyFunction<STATE> (funcE2, 1);
    material1->SetElasticityFunction(func);

    cmesh->InsertMaterialObject(material1);


    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    //BC -1
    val1(0,0) = 0.;
    val2.Zero();
    val1(0,0) = 0;
    val1(1,1) = 1.e9;
    val2(0,0) = 0.;
    TPZMaterial * BCondD1 = material1->CreateBC(material1, bc1,2, val1, val2);
    cmesh->InsertMaterialObject(BCondD1);
    //BC -2
    val1.Zero();
    val2(0,0) = 10.;
    TPZMaterial * BCondD2 = material1->CreateBC(material1, bc2,neumann, val1, val2);
    cmesh->InsertMaterialObject(BCondD2);

    //BC -3
    val1.Zero();
    val2.Zero();
    val2(1,0) = 20.;
    TPZMaterial * BCondD3 = material1->CreateBC(material1, bc3,neumann, val1, val2);
     cmesh->InsertMaterialObject(BCondD3);

    //BC -4
    val2.Zero();
    val1(0,0) = 1.e9;
    val1(1,1) = 0.;
    val2(0,0) = 0.;
    TPZMaterial * BCondD4 = material1->CreateBC(material1, bc4,2, val1, val2);
    cmesh->InsertMaterialObject(BCondD4);

    //BC -5: dirichlet nulo
    TPZMaterial * BCondD5 = material1->CreateBC(material1, bc5,dirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCondD5);

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
    InsertMaterialObjects(cmeshauto.operator->());
    
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

