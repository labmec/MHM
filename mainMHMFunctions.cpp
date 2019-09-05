
#include "pzlog.h"
#include "pzgmesh.h"
#include "pzmanvector.h"
#include "TPZVTKGeoMesh.h"

#include "TPZGmshReader.h"

#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "meshgen.h"

#ifndef USING_MKL
#include "pzskylstrmatrix.h"
#endif

TPZCompMesh *BuildComputationalMesh(TPZGeoMesh *gmesh, int pOrder);

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mainskeleton"));
#endif

int const darcy = 1;
int const frac = 2;
int const bc1=-1;
int const bcfracpoint=-2;
int const bcpoint = -6;

struct BoundaryForce
{
    int sidecount = 4;
    
    void operator()(const TPZVec<REAL> &x, TPZVec<STATE> &func)
    {
        switch(sidecount)
        {
            case 0:
                if(fabs(x[1]) < 1.e-6 && x[0] < 0.5)
                {
                    func[0] = -4.*(1.-2.*x[0]);
                }
                break;
            case 1:
                if(fabs(x[1]-1.) < 1.e-6)
                {
                    func[0] = -1.;
                }
                break;
            case 2:
                if(fabs(x[1]) < 1.e-6 && x[0] < 0.5)
                {
                    func[0] = -2.;
                }
                break;
            case 3:
                if(fabs(x[1]) < 1.e-6 && x[0] > 0.5)
                {
                    func[0] = -2.;
                }
                break;
            case 4:
                func[0] = 0.;
                break;
            case 5:
                if(fabs(x[0]-1.) < 1.e-6)
                {
                    func[0] = -1.;
                }
                break;

            default:
                DebugStop();
        }
    }
};

BoundaryForce bound;

void StupidFunc(const TPZVec<REAL> &x, TPZVec<STATE> &force)
{
    bound(x,force);
}

TPZAnalyticSolution *example;
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    TExceptionManager except;
    
    TPZGeoMesh *gmesh = 0;
    
    {
        TPZGmshReader gmsh;
        
        // Assigns IDs of 1D and 2D elements defining boundary conditions.
        gmsh.GetDimNamePhysical()[1]["BC"] = bc1;
        
        gmsh.GetDimNamePhysical()[0]["BC_POINT"] = bcpoint;
        gmsh.GetDimNamePhysical()[0]["BC_FRACPOINT"] = bcfracpoint;

        // Assigns IDs of 2D and 3D elements defining the problem domain.
        gmsh.GetDimNamePhysical()[2]["DARCY"] = 1;
        gmsh.GetDimNamePhysical()[1]["FRAC"] = 2;
        
    #ifdef MACOSX
        gmesh = gmsh.GeometricGmshMesh("../MHMDFN_Functions.msh");
    #else
        gmesh = gmsh.GeometricGmshMesh("MHMDFN_Functions.msh");
    #endif
    }
#ifdef PZDEBUG
    if(1)
    {
        std::ofstream file("GMeshDFN.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }
#endif

    TPZCompMesh *cmesh = BuildComputationalMesh(gmesh, 2);
    
    
    //calculo solution
    bool shouldrenumber = true;
    TPZAnalysis an(cmesh,shouldrenumber);
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh);
    strmat.SetNumThreads(0);
#else
    TPZSkylineStructMatrix strmat(cmesh);
    strmat.SetNumThreads(0);
#endif
    
    
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    
    for(int count = 0; count < 6; count++)
    {
        bound.sidecount = count;
        
        if(bound.sidecount == 4)
        {
            TPZMaterial *mat = cmesh->FindMaterial(bcfracpoint);
            TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
            if(!bnd) DebugStop();
            bnd->Val2()(0,0) = -1.;
        }
        an.Assemble();

        an.Solve();
        
        an.SetStep(count);
        std::string plotfile2("MHMFunc2D.vtk"),plotfile1("MHMFunc1D.vtk");
        TPZStack<std::string> scalnames,vecnames;
        scalnames.Push("Pressure");
        vecnames.Push("Flux");

        an.DefineGraphMesh(cmesh->Dimension()-1, scalnames, vecnames, plotfile1);
        an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotfile2);
        int resolution = 0;
    //    an.PostProcess(resolution,cmesh->Dimension()-1);
        an.PostProcess(resolution,cmesh->Dimension());
    }
    return 0;
    
}

/// Insert material objects for the MHM Mesh solution
TPZCompMesh *BuildComputationalMesh(TPZGeoMesh *gmesh, int pOrder)
{
    TPZCompMesh *cmeshptr = new TPZCompMesh(gmesh);
    TPZCompMesh &cmesh = *cmeshptr;
    int dim = gmesh->Dimension();
    int dirichlet = 0;
    int neumann = 1;
    
    // Creates Poisson material
    TPZMatPoisson3d *material = new TPZMatPoisson3d(darcy, dim);
    material->SetInternalFlux(-1.);
    cmesh.SetDimModel(dim);
    cmesh.InsertMaterialObject(material);
    
    {
        TPZMatPoisson3d *material = new TPZMatPoisson3d(frac, dim-1);
        TPZManVector<REAL,3> convdir(3,0.);
        material->SetParameters(100., 0., convdir);
        cmesh.InsertMaterialObject(material);
    }
    //    TPZMaterial * mat(material);
    //    cmesh->InsertMaterialObject(mat);
    
    // Inserts boundary conditions
    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
    {
        TPZMaterial * BCond = material->CreateBC(material, bc1, neumann, val1, val2);
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(StupidFunc,2);
        TPZAutoPointer<TPZFunction<STATE> > autofunc(dummy);
        BCond->SetForcingFunction(autofunc);
        cmesh.InsertMaterialObject(BCond);
    }
    {
        TPZMaterial * BCond = material->CreateBC(material, bcfracpoint, neumann, val1, val2);
        cmesh.InsertMaterialObject(BCond);
    }
    {
        TPZMaterial * BCond = material->CreateBC(material, bcpoint, dirichlet, val1, val2);
        cmesh.InsertMaterialObject(BCond);
    }
    
    cmesh.SetDefaultOrder(pOrder);
    cmesh.SetAllCreateFunctionsContinuous();
    
    // Adjusts computational data structure
    cmesh.AutoBuild();

    return cmeshptr;
}

