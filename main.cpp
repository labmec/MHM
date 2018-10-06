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
#include "TPZHybridizeHDiv.h"
#include "meshgen.h"

#include <iostream>
#include <string>
#include <opencv2/opencv.hpp>
#include <math.h>
#include <set>



using namespace std;
using namespace cv;

// Creating the computational H1 mesh
TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int p_order);
TPZCompMesh *CMeshFlux(TPZGeoMesh * gmesh,int pOrder);
TPZCompMesh *CMeshPressure(TPZGeoMesh * gmesh, int pOrder);
TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);
TPZGeoMesh *GeoMeshFromPng(string name);
TPZGeoMesh *GeoMeshCoarse(TPZGeoMesh *GeoMesh, double l, double h, int nx, int ny);
TPZCompMesh *CmeshMHM(TPZGeoMesh *GeoMesh,TPZGeoMesh *GeoMesh2, int refin);

int H1Test();
int MixedTest();
int MHMTest();

int main(){
    InitializePZLOG();
    TPZGeoMesh *gmsh;
    TPZGeoMesh *gmsh2;
    TPZCompMesh *cmesh;

    
//    MHMTest();
//    H1Test();
    MixedTest();
    return 0;
}

int MixedTest(){
    
    TPZGeoMesh *gmesh = GeoMeshFromPng("../small.png");
    int flux_order = 1;
    int p_order = 1;
    
    {
#ifdef PZDEBUG
        std::ofstream file("maze.txt");
        gmesh->Print(file);
        
        std::ofstream out("maze.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
#endif
    }
    
    TPZCompMesh *cmesh_flux =CMeshFlux(gmesh,flux_order);
    TPZCompMesh *cmesh_presure =CMeshPressure(gmesh,p_order);
    
    
    TPZVec<TPZCompMesh *> fmeshvec(2);
    fmeshvec[0]=cmesh_flux;
    fmeshvec[1]=cmesh_presure;
    gmesh->ResetReference();
    
    TPZCompMesh *MixedMesh = CMeshMultphysics(gmesh,fmeshvec);
    
    std::ofstream file("MixedCMesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(MixedMesh, file);
    
    std::ofstream out("MixedCMesh.txt");
    //    MixedMesh->Print(out);
    
    std::cout << "number of equations = " << MixedMesh->NEquations() << std::endl;
    
    
    
    //Solving the system:
    bool optimizeBandwidth = true;
    MixedMesh->InitializeBlock();
    
    TPZCompMesh * cmesh_m_Hybrid;
    TPZManVector<TPZCompMesh*, 3> meshvector_Hybrid(3);
    TPZHybridizeHDiv hybridizer;
    tie(cmesh_m_Hybrid, meshvector_Hybrid) = hybridizer.Hybridize(MixedMesh, fmeshvec, true, -1.);
    cmesh_m_Hybrid->InitializeBlock();
    
    bool must_opt_band_width_Q = true;
    int number_threads = 4;
    TPZAnalysis *an = new TPZAnalysis(cmesh_m_Hybrid,must_opt_band_width_Q);
    
    //
    TPZSymetricSpStructMatrix sparse_matrix(cmesh_m_Hybrid);
    TPZStepSolver<STATE> step;
    sparse_matrix.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    an->SetStructuralMatrix(sparse_matrix);
    an->SetSolver(step);
    
    // TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, MixedMesh);
    
    
    
    // Solving the LS
    an->Assemble();
    
    //    an->Solver().Matrix()->Print("j = ", std::cout,EMathematicaInput);
    //    an->Rhs().Print("r = ", std::cout,EMathematicaInput);
    
    an->Solve();
    
    
    //POS
    TPZManVector<std::string,10> scalnames(2), vecnames(1);
    vecnames[0]  = "Flux";
    
    scalnames[0] = "Pressure";
    scalnames[1] = "Permeability";
    
    
    
    const int dim = an->Mesh()->Dimension();
    int div = 0;
    
    an->DefineGraphMesh(dim,scalnames,vecnames,"hdiv.vtk");
    an->PostProcess(div,dim);
    std::cout << "Standard post-processing finished." << std::endl;
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector_Hybrid, cmesh_m_Hybrid);

    {
        std::ofstream out("cmeshHybrid.txt");
        cmesh_m_Hybrid->Print(out);
    }
//    meshvector_Hybrid[0]->Solution().Print("flux",std::cout);
    {
        std::ofstream FileStreamLines("StreamLines.nb");
        TPZStack<TPZFMatrix<REAL>> streamlines;
        TPZManVector<REAL,3> x(3,0.);
        gmesh->ResetReference();
        meshvector_Hybrid[0]->LoadReferences();
        TPZGeoEl *gelentry = FindEntry(gmesh);
        TPZManVector<REAL,3> xi(1,0.), xco(3,0.);
        for(int i=0; i<5; i++)
        {
            for(int j=0; j<2; j++)
            {
                if(j==1 && i==0) continue;
                xi[0] = pow(i/5.,0.25);
                if(j==1) xi[0] *= -1.;
                TPZStack<TStreamLine> streamstack;
                gelentry->X(xi, xco);
                TStreamLineData result;
                result = ComputeStreamLine(meshvector_Hybrid[0], xco);
                TPZFMatrix<REAL> res = result;
                streamlines.Push(res);
            }
        }
        FileStreamLines << "AllStream = Table[0,{" << streamlines.size() << "}];\n";
        for (int i=0; i<streamlines.size(); i++) {
            FileStreamLines << "AllStream[[" << i+1 << "]] =";
            streamlines[i].Print("",FileStreamLines,EMathematicaInput);
        }
    }

    
    return 0;
}

int H1Test()
{
    
    TPZGeoMesh *gmesh = GeoMeshFromPng("maze8x8.png");
    {
#ifdef PZDEBUG
        std::ofstream file("mazeh1.txt");
        gmesh->Print(file);
        
        std::ofstream out("mazeh1.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
        
#endif
    }
    
    //Creando a malla computacional
    int p_order = 1 ;
    int number_threads = 4;
    bool must_opt_band_width_Q = true;
    TPZCompMesh *cmesh = CMeshH1(gmesh,p_order);
    TPZAnalysis *an = new TPZAnalysis(cmesh,must_opt_band_width_Q);
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix struct_mat(cmesh);
    struct_mat.SetNumThreads(number_threads);
    an->SetStructuralMatrix(struct_mat);
#else
    TPZSkylineStructMatrix struct_mat(cmesh);
    struct_mat.SetNumThreads(number_threads);
    an->SetStructuralMatrix(struct_mat);
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an->SetSolver(step);
    
    
    // Solving the LS
    an->Assemble();
    an->Solve();
    
    // post-processing step
    {
        const int dim = an->Mesh()->Dimension();
        int div = 0;
        std::string plotfile = "h1_approximation.vtk";
        TPZStack<std::string> scalar_names, vec_names;
        
        scalar_names.push_back("Solution");
        vec_names.push_back("MinusKGradU");
        an->DefineGraphMesh(dim,scalar_names,vec_names,plotfile);
        an->PostProcess(div,dim);
        std::cout << "Standard post-processing finished." << std::endl;
    }
    
    return 0;
}
TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int p_order){
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1.0e-3;
    REAL perm_1 = 1000.0;
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(impervious_mat,dim);
    TPZMatPoisson3d *mat_1 = new TPZMatPoisson3d(permeable_mat,dim);
    mat_0->SetParameters(perm_0, conv, convdir);
    mat_1->SetParameters(perm_1, conv, convdir);
    
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
    cmesh->InsertMaterialObject(mat_1);
    
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    int right_bc_id = -2;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_N, val1, val2);
    cmesh->InsertMaterialObject(right_bc);
    
    int left_bc_id = -4;
    val2(0,0) = 0.0;
    TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_N, val1, val2);
    cmesh->InsertMaterialObject(left_bc);
    
    int bottom_bc_1id = -1;
    val2(0,0) = 0.0;
    TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc_1);
    
    int top_bc_1id = -3;
    val2(0,0) = 0.0;
    TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(top_bc_1);
    
    
    //Dirichlet Conditions (p=1 in, p=0 out)
    int bottom_bc_id = -5;
    val2(0,0) = 100.0;
    TPZMaterial * bottom_bc = mat_0->CreateBC(mat_0, bottom_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc);
    
    int top_bc_id = -6;
    val2(0,0) = 10.0;
    TPZMaterial * top_bc = mat_0->CreateBC(mat_0, top_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(top_bc);
    
    cmesh->SetName("LaberintoTest");
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(p_order);
    cmesh->AutoBuild();
    
    
    
#ifdef PZDEBUG
    std::ofstream file("cmesh_h.txt");
    cmesh->Print(file);
#endif
    
    return cmesh;
}

TPZCompMesh *CMeshFlux(TPZGeoMesh * gmesh,int pOrder){
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1.0e-3;
    REAL perm_1 = 10;
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    
    //Definition of the approximation space:
    
    TPZVecL2 *mat_0 = new TPZVecL2(impervious_mat);
    TPZVecL2 *mat_1 = new TPZVecL2(permeable_mat);
    
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
    cmesh->InsertMaterialObject(mat_1);
    
    cmesh->SetAllCreateFunctionsHDiv(); //Creating H(div) functions
    
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    int right_bc_id = -2;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_N, val1, val2);
    cmesh->InsertMaterialObject(right_bc);
    
    int left_bc_id = -4;
    TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_N, val1, val2);
    cmesh->InsertMaterialObject(left_bc);
    
    int bottom_bc_1id = -1;
    TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc_1);
    
    int top_bc_1id = -3;
    TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_N, val1, val2);
    cmesh->InsertMaterialObject(top_bc_1);
    
    
    //Dirichlet Conditions (p=1 in, p=0 out)
    int bottom_bc_id = -5;
    TPZMaterial * bottom_bc = mat_0->CreateBC(mat_0, bottom_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc);
    
    int top_bc_id = -6;
    TPZMaterial * top_bc = mat_0->CreateBC(mat_0, top_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(top_bc);
    
    cmesh->SetName("LaberintoTest");
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream file("cmesh_flux.txt");
    cmesh->Print(file);
#endif
    
    return cmesh;
    
}
TPZCompMesh *CMeshPressure(TPZGeoMesh * gmesh, int pOrder){
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1.0e-3;
    REAL perm_1 = 10;
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    
    TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(impervious_mat,dim);
    TPZMatPoisson3d *mat_1 = new TPZMatPoisson3d(permeable_mat,dim);
    mat_0->SetParameters(perm_0, conv, convdir);
    mat_1->SetParameters(perm_1, conv, convdir);
    
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
    cmesh->InsertMaterialObject(mat_1);
    
    
    //    int type_D = 0;
    //    int type_N = 1;
    //    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    //
    //    // Insert boundary conditions
    //    //Neumann boundary conditions (flux = 0)
    //    int right_bc_id = -2;
    //    val2(0,0) = 0.0;
    //    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_N, val1, val2);
    //    cmesh->InsertMaterialObject(right_bc);
    //
    //    int left_bc_id = -4;
    //    val2(0,0) = 0.0;
    //    TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_N, val1, val2);
    //    cmesh->InsertMaterialObject(left_bc);
    //
    //    int bottom_bc_1id = -1;
    //    val2(0,0) = 0.0;
    //    TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
    //    cmesh->InsertMaterialObject(bottom_bc_1);
    //
    //    int top_bc_1id = -3;
    //    val2(0,0) = 0.0;
    //    TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_N, val1, val2);
    //    cmesh->InsertMaterialObject(top_bc_1);
    //
    //
    //    //Dirichlet Conditions (p=1 in, p=0 out)
    //    int bottom_bc_id = -5;
    //    val2(0,0) = 100.0;
    //    TPZMaterial * bottom_bc = mat_0->CreateBC(mat_0, bottom_bc_id, type_D, val1, val2);
    //    cmesh->InsertMaterialObject(bottom_bc);
    //
    //    int top_bc_id = -6;
    //    val2(0,0) = 10.0;
    //    TPZMaterial * top_bc = mat_0->CreateBC(mat_0, top_bc_id, type_D, val1, val2);
    //    cmesh->InsertMaterialObject(top_bc);
    
    cmesh->SetName("Pressure");
    cmesh->AutoBuild();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
#ifdef PZDEBUG
    std::ofstream out("cmeshPress.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
    
}
TPZGeoMesh *GeoMeshFromPng(string name){
    const int bcDL = -1;
    const int bcB = -2;
    const int bcDR = -3;
    const int bcDT = -4;
    
    
    //  Mat image = imread("normal.png",IMREAD_GRAYSCALE);
    Mat image = imread(name,IMREAD_GRAYSCALE);
    //    Mat image = imread("single_quad.png",IMREAD_GRAYSCALE);
    int k=0;
    int px=image.size[0];
    int py=image.size[1];
    int p =px*py;
    if(p==0){
        DebugStop();
    }
    vector<int> vec(p,0);
    
    for (int i = 0; i<px; ++i) {
        for (int j = py  ; j>0; --j) {
            int val =(int)image.at<uchar>(Point(j, i));
            if (val>200){
                val=255;
            }
            int pix =val/255;
            vec[p-k]=pix;
            
            
            
            k++;
        }
    }
    
    
    
    // Creating the Geo mesh
    TPZManVector<REAL,3> x0(3,0.),x1(3,px);
    x1[2] = 0.;
    TPZManVector<int,2> nelx(2,py);
    nelx[0] = px;
    TPZGenGrid gengrid(nelx,x0,x1);
    gengrid.SetElementType(EQuadrilateral);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gengrid.Read(gmesh);
    //gengrid.Read(gmesh,2);
    
    //MatsID
    int nels = gmesh->NElements();
    TPZGeoEl *gel_in;
    TPZGeoEl *gel_out;
    TPZGeoEl *gel_in1D;
    TPZGeoEl *gel_out1D;
    
    
    for (int i=0; i<nels; i++) {
        TPZGeoEl *gel =gmesh->Element(i);
        gel->SetMaterialId(vec[i]+ 1 );
        
        if (i<= px) {
            if((vec[i]+1)==2){
                gel_in =gel;
            }
        }
        
        if (i >= (px)*(py-1)) {
            if((vec[i]+1)==2){
                gel_out=gel;
            }
        }
        
    }
    
    //gengrid.SetBC(TPZGeoMesh *gr, int side, int bc)
    gengrid.SetBC(gmesh, 4, bcDL);
    gengrid.SetBC(gmesh, 5, bcB);
    gengrid.SetBC(gmesh, 6, bcDR);
    gengrid.SetBC(gmesh, 7, bcDT);
    
    
    int gel_in_index = gel_in->Index();
    gel_in1D = gmesh->Element(gel_in_index)->Neighbour(4).Element();
    gel_in1D->SetMaterialId(-5);
    
    int gel_out_index = gel_out->Index();
    gel_out1D = gmesh->Element(gel_out_index)->Neighbour(6).Element();
    gel_out1D->SetMaterialId(-6);
    
    
    gmesh->BuildConnectivity();
    return gmesh;
}
TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec){
    
    //Creating computational mesh for multiphysic elements
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    REAL perm_0 = 1.0e-3;
    REAL perm_1 = 1000.0;
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
    std::cout<<mphysics->NMaterials();
    
    TPZMixedPoisson *mat_0 = new TPZMixedPoisson(impervious_mat,dim);
    mat_0->SetPermeability(perm_0);
    
    
    TPZMixedPoisson *mat_1 = new TPZMixedPoisson(permeable_mat,dim);
    mat_1->SetPermeability(perm_1);
    
    mat_0->SetParameters(perm_0, conv, convdir);
    mat_1->SetParameters(perm_1, conv, convdir);
    
    mphysics->InsertMaterialObject(mat_0);
    mphysics->InsertMaterialObject(mat_1);
    
    //Inserir condicoes de contorno
    
    
    int type_D = 0;
    int type_N = 1;
    TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
    
    // Insert boundary conditions
    //Neumann boundary conditions (flux = 0)
    int right_bc_id = -2;
    val2(0,0) = 0.0;
    TPZMaterial * right_bc = mat_0->CreateBC(mat_0, right_bc_id, type_N, val1, val2);
    mphysics->InsertMaterialObject(right_bc);
    
    
    int left_bc_id = -4;
    val2(0,0) = 0.0;
    TPZMaterial * left_bc = mat_0->CreateBC(mat_0, left_bc_id, type_N, val1, val2);
    mphysics->InsertMaterialObject(left_bc);
    
    int bottom_bc_1id = -1;
    val2(0,0) = 0;
    TPZMaterial * bottom_bc_1 = mat_0->CreateBC(mat_0, bottom_bc_1id, type_N, val1, val2);
    mphysics->InsertMaterialObject(bottom_bc_1);
    
    int top_bc_1id = -3;
    val2(0,0) = 0.0;
    TPZMaterial * top_bc_1 = mat_0->CreateBC(mat_0, top_bc_1id, type_N, val1, val2);
    mphysics->InsertMaterialObject(top_bc_1);
    
    //Dirichlet Conditions (p=1 in, p=0 out)
    int bottom_bc_id = -5;
    val2(0,0) = 100.0;
    TPZMaterial * bottom_bc = mat_0->CreateBC(mat_0, bottom_bc_id, type_D, val1, val2);
    mphysics->InsertMaterialObject(bottom_bc);
    
    int top_bc_id = -6;
    val2(0,0) = 10.0;
    TPZMaterial * top_bc = mat_0->CreateBC(mat_0, top_bc_id, type_D, val1, val2);
    mphysics->InsertMaterialObject(top_bc);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->SetDimModel(gmesh->Dimension());
    mphysics->AutoBuild();
    
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
#ifdef PZDEBUG
    std::ofstream file("cmesh_mphysics.txt");
    mphysics->Print(file);
#endif
    
    return mphysics;
}
TPZGeoMesh *GeoMeshCoarse(TPZGeoMesh *GeoMesh, double l, double h, int nx, int ny){
    const int bcDL = -1;
    const int bcB = -2;
    const int bcDR = -3;
    const int bcDT = -4;
    
    // Creating the Geo mesh
    TPZManVector<REAL,3> x0(3,0.),x1(3,l);
    x1[1] = h;
    x1[2] = 0.;
    TPZManVector<int,2> nelx(2,ny);
    nelx[0] = nx;
    TPZGenGrid gengrid(nelx,x0,x1);
    gengrid.SetElementType(EQuadrilateral);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gengrid.Read(gmesh);
    //gengrid.Read(gmesh,2);

//    TPZVec<TPZGeoEl *> new_elems;
//    gmesh->Element(0)->Divide(new_elems);
//    for (auto &elem : new_elems) {
//        TPZVec<TPZGeoEl *> new_elems2;
//        elem->Divide(new_elems2);
//    }
    gmesh->BuildConnectivity();
    return gmesh;
}
TPZCompMesh *CmeshMHM(TPZGeoMesh *GeoMeshCoarse, TPZGeoMesh *GeoMeshFine, int nref){
    int dim = GeoMeshCoarse->Dimension();
    
    TPZCompMesh *cmesh = new TPZCompMesh(GeoMeshCoarse);
    cmesh->SetDimModel(dim);
    
    int impervious_mat = 1;
    int permeable_mat = 2;

    
    REAL perm_0 = 1.0e-3;
    REAL perm_1 = 1000.0;
    
    REAL conv=0;
    TPZVec<REAL> convdir(dim,0.0);
    
  
    cmesh->SetDimModel(dim);
    
    TPZMatPoisson3d *mat_0 = new TPZMatPoisson3d(impervious_mat,dim);
    TPZMatPoisson3d *mat_1 = new TPZMatPoisson3d(permeable_mat,dim);
    mat_0->SetParameters(perm_0, conv, convdir);
    mat_1->SetParameters(perm_1, conv, convdir);
    
    //  inserting volumetric materials objects
    cmesh->InsertMaterialObject(mat_0);
    cmesh->InsertMaterialObject(mat_1);
    cmesh->AutoBuild();
    TPZVec<REAL> qsi(3,0);
    TPZVec<REAL> result(3,0);
    
   

    TPZStack<TPZVec<int64_t>> vecs;
    TPZVec<int64_t> indexf;
    for(int i=0; i<nref; i++){
         int nel = cmesh->NElements();
        for(int i=0; i<nel; i++){
            if (cmesh->Element(i)){
            cmesh->Element(i)->Divide(i,indexf );
            vecs.push_back(indexf);
            }
        }
    }
    

    
   int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        if (cmesh->Element(i)){
            TPZGeoEl *gel = cmesh->Element(i)->Reference();
            TPZFMatrix<REAL> cooridnates1(3,4);
            TPZVec<REAL> qsi(3,0);
            TPZVec<REAL> result(3,0);
            gel->X(qsi,result);
            int flor =floor(result[0]);
            int y =floor(result[1])*8;
            int pos = flor + y;
            for(int j=0; j<nel; j++){
                TPZGeoEl *gel2 = GeoMeshFine->Element(pos);
                int matid= gel2->MaterialId();
                gel->SetMaterialId(matid);
            }
        }
    }
    
    TPZVec<int64_t> subindex;
//    cmesh->Element(0)->Divide(0,subindex );
    int nel2 = cmesh->NElements();
  //  cmesh->AutoBuild();
    std::ofstream out("mazehmhmCOMP.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out);
    

//
//
//    std::cout << coord(0,0) << " " << coord(1,0) << std::endl;
//    std::cout << coord(0,1) << " " << coord(1,1) << std::endl;
//    std::cout << coord(0,2) << " " << coord(1,2) << std::endl;
//    std::cout << coord(0,3) << " " << coord(1,3) << std::endl;
    //cmesh->Print();
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out);
    
    
    return 0;
}
int MHMTest(){
    
    TPZGeoMesh *gmesh = GeoMeshFromPng("maze8x8.png");
    TPZGeoMesh *gmeshMHM = GeoMeshCoarse(gmeshMHM, 8, 8, 2, 2);
    TPZCompMesh *cmesh = CmeshMHM(gmeshMHM,gmesh,2);
    std::ofstream out("mazehmhm.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmeshMHM, out, true);
    std::ofstream out2("mazeh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2, true);
    return 0;
}
