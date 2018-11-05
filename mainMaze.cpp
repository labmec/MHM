#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <fstream>
#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"
#include "pzreal.h"
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
#include "TPZMHMixedMeshChannelControl.h"

#include "TPZMHMixedHybridMeshControl.h"
#include "TPZHybridizeHDiv.h"
#include "meshgen.h"
#include "ConfigCasesMaze.h"
#include <iostream>
#include <string>
#include <opencv2/opencv.hpp>
#include <math.h>
#include <set>

#define new_identifier


using namespace std;
using namespace cv;

// Creating the computational H1 mesh
TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int p_order, ConfigCasesMaze Conf);

// Creating the computational flux mesh
TPZCompMesh *CMeshFlux(TPZGeoMesh * gmesh,int pOrder);

// Creating the computational pressure mesh
TPZCompMesh *CMeshPressure(TPZGeoMesh * gmesh, int pOrder,ConfigCasesMaze Conf);

// Creating the computational multphysics mesh
TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec,ConfigCasesMaze Conf);

// Read a mesh from a png file. The size of the domain will be npix_x by npix_y (read from the image) . Return l=nx; h=ny.
TPZGeoMesh *GeoMeshFromPng(string name, double &l, double &h);

// Create a geometric mesh with the given parameters, nx and ny are the coarse elements number. The total number of elements are defined by the image read.
TPZGeoMesh *GenerateGeoMesh(string name, int nx, int ny);

// Create a geoElSide map to be open because it has flux.
// cmesh is the flux mesh obtained from a Hdiv previus simulation.
std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> IdentifyChanel (TPZCompMesh *cmesh);


// Compute the geometric mesh coarse indices
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);

// Insert the necessary objects material in the computational mesh
void InsertMaterialObjects(TPZMHMixedMeshControl &control);

// Solve the H1 problem with "Conf" configuration
// Conf contains the maze information and the problem boundary conditions
int H1Test(ConfigCasesMaze Conf);

// Solve the mixed problem with "Conf" configuration
// Conf contains the maze information and the problem boundary conditions
TPZCompMesh* MixedTest(ConfigCasesMaze Conf);

// Solve the maze using MHM. By default (2x2 coarse elements)
// Conf contains the maze information and the problem boundary conditions
int MHMTest(ConfigCasesMaze Conf);

int main(){
    InitializePZLOG();
    
    ConfigCasesMaze ConfCasesMeze;
    ConfCasesMeze.SetImageName("MAZE4X4C.png");
    ConfCasesMeze.SetImperviousMatPermeability(1);
    ConfCasesMeze.SetPermeableMatPermeability(1000000);
    ConfCasesMeze.SetFluxOrder(1);
    ConfCasesMeze.SetPressureOrder(1);
    ConfCasesMeze.SetCCPressureIn(100);
    ConfCasesMeze.SetCCPressureOut(1);
    ConfCasesMeze.SetMHMOpenChannel(false);
    ConfCasesMeze.SetVTKName("MHMOmar.vtk");
 
//    H1Test(ConfCasesMeze);
//    MixedTest(ConfCasesMeze);
      MHMTest(ConfCasesMeze);

      return 0;
    
    
}


TPZCompMesh* MixedTest(ConfigCasesMaze Conf){
  
    TPZGeoMesh *gmesh = GenerateGeoMesh(Conf.GetImageName(),2,2);
    int flux_order = Conf.GetFluxOrder();
    int p_order = Conf.GetPressureOrder();
    
    {
#ifdef PZDEBUG
        std::ofstream file("maze.txt");
        gmesh->Print(file);
        
        std::ofstream out("maze.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
#endif
    }
    
    TPZCompMesh *cmesh_flux =CMeshFlux(gmesh,flux_order);
    TPZCompMesh *cmesh_presure =CMeshPressure(gmesh,p_order,Conf);
    
    
    TPZVec<TPZCompMesh *> fmeshvec(2);
    fmeshvec[0]=cmesh_flux;
    fmeshvec[1]=cmesh_presure;
    gmesh->ResetReference();
    
    TPZCompMesh *MixedMesh = CMeshMultphysics(gmesh,fmeshvec,Conf);
    
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
    an->Assemble();
    an->Solve();
    
   
    
    //POS
    TPZManVector<std::string,10> scalnames(2), vecnames(1);
    vecnames[0]  = "Flux";
    
    scalnames[0] = "Pressure";
    scalnames[1] = "Permeability";
    
    
    const int dim = an->Mesh()->Dimension();
    int div = 0;
    
    an->DefineGraphMesh(dim,scalnames,vecnames,Conf.GetVTKName());
    an->PostProcess(div,dim);
    std::cout << "Standard post-processing finished." << std::endl;
    
    return cmesh_flux;
}


int H1Test(ConfigCasesMaze Conf)
{
    double l;
    double h;
    TPZGeoMesh *gmesh = GeoMeshFromPng(Conf.GetImageName(),l,h);
    {
#ifdef PZDEBUG
        std::ofstream file("mazeh1.txt");
        gmesh->Print(file);
        
        std::ofstream out("mazeh1.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
        
#endif
    }
    
    //Creando a malla computacional
    int p_order = Conf.GetPressureOrder() ;
    int number_threads = 4;
    bool must_opt_band_width_Q = true;
    TPZCompMesh *cmesh = CMeshH1(gmesh,p_order,Conf);
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
        std::string plotfile = Conf.GetVTKName();
        TPZStack<std::string> scalar_names, vec_names;
        
        scalar_names.push_back("Solution");
        vec_names.push_back("MinusKGradU");
        an->DefineGraphMesh(dim,scalar_names,vec_names,plotfile);
        an->PostProcess(div,dim);
        std::cout << "Standard post-processing finished." << std::endl;
    }
    
    return 0;
}
TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int p_order, ConfigCasesMaze Conf){
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    REAL perm_0 = Conf.GetImperviousMatPermeability();
    REAL perm_1 = Conf.GetPermeableMatPermeability();
    
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
    val2(0,0) = Conf.GetCCPressureIn();
    TPZMaterial * bottom_bc = mat_0->CreateBC(mat_0, bottom_bc_id, type_D, val1, val2);
    cmesh->InsertMaterialObject(bottom_bc);
    
    int top_bc_id = -6;
    val2(0,0) = Conf.GetCCPressureOut();
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
TPZCompMesh *CMeshPressure(TPZGeoMesh * gmesh, int pOrder,ConfigCasesMaze Conf){
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    REAL perm_0 = Conf.GetImperviousMatPermeability();
    REAL perm_1 = Conf.GetPermeableMatPermeability();
    
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
TPZGeoMesh *GeoMeshFromPng(string name, double &l, double &h){
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
    l=px;
    h=py;
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
TPZCompMesh *CMeshMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, ConfigCasesMaze Conf ){
    
    //Creating computational mesh for multiphysic elements
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int impervious_mat = 1;
    int permeable_mat = 2;
    int dim = gmesh->Dimension();
    
    REAL perm_0 = Conf.GetImperviousMatPermeability();
    REAL perm_1 = Conf.GetPermeableMatPermeability();
    
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
    val2(0,0) = Conf.GetCCPressureIn();
    TPZMaterial * bottom_bc = mat_0->CreateBC(mat_0, bottom_bc_id, type_D, val1, val2);
    mphysics->InsertMaterialObject(bottom_bc);
    
    int top_bc_id = -6;
    val2(0,0) = Conf.GetCCPressureOut();
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
int MHMTest(ConfigCasesMaze Conf){
    
    TRunConfig Configuration;
  
    TPZGeoMesh *gmeshcoarse =GenerateGeoMesh(Conf.GetImageName(), 2, 2);
//    std::ofstream file("mesh4x4.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmeshcoarse, file);
//
//
//    std::ofstream out2("mesh4x4.txt");
//    gmeshcoarse->Print(out2);
    
    int interface_mat_id = 600;
    int flux_order = 1;
    int p_order = 1;
    bool OpenChannel = Conf.GetMHMOpenChannel();
    
    TPZAutoPointer<TPZMHMixedMeshChannelControl> MHMixed;
  
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmeshcoarse);
        {
            std::ofstream out("gmeshauto.txt");
            gmeshauto->Print(out);
        }
        TPZMHMixedMeshChannelControl *mhm = new TPZMHMixedMeshChannelControl(gmeshauto);
        TPZVec<int64_t> coarseindices;
        ComputeCoarseIndices(gmeshauto.operator->(), coarseindices);
//        for(int i =0; i < coarseindices.size(); i++){
//            std::cout << "i = " << coarseindices[i] << std::endl;
//        }
        gmeshauto->AddInterfaceMaterial(1, 2, interface_mat_id);
        gmeshauto->AddInterfaceMaterial(2, 1, interface_mat_id);
        
        
        // criam-se apenas elementos geometricos
        mhm->DefinePartitionbyCoarseIndices(coarseindices);
        //        MHMMixedPref << "MHMixed";
        MHMixed = mhm;
        
        TPZMHMixedMeshChannelControl &meshcontrol = *mhm;
        {
            std::set<int> matids;
            matids.insert(1);
            matids.insert(2);
            mhm->fMaterialIds = matids;
            matids.clear();
            matids.insert(-1);
            matids.insert(-2);
            matids.insert(-3);
            matids.insert(-4);
            matids.insert(-5);
            matids.insert(-6);
            mhm->fMaterialBCIds = matids;
        }
        
        
        InsertMaterialObjects(*mhm);
        
#ifdef PZDEBUG2
        if(1)
        {
            std::ofstream out("MixedMeshControlHDiv.txt");
            meshcontrol.Print(out);
        }
#endif
        meshcontrol.SetInternalPOrder(2);
        meshcontrol.SetSkeletonPOrder(1);
        
        meshcontrol.DivideSkeletonElements(0);
        meshcontrol.DivideBoundarySkeletonElements();
        
//        std::ofstream file_geo("geometry.txt");
//        meshcontrol.GMesh()->Print(file_geo);
//
        bool substructure = true;
//        std::ofstream filee("Submesh.txt");
//        meshcontrol.CMesh()->Print(filee);
        std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> test;
        if(OpenChannel){
        TPZCompMesh *flux_temp = MixedTest(Conf);
        test = IdentifyChanel(flux_temp);
        }
    
     meshcontrol.BuildComputationalMesh(substructure,OpenChannel,test);
        
        //TEST
        
        
        
        //neoTest
        //numero de elementos comp
      
        
//
//        TPZVec<TPZSubCompMesh *> submeshes;
//        submeshes.Resize(nsubmeshes);
//
//        std::map<int64_t,int64_t>::iterator it;
//        int64_t count = 0;
//        for (it = meshcontrol.Coarse_to_Submesh().begin(); it != meshcontrol.Coarse_to_Submesh().end(); it++)
//        {
//
//            int64_t skelindex = it->first;
//            int64_t MHM_index = it->second;
//            if (meshcontrol.Coarse_to_Submesh().find(skelindex) == meshcontrol.Coarse_to_Submesh().end()) {
//                DebugStop();
//            }
//            int64_t MHMixed_index = meshcontrol.Coarse_to_Submesh()[skelindex];
//            int64_t vals= meshcontrol.Coarse_to_Submesh()[skelindex];
//            std::cout<<vals<<"\n";
//            //            if (MHM_index == -1 || MHMixed_index == -1) {
//            //                continue;
//            //            }
//
//            TPZSubCompMesh *submhm = dynamic_cast<TPZSubCompMesh *>(meshcontrol.CMesh()->Element(MHMixed_index));
//
//            submeshes[count]=submhm;
//
//            {
//                std::ofstream file("Submesh.vtk");
//                TPZVTKGeoMesh::PrintCMeshVTK(submhm, file);
//                //    std::ofstream out("Submesh0.txt");
//                submhm->Print(std::cout );
//            }
//            count++;
//
//
//        }
//
//        //    ComputationalMesh->Reference()->ResetReference();
//        //    ComputationalMesh->LoadReferences();
//        //    return 0;
//
//        std::map<int64_t, std::pair<int64_t, int64_t>> g_indexes,c_indexes;
        //IdentifySkeletonChanelIntersections(submeshes,g_indexes,c_indexes);
       
        
        //GlueChannelConnects(submeshes, g_indexes ,c_indexes);
        
        
        
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream file("GMeshControlHDiv.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file);
        }
#endif
#ifdef PZDEBUG2
        if(1)
        {
            std::ofstream out("MixedMeshControlHDiv.txt");
            meshcontrol.Print(out);
        }
#endif
        
        std::cout << "MHM Hdiv Computational meshes created\n";
#ifdef PZDEBUG2
        if(1)
        {
            std::ofstream gfile("geometryMHMHdiv.txt");
            gmeshauto->Print(gfile);
            std::ofstream out_mhm("MHM_hdiv.txt");
            meshcontrol.CMesh()->Print(out_mhm);
            
        }
#endif
        
        std::cout << "Number of equations MHMixed " << MHMixed->CMesh()->NEquations() << std::endl;
        
        
    }
    
    TPZCompMesh *MixedMesh = MHMixed->CMesh().operator->();
    
    //    MixedMesh->Print(out);
    
//    std::cout << "number of equations = " << MixedMesh->NEquations() << std::endl;
    
    SolveProblem(MHMixed->CMesh(), MHMixed->GetMeshes(), 0,  Conf.GetVTKName(), Configuration);
    
    
    
    /*
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
     
     */
    return 0;
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
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat);
    
    TPZMixedPoisson * mat_2 = new TPZMixedPoisson(2,dim);
    mat_2->SetSymmetric();
    mat_2->SetPermeability(1000000.);
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat_2);
    
   // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, -1, typeFlux, val1, val2Flux);
    //    bcN->SetForcingFunction(0, force);
    
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    bcN = mat->CreateBC(mat, -3, typeFlux, val1, val2Flux);
    //    bcN->SetForcingFunction(0, force);
    
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, -2, typeFlux, val1, val2Flux);
    
    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    bcS = mat->CreateBC(mat, -4, typeFlux, val1, val2Flux);
    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    
    TPZBndCond * bcIn = mat->CreateBC(mat, -5, typePressure, val1, val2Pressure);
    
    MixedFluxPressureCmesh->InsertMaterialObject(bcIn);

    TPZBndCond * bcOut = mat->CreateBC(mat, -6, typePressure, val1, val2Flux);
    
    MixedFluxPressureCmesh->InsertMaterialObject(bcOut);
    
}

// compute the coarse indices of the geometric mesh
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices)
{
//    {
//        std::ofstream out("gmeshref.txt");
//        gmesh->Print(out);
//    }
    coarseindices.Resize(gmesh->NElements());
    int count = 0;
    for (int64_t el=0; el<gmesh->NElements(); el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->Dimension() != gmesh->Dimension()) continue;
        if(gel->Father()) continue;
        coarseindices[count] = el;
        count++;
    }
    coarseindices.Resize(count);
}
TPZGeoMesh *GenerateGeoMesh(string name, int nx, int ny){
    
    double l;
    double h;
    TPZGeoMesh *FineMesh = GeoMeshFromPng(name,l,h);
    
    std::ofstream file_base_vtk("base.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(FineMesh, file_base_vtk);
    
    std::ofstream file_base_txt("base.txt");
    FineMesh->Print(file_base_txt);
    
    // Creating the Geo mesh
    TPZManVector<REAL,3> x0(3,0.),x1(3,l);
    x1[1] = h;
    x1[2] = 0.;
    TPZManVector<int,2> nelx(2,ny);
    nelx[0] = nx;
    TPZGenGrid gengrid(nelx,x0,x1);
    gengrid.SetElementType(EQuadrilateral);
    TPZGeoMesh *gmeshcoarse = new TPZGeoMesh;
    gmeshcoarse->SetDimension(2);
        gengrid.SetRefpatternElements(true);
    gengrid.Read(gmeshcoarse);
    //gengrid.Read(gmesh,2);

    
    //gengrid.SetBC(TPZGeoMesh *gr, int side, int bc)
    gengrid.SetBC(gmeshcoarse, 4, -1);
    gengrid.SetBC(gmeshcoarse, 5, -2);
    gengrid.SetBC(gmeshcoarse, 6, -3);
    gengrid.SetBC(gmeshcoarse, 7, -4);
    
    gmeshcoarse->BuildConnectivity();
    
    std::ofstream file_base_c_vtk("base_c.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmeshcoarse, file_base_c_vtk);
    
    std::ofstream file_base_c_txt("base_c.txt");
    gmeshcoarse->Print(file_base_c_txt);
    
    //Refine
    TPZVec<REAL> qsi(3,0);
    TPZVec<REAL> result(3,0);
    TPZStack<TPZVec<int64_t>> vecs;
    TPZVec<TPZGeoEl*> indexf;
//    int nref=6; para 128x128 com coarse 2x2
    int nref = log2(l) - log2(nx);
    for(int i=0; i<nref; i++){
        int nel = gmeshcoarse->NElements();
        for(int i=0; i<nel; i++){
            TPZGeoEl * gel = gmeshcoarse->Element(i);
            if (!gel or gel->HasSubElement()) {
                continue;
            }
            gel->Divide(indexf);
        }
    }
  
    
    //
    int nel = gmeshcoarse->NElements();
    
    for(int i=0; i<nel; i++){
        TPZGeoEl *gel = gmeshcoarse->Element(i);
        if (!gel or gel->HasSubElement()) {
            continue;
        }
        
        if(gel->Dimension()==2){
            TPZFMatrix<REAL> cooridnates1(3,4);
            TPZVec<REAL> qsi(3,0);
            TPZVec<REAL> result(3,0);
            gel->X(qsi,result);
            int flor =floor(result[0]);
            int y =floor(result[1])*l;
            int pos = flor + y;
            TPZGeoEl *gel2 = FineMesh->Element(pos);
            if(!gel2){
                DebugStop();
            }
            int matid= gel2->MaterialId();
            gel->SetMaterialId(matid);
            
            if(y==0 && matid==2){
                TPZGeoEl *el1D = gel->Neighbour(4).Element();
                el1D->SetMaterialId(-5);
            }
            int niv =y/l;
            if(niv==(l-1) && matid==2){
                TPZGeoEl *el1D = gel->Neighbour(6).Element();
                el1D->SetMaterialId(-6);
            }
            
        }
    }
    
    std::ofstream out("mazefine.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(FineMesh, out, true);
    
    std::ofstream out2("mazehcoarse.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmeshcoarse, out2, true);
    
    return gmeshcoarse;
};

std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> IdentifyChanel (TPZCompMesh *cmesh){
    
    int nelements = cmesh->NElements();
    TPZGeoElSide first_gelside;
    for(int iel=0; iel<=nelements; iel++){
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel){continue;}
        TPZGeoEl *gel = cel->Reference();
        TPZGeoElSide gelside(gel,4);
        TPZGeoEl *neig = gelside.Neighbour().Element();
        if(neig->MaterialId()==-5){
            first_gelside.SetElement(gel);
            first_gelside.SetSide(4);
            break;
        }
    }
    
    int count=0;
    int count_chain=0;
    double Flux_Max = 0.0;
    std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> chain;
    bool exit=false;
    
    while(exit==false){
            
        TPZGeoElSide exit_test = first_gelside;
        exit_test.SetSide(6);
        TPZGeoElSide candidate_exist =exit_test.Neighbour();
        if(candidate_exist.Element()->MaterialId()==-6){
            std::cout<<"Cadena encontrada con exito";
            exit = true;
            break;
        }
        
        
//    while(candidate.Element()->Dimension()!=2){
//
//            candidate=candidate.Neighbour();
//
//    }
    int side_in = first_gelside.Side();
    for(int ican=4; ican<8; ican++){
        first_gelside.SetSide(ican);
        TPZGeoElSide candidate = first_gelside.Neighbour();
        while(candidate.Element()->Dimension()!=2){
           
            candidate=candidate.Neighbour();
        }
        if(side_in == ican ){continue;}
        
        //calcula los 3 geoelement side y el fluxo
        TPZVec<REAL> qsi(2);
        qsi[0]=0;
        qsi[1]=0;
        int var=31;
        TPZVec<STATE> sol;
        TPZCompEl *cel = candidate.Element()->Reference();
        cel->Solution(qsi, var,sol);
        double Flux_can_Mag = sqrt(sol[0]*sol[0] + sol[1]*sol[1]);
      
        if(Flux_can_Mag > Flux_Max){
            Flux_Max =Flux_can_Mag;
            first_gelside.SetSide(ican);
            chain[count].first= first_gelside;
            chain[count].second =candidate;
        }
    }
            Flux_Max=0.0;
            first_gelside =chain[count].second;
            count++;
}
    
    //here
    
        cmesh->LoadReferences();
        TPZGeoMesh *gmesh =cmesh->Reference();
        for(auto it:gmesh->ElementVec()){
            int father_index = it->FatherIndex();
            int element = it->Index();
            std::cout<<"Element: "<<element<<" father: "<<father_index<<std::endl;
    
        }
    
    
        int n_el_chain = chain.size();
        std::map<int,pair<TPZGeoElSide, TPZGeoElSide>> skelchanel;
        int count_skel_chanel=0;
        int matId_skel_chanel = 10;
        count =0;
        for(auto it : chain){
            TPZGeoEl *first_element = it.second.first.Element();
            TPZGeoEl *second_element = it.second.second.Element();
    
            int first_father_index = first_element->LowestFather()->Index();
            int second_father_index = second_element->LowestFather()->Index();

            if(first_father_index!=second_father_index){
                TPZGeoElBC(it.second.first, matId_skel_chanel);
                skelchanel[count_skel_chanel] = it.second;
                count_skel_chanel++;
            }
    
        }
        std::ofstream out("TestMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
        return skelchanel;
    
}
