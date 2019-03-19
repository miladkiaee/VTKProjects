#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>

#include <vtkSTLWriter.h>
#include <vtkSTLReader.h>

int main(int argc, char* argv[])
{
    // PolyData to process
    std::string input_name1(argv[1]);
    std::cout << "Reading stl file : " << input_name1 << std::endl;
    vtkSmartPointer<vtkSTLReader> stlReader1 =
    vtkSmartPointer<vtkSTLReader>::New();
    stlReader1->SetFileName(input_name1.c_str());
    stlReader1->Update();
    vtkSmartPointer<vtkPolyData> pd1;
    pd1 = stlReader1->GetOutput();
    
    double ni = 2000;
    double rx = 0.1;
    double ea = 110;
    double fa = 110;
    
    std::cout << "smoothing .." << std::endl;
    
    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter =
        vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputData(pd1);
    smoothFilter->SetNumberOfIterations(ni);
    smoothFilter->SetRelaxationFactor(rx);
    smoothFilter->FeatureEdgeSmoothingOn();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->SetEdgeAngle(ea);
    smoothFilter->SetFeatureAngle(fa);
    smoothFilter->Update();
    
    std::cout << "iter = " << ni << ", relax = " << rx << ", edge angle = " << ea << ", feature angle = " << fa << std::endl;

    // Update normals on newly smoothed polydata
    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalGenerator->SetInputConnection(smoothFilter->GetOutputPort());
    normalGenerator->ComputePointNormalsOn();
    normalGenerator->ComputeCellNormalsOn();
    normalGenerator->Update();

    /////
    std::string outname(argv[2]);
    std::string outSTL = outname + ".stl";
    vtkSmartPointer<vtkSTLWriter> sw2 = vtkSmartPointer<vtkSTLWriter>::New();
    sw2->SetFileName(outSTL.c_str());
    std::cout << "writing stl .. " << std::endl;
    sw2->SetInputData(normalGenerator->GetOutput());
    sw2->Write();

    return EXIT_SUCCESS;
}
