#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include <vtkWindowedSincPolyDataFilter.h>

#include <vtkCleanPolyData.h>

#include <vtkSTLWriter.h>
#include <vtkSTLReader.h>

#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

int main(int argc, char * argv[])
{
    // PolyData to process
    std::string input_name1(argv[1]);
    std::cout << "reading stl file : " << input_name1 << std::endl;
    vtkSmartPointer<vtkXMLPolyDataReader> pd1 =
    vtkSmartPointer<vtkXMLPolyDataReader>::New();
    pd1->SetFileName(argv[1]);
    pd1->Update();
    
    double ni = 15;
    double ea = 180;
    double fa = 120;
    double pb = 0.001;
    
    vtkSmartPointer<vtkCleanPolyData> cpd =
    vtkSmartPointer<vtkCleanPolyData>::New();
    cpd->SetInputData(pd1->GetOutput());
    cpd->ConvertLinesToPointsOn();
    cpd->ConvertPolysToLinesOn();
    cpd->ConvertStripsToPolysOn();
    cpd->PointMergingOn();
    cpd->Update();
  
  vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
    vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  smoother->SetInputData(cpd->GetOutput());
  smoother->SetNumberOfIterations(ni);
  smoother->BoundarySmoothingOn();
  smoother->FeatureEdgeSmoothingOn();
  smoother->SetFeatureAngle(fa);
  smoother->SetPassBand(pb);
  smoother->NonManifoldSmoothingOn();
  smoother->NormalizeCoordinatesOn();
  smoother->Update();

    /////
    std::string outname(argv[2]);
    std::string outSTL = outname;
    vtkSmartPointer<vtkSTLWriter> sw2 = vtkSmartPointer<vtkSTLWriter>::New();
    sw2->SetFileName(outSTL.c_str());
    std::cout << "writing stl .. " << std::endl;
    sw2->SetInputData(smoother->GetOutput());
    sw2->Write();

  return EXIT_SUCCESS;
}
