#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkPointData.h>
#include <vtkReverseSense.h>
#include <vtkFloatArray.h>

#include <vtkSTLWriter.h>
#include <vtkSTLReader.h>

int main(int argc, char * argv[])
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


  vtkSmartPointer<vtkReverseSense> reverseSense =
    vtkSmartPointer<vtkReverseSense>::New();
  reverseSense->SetInputData(pd1);
  reverseSense->ReverseNormalsOn();
  reverseSense->Update();


    /////
    std::string outname(argv[2]);
    std::string outSTL = outname;
    vtkSmartPointer<vtkSTLWriter> sw2 = vtkSmartPointer<vtkSTLWriter>::New();
    sw2->SetFileName(outSTL.c_str());
    std::cout << "writing stl .. " << std::endl;
    sw2->SetInputData(reverseSense->GetOutput());
    sw2->Write();

  return EXIT_SUCCESS;
}
