#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>

#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPLYReader.h>

#include <vtkCurvatures.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkScalarBarActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkColorSeries.h>
#include <vtkPointData.h>

int main(int argc, char *argv[])
{
    
    vtkSmartPointer<vtkPLYReader> reader =
    vtkSmartPointer<vtkPLYReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();
    
    vtkSmartPointer<vtkCurvatures> curvaturesFilter =
        vtkSmartPointer<vtkCurvatures>::New();
    curvaturesFilter->SetInputData(reader->GetOutput());
    curvaturesFilter->SetCurvatureType(0); // gauss 0, mean 1, max 2, min 3
    curvaturesFilter->Update();

    std::cout << "curvature for " << argv[1] << std::endl;
    
    vtkSmartPointer<vtkPolyDataWriter> writer =
    vtkSmartPointer<vtkPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(curvaturesFilter->GetOutput());
#else
    writer->SetInputData(curvaturesFilter->GetOutput());
#endif
    writer->SetFileName(argv[2]);
    writer->Write();
    
    
    return EXIT_SUCCESS;
}
