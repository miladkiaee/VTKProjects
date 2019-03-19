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

#include <vtkConnectedPointsFilter.h>

int main(int argc, char *argv[])
{
    
    vtkSmartPointer<vtkPLYReader> reader =
    vtkSmartPointer<vtkPLYReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();

   
    vtkSmartPointer<vtkConnectedPointsFilter> pointsCon =
    vtkSmartPointer<vtkConnectedPointsFilter>::New();


    pointsCon->SetInputData(reader->GetOutput());
    double r = atof(argv[3]);
    pointsCon->SetRadius(r);
    // pointsCon->SetExtractionModeToLargestRegion();
    pointsCon->Update();
    std::cout << "point radius: " << pointsCon->GetRadius() << std::endl;
    std::cout << "mode : " << pointsCon->GetExtractionModeAsString() << std::endl;
    std::cout << "num regions : " << pointsCon->GetNumberOfExtractedRegions() << std::endl;
    

    vtkSmartPointer<vtkPolyDataWriter> writer =
    vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(pointsCon->GetOutput());
    writer->SetFileName(argv[2]);
    writer->Write();
    
    return EXIT_SUCCESS;
}
