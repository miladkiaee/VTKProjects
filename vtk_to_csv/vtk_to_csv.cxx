#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>

#include <vtkDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkTableReader.h>
#include <vtkPLYReader.h>

#include <vtkCurvatures.h>

#include <vtkLookupTable.h>
#include <vtkPointData.h>

#include <vtkDataObjectToTable.h>

#include <vtkTable.h>
#include <vtkPassArrays.h>

#include <vtkDelimitedTextWriter.h>

int main(int argc, char *argv[])
{
    
    vtkSmartPointer<vtkPLYReader> reader =
    	vtkSmartPointer<vtkPLYReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();
    
    vtkSmartPointer<vtkPassArrays> passArrays =
        vtkSmartPointer<vtkPassArrays>::New();
    passArrays->SetInputData(reader->GetOutput());
    passArrays->AddArray(vtkDataObject::POINT, "Gauss_Curvature");    
    passArrays->Update();

    std::cout << "curvature for " << argv[1] << std::endl;
    
    vtkSmartPointer<vtkDataObjectToTable> vtkTable =
    	vtkSmartPointer<vtkDataObjectToTable>::New();
    vtkTable->SetRowData(passArrays->GetOutput());
    vtkTable->SetFieldType(VTK::POINT_DATA);
    vtkTable->Update();
    
    vtkSmartPointer<vtkDelimitedTextWriter> writer =
    vtkSmartPointer<vtkDelimitedTextWriter>::New();
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(vtkTable->GetPoints());
#else
    writer->SetInputData(vtkTable->GetOutput());
#endif
    writer->SetFileName(argv[2]);
    writer->SetFieldDelimiter(" ");
    writer->Write();
    
    
    return EXIT_SUCCESS;
}
