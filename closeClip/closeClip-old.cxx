/**
***********************************************
author: Milad Kiaee, kiaeedar@ualberta

discription: vtk code to close-clip a surface mesh (STL)
the resulting features are being connected to provide
a closed surtace.
usuage: exec ./closeClip input.stl clipper.vtk
***********************************************
**/

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPolyData.h>
#include <vtkStripper.h>
#include <vtkTriangleStrip.h>
#include <vtkFeatureEdges.h>
#include <vtkClipPolyData.h>
#include <vtkBox.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPLYWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkReverseSense.h>
#include <vtkIdFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkImplicitDataSet.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDataSetWriter.h>
#include <vtkCellArray.h>
#include <vtkRuledSurfaceFilter.h>

#include <vtkLine.h>
#include <vtkPolygon.h>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkCamera.h>

void WritePolyData(vtkPolyData* const polyData, const std::string& filename);
void WriteDataSet(vtkDataSet* const dataSet, const std::string& filename);
void writeBoundaryPoints( vtkSmartPointer<vtkPoints> p, std::string name );
void writeBoundaryPolyData( vtkSmartPointer<vtkPolyData> d , std::string name );
void writeClosedClipSTL( vtkSmartPointer<vtkPolyData> d, std::string name );

int main (int argc, char *argv[])
{
  // clip box and polygons
  vtkSmartPointer<vtkBox> box = vtkSmartPointer<vtkBox>::New();
  vtkSmartPointer<vtkPolygon> leftPolygon = vtkSmartPointer<vtkPolygon>::New();
  vtkSmartPointer<vtkPolygon> rightPolygon = vtkSmartPointer<vtkPolygon>::New();
  vtkSmartPointer<vtkPolygon> topPolygon = vtkSmartPointer<vtkPolygon>::New();
  vtkSmartPointer<vtkPolygon> bottomPolygon = vtkSmartPointer<vtkPolygon>::New();

  vtkSmartPointer<vtkPoints> leftPoints = vtkSmartPointer<vtkPoints>::New();
  leftPoints->InsertNextPoint( 0.02, 0.03, 0.02);
  leftPoints->InsertNextPoint( 0.02, 0.03, 0.04);
  leftPoints->InsertNextPoint( -0.02, 0.03, 0.04);
  leftPoints->InsertNextPoint( -0.02, 0.03, 0.02);

  leftPolygon->GetPoints()->DeepCopy(leftPoints);
  leftPolygon->GetPointIds()->SetNumberOfIds(4); // 4 corners of the square
  leftPolygon->GetPointIds()->SetId(0, 0);
  leftPolygon->GetPointIds()->SetId(1, 1);
  leftPolygon->GetPointIds()->SetId(2, 2);
  leftPolygon->GetPointIds()->SetId(3, 3);


  box->SetBounds(-0.02, 0.02, 0.03, 0.05, 0.02, 0.04);


  // PolyData to process
  std::string input_name(argv[1]);
  std::cout << "close clip : " << input_name << std::endl;
  vtkSmartPointer<vtkSTLReader> stlReader =
    vtkSmartPointer<vtkSTLReader>::New();
  stlReader->SetFileName(input_name.c_str());
  stlReader->Update();;
  vtkSmartPointer<vtkPolyData> poly =
    vtkSmartPointer<vtkPolyData>::New();;
  poly = stlReader->GetOutput(); //or try shallowcopy

  // clean polyData
  vtkSmartPointer<vtkCleanPolyData> cpd =
    vtkSmartPointer<vtkCleanPolyData>::New();
  cpd->SetInputData(poly);
  cpd->Update();

  std::cout << "-- The geometry has " << cpd->GetOutput()->GetNumberOfPoints() <<
  " points." << " and " << cpd->GetOutput()->GetNumberOfCells() << " cells. " <<
  std::endl;

  // Add ids to the points and cells of the sphere
  vtkSmartPointer<vtkIdFilter> cellIdsFilter =
    vtkSmartPointer<vtkIdFilter>::New();
  cellIdsFilter->SetInputConnection(cpd->GetOutputPort());
  cellIdsFilter->SetCellIds(true);
  cellIdsFilter->SetPointIds(false);
  cellIdsFilter->SetIdsArrayName("CellIds");
  cellIdsFilter->Update();
  WriteDataSet( cellIdsFilter->GetOutput(), "CellIdFilter.vtp" );

  vtkSmartPointer<vtkIdFilter> pointIdsFilter =
    vtkSmartPointer<vtkIdFilter>::New();
  pointIdsFilter->SetInputConnection(cellIdsFilter->GetOutputPort());
  pointIdsFilter->SetCellIds(false);
  pointIdsFilter->SetPointIds(true);
  pointIdsFilter->SetIdsArrayName("PointIds");
  pointIdsFilter->Update();
  vtkDataSet* cpdWithIds = pointIdsFilter->GetOutput();
  WriteDataSet( cpdWithIds, "BothIdsData.vtp" );

  // vtk clip polydata
  vtkSmartPointer<vtkClipPolyData> clip =
    vtkSmartPointer<vtkClipPolyData>::New();
  clip->SetInputData( cpdWithIds ); // ***** adds the data with id instead of cpd
  clip->SetClipFunction(box);
  clip->SetValue(0);
  clip->Update();

  WriteDataSet( clip->GetOutput(), "clip.vtp" );
  vtkPolyData* pClip = clip->GetOutput();
  std::cout << "-- clipped has " << pClip->GetNumberOfPoints() <<
    " points." << std::endl;
  std::cout << "-- clipped has " << pClip->GetNumberOfCells() <<
    " cells." << std::endl;

  vtkIdTypeArray* pClipCellIds =
    vtkIdTypeArray::SafeDownCast(pClip->GetCellData()->GetArray("CellIds"));

  // ********************************************* //
  // extract feature edges
  vtkSmartPointer<vtkFeatureEdges> fEdges =
    vtkSmartPointer<vtkFeatureEdges>::New();
  fEdges->SetInputData( pClip ); // input is already labled
  fEdges->BoundaryEdgesOn();
  fEdges->FeatureEdgesOff();
  fEdges->NonManifoldEdgesOff();
  fEdges->ManifoldEdgesOff();

  // strips (polyline + points) hierarchy
  vtkSmartPointer<vtkStripper> stripper =
    vtkSmartPointer<vtkStripper>::New();
  stripper->SetInputConnection( fEdges->GetOutputPort() );
  stripper->Update();

  // change the polylines into boundary polygons polydata bpd
  vtkSmartPointer<vtkPolyData> bpd =
    vtkSmartPointer<vtkPolyData>::New();
  bpd->SetPoints(stripper->GetOutput()->GetPoints());
  bpd->SetPolys(stripper->GetOutput()->GetLines());

  // write
  writeBoundaryPoints( fEdges->GetOutput()->GetPoints(), "pts.vtk" );
  WritePolyData( fEdges->GetOutput(), "fedge.vtk");

  std::cout << "-- boundary polydata has " << bpd->GetNumberOfPoints() <<
    " points." << std::endl;
  std::cout << "-- boundary polydata has " << bpd->GetNumberOfCells() <<
    " cells." << std::endl;

  // boundary data id filter
  vtkSmartPointer<vtkIdFilter> bCellIdsFilter =
    vtkSmartPointer<vtkIdFilter>::New();
  bCellIdsFilter->SetInputData( bpd );
  bCellIdsFilter->SetCellIds(true);
  bCellIdsFilter->SetPointIds(false);
  bCellIdsFilter->SetIdsArrayName("bCellIds");
  bCellIdsFilter->Update();
  WriteDataSet( bCellIdsFilter->GetOutput(), "bCellIdsFilter.vtp" );

  vtkSmartPointer<vtkIdFilter> bPointIdsFilter =
    vtkSmartPointer<vtkIdFilter>::New();
  bPointIdsFilter->SetInputConnection(bCellIdsFilter->GetOutputPort());
  bPointIdsFilter->SetCellIds(false);
  bPointIdsFilter->SetPointIds(true);
  bPointIdsFilter->SetIdsArrayName("bPointIds");
  bPointIdsFilter->Update();
  vtkDataSet* bpdWithIds = bPointIdsFilter->GetOutput();
  WriteDataSet( bpdWithIds, "bBothIdsFilter.vtp" );

  vtkIdTypeArray* bCellIds =
    vtkIdTypeArray::SafeDownCast( bpdWithIds->GetCellData()->GetArray("bCellIds") );
  std::cout << "-- boundary cell ids: ";
  for (vtkIdType i = 0; i < bCellIds->GetNumberOfTuples(); i++)
  {
    std::cout << bCellIds->GetValue(i) << ", ";
  }
  std::cout << std::endl;

  vtkIdTypeArray* bPointIds =
    vtkIdTypeArray::SafeDownCast( bpdWithIds->GetPointData()->GetArray("bPointIds") );
/*
  std::cout << "boundary points ids: " << std::endl;
  for(vtkIdType i = 0; i < bPointIds->GetNumberOfTuples(); i++)
  {
    std::cout << bPointIds->GetValue(i) << ", ";
  }
  std::cout << std::endl;
*/

  vtkSmartPointer<vtkPolyData> bpdFinal =
    vtkSmartPointer<vtkPolyData>::New();
  bpdFinal->DeepCopy(bpdWithIds); // ***** adds the data with id instead of bpd

  // ************************************************************** //
  // data is labeled with ids there are two cells with many points.
  // we want to make surfaces between the points of two cells

  vtkSmartPointer<vtkIdList> cell0_PointIds =
    vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cell1_PointIds =
    vtkSmartPointer<vtkIdList>::New();

  vtkIdType cell0_Id = 0;
  vtkIdType cell1_Id = 1;
  bpdFinal->GetCellPoints( cell0_Id, cell0_PointIds );
  bpdFinal->GetCellPoints( cell1_Id, cell1_PointIds );

  std::cout << "-- cell 0 point ids : ";
  for(vtkIdType i = 0; i < cell0_PointIds->GetNumberOfIds(); i++)
  {
    std::cout << cell0_PointIds->GetId(i) << ",";
  }
  std::cout << std::endl << std::endl;

  std::cout << "-- cell 1 point ids : ";
  for(vtkIdType i = 0; i < cell0_PointIds->GetNumberOfIds(); i++)
  {
    std::cout << cell1_PointIds->GetId(i) << ",";
  }
  std::cout << std::endl << std::endl;

  // Create line pairs on two cells and add to lines.
  vtkSmartPointer<vtkCellArray> lines =
    vtkSmartPointer<vtkCellArray>::New();

  for (vtkIdType i = 0; i < 10; i++) {

    vtkSmartPointer<vtkLine> line0 =
      vtkSmartPointer<vtkLine>::New();
    line0->GetPointIds()->SetId(0, cell0_PointIds->GetId(i));
    line0->GetPointIds()->SetId(1, cell0_PointIds->GetId(i+1));
    lines->InsertNextCell(line0);

    vtkSmartPointer<vtkLine> line1 =
      vtkSmartPointer<vtkLine>::New();
    line1->GetPointIds()->SetId(0, cell1_PointIds->GetId(i));
    line1->GetPointIds()->SetId(1, cell1_PointIds->GetId(i+1));
    lines->InsertNextCell(line1);
  }

  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  //polydata->SetPoints(bpdFinal->GetPoints());
  polydata->DeepCopy(bpdFinal);
  polydata->SetLines(lines);

  WritePolyData(polydata, "polyline.vtk");

  vtkSmartPointer<vtkRuledSurfaceFilter> ruledSurfaceFilter =
    vtkSmartPointer<vtkRuledSurfaceFilter>::New();
  ruledSurfaceFilter->SetInputData(polydata);
  ruledSurfaceFilter->SetResolution(20, 20);
  ruledSurfaceFilter->SetRuledModeToResample();
  ruledSurfaceFilter->SetPassLines(true);

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  writer->SetInputConnection(ruledSurfaceFilter->GetOutputPort());

  writer->SetFileName("ruled.vtk");
  writer->Write();

/*
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();

  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  vtkSmartPointer<vtkRenderWindowInteractor> interactor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetRenderWindow(renderWindow);

  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(ruledSurfaceFilter->GetOutputPort());

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(0.89, 0.81, 0.34);

  // Add the actors to the renderer, set the background and size
  renderer->AddActor(actor);
  renderer->SetBackground(.1, .5, .5);

  renderer->ResetCamera();
  renderer->GetActiveCamera()->Azimuth(60);
  renderer->GetActiveCamera()->Elevation(60);
  renderer->GetActiveCamera()->Dolly(1.2);
  renderer->ResetCameraClippingRange();

  renderWindow->Render();
  interactor->Start();
*/

  return EXIT_SUCCESS;
}

// ***************************************
// ************ global functions *********
// ***************************************

void writeBoundaryPoints(vtkSmartPointer<vtkPoints> p, std::string name){
  // write the points on boundary
  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(p);
  vtkSmartPointer<vtkXMLPolyDataWriter> pWriter =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  pWriter->SetFileName(name.c_str());
  pWriter->SetInputData(polydata);
  pWriter->SetDataModeToAscii();
  std::cout << "writing stuff .." << std::endl;
  pWriter->Write();
}

void writeBoundaryPolyData(vtkSmartPointer<vtkPolyData> data, std::string name){
  //sf
  vtkSmartPointer<vtkDataSetSurfaceFilter> sf =
    vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  sf->SetInputData(data);
  sf->Update();
  // stl writer
  vtkSmartPointer<vtkSTLWriter> sw =
    vtkSmartPointer<vtkSTLWriter>::New();
  sw->SetFileName(name.c_str());
  sw->SetInputConnection(sf->GetOutputPort());
  sw->SetFileTypeToBinary();
  std::cout << "writing .. " << std::endl;
  sw->Write();
}

void writeClosedClipSTL(vtkSmartPointer<vtkPolyData> data, std::string name){
  // write the detected boundary edges
  vtkSmartPointer<vtkSTLWriter> sw2 =
    vtkSmartPointer<vtkSTLWriter>::New();
  sw2->SetFileName(name.c_str());
  std::cout << "writing boundary poly data .. " << std::endl;
  sw2->SetInputData(data);
  sw2->Write();
}

void WritePolyData(vtkPolyData* const polyData, const std::string& filename)
{
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(polyData);
#else
    writer->SetInputData(polyData);
#endif
    writer->SetFileName(filename.c_str());
    writer->Write();
}

void WriteDataSet(vtkDataSet* const dataSet, const std::string& filename)
{
    vtkSmartPointer<vtkDataSetWriter> writer =
      vtkSmartPointer<vtkDataSetWriter>::New();
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(dataSet);
#else
    writer->SetInputData(dataSet);
#endif
    writer->SetFileName(filename.c_str());
    writer->Write();
}
