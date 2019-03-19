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
#include <vtkXMLPolyDataWriter.h>

#include <vtkDataSetSurfaceFilter.h>
#include <vtkPLYWriter.h>
#include <vtkPolyDataWriter.h>

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
#include <vtkTriangle.h>
#include <vtkIntersectionPolyDataFilter.h>

#include <vtkDelaunay2D.h>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkVertexGlyphFilter.h>

#include <vtkConvexHull2D.h>

#include <vtkGeometryFilter.h>

#include <vtkContourFilter.h>

#include <vtkContourTriangulator.h>

#include <vtkPointLocator.h>
#include <vtkIncrementalOctreePointLocator.h>

void WritePolyData(vtkPolyData* const polyData, const std::string& filename);
void WriteDataSet(vtkDataSet* const dataSet, const std::string& filename);
void writeBoundaryPoints( vtkSmartPointer<vtkPoints> p, std::string name );
void writeBoundaryPolyData( vtkSmartPointer<vtkPolyData> d , std::string name );
void writeClosedClipSTL( vtkSmartPointer<vtkPolyData> d, std::string name );

int main (int argc, char *argv[])
{

  vtkSmartPointer<vtkPolygon> leftPolygon =
    vtkSmartPointer<vtkPolygon>::New();

  vtkSmartPointer<vtkPoints> leftPoints =
    vtkSmartPointer<vtkPoints>::New();
  leftPoints->InsertNextPoint(0.02, 0.03, 0.02);
  leftPoints->InsertNextPoint(0.02, 0.03, 0.04);
  leftPoints->InsertNextPoint(-0.02, 0.03, 0.04);
  leftPoints->InsertNextPoint(-0.02, 0.03, 0.02);
  leftPolygon->GetPoints()->DeepCopy(leftPoints);
  leftPolygon->GetPointIds()->SetNumberOfIds(4);
  leftPolygon->GetPointIds()->SetId(0, 0);
  leftPolygon->GetPointIds()->SetId(1, 1);
  leftPolygon->GetPointIds()->SetId(2, 2);
  leftPolygon->GetPointIds()->SetId(3, 3);

  // plays the role of topology
  vtkSmartPointer<vtkCellArray> leftPolygons =
    vtkSmartPointer<vtkCellArray>::New();
  leftPolygons->InsertNextCell(leftPolygon);

  vtkSmartPointer<vtkPolyData> leftPoly =
    vtkSmartPointer<vtkPolyData>::New();
  // Add the geometry and topology to the polydata
  leftPoly->SetPoints ( leftPoints );
  leftPoly->SetPolys ( leftPolygons );

  vtkSmartPointer<vtkTriangleFilter> leftPolyTriangulated =
    vtkSmartPointer<vtkTriangleFilter>::New();
  leftPolyTriangulated->SetInputData(leftPoly);
  leftPolyTriangulated->Update();

  vtkSmartPointer<vtkPLYWriter> plyWriter0 =
    vtkSmartPointer<vtkPLYWriter>::New();
  plyWriter0->SetFileName("leftRectangle.ply");
  plyWriter0->SetInputData(leftPolyTriangulated->GetOutput());
  plyWriter0->Write();

  // PolyData to process
  std::string input_name(argv[1]);
  std::cout << "Reading stl file : " << input_name << std::endl;
  vtkSmartPointer<vtkSTLReader> stlReader =
    vtkSmartPointer<vtkSTLReader>::New();
  stlReader->SetFileName(input_name.c_str());
  stlReader->Update();;
  vtkSmartPointer<vtkPolyData> mainPoly;
  mainPoly = stlReader->GetOutput(); //or try shallowcopy

  // clean polyData
  vtkSmartPointer<vtkCleanPolyData> cpd =
    vtkSmartPointer<vtkCleanPolyData>::New();
  cpd->SetInputData(mainPoly);
  cpd->Update();
//  cpd->Print(std::cout);

  // calculating intersection
  vtkSmartPointer<vtkIntersectionPolyDataFilter> leftIntersection =
   vtkSmartPointer<vtkIntersectionPolyDataFilter>::New();
  leftIntersection->SetInputData( 0, leftPolyTriangulated->GetOutput() );
  leftIntersection->SetInputData( 1, cpd->GetOutput() );
  leftIntersection->Update();
//  leftIntersection->Print(std::cout);

  // ******************************************************************
  // Add ids to the points and cells of the sphere
  vtkSmartPointer<vtkIdFilter> iIdsFilter = vtkSmartPointer<vtkIdFilter>::New();
  iIdsFilter->SetInputConnection(leftIntersection->GetOutputPort());
  iIdsFilter->SetCellIds(true);
  iIdsFilter->SetPointIds(false);
  iIdsFilter->SetIdsArrayName("iIds");
  iIdsFilter->Update();

  vtkDataSet* iIds = iIdsFilter->GetOutput();
  WriteDataSet( iIdsFilter->GetOutput(), "iIdFilter.vtp" );
  //iIdsFilter->Print( std::cout );

  vtkSmartPointer<vtkIdFilter> ipointIdsFilter = vtkSmartPointer<vtkIdFilter>::New();
  ipointIdsFilter->SetInputConnection(iIdsFilter->GetOutputPort());
  ipointIdsFilter->SetCellIds(false);
  ipointIdsFilter->SetPointIds(true);
  ipointIdsFilter->SetIdsArrayName("iPointIds");
  ipointIdsFilter->Update();

  vtkDataSet* iWithIds = ipointIdsFilter->GetOutput();

  WriteDataSet( iWithIds, "iBothIdsData.vtp" );
  //iWithIds->Print(std::cout);
  // ******************************************************************

  vtkSmartPointer<vtkPoints> testPoints =
    vtkSmartPointer<vtkPoints>::New();
  testPoints->DeepCopy( leftIntersection->GetOutput()->GetPoints() );
//  testPoints->Print(std::cout);

  vtkSmartPointer<vtkPointLocator> pointLocator =
    vtkSmartPointer<vtkPointLocator>::New();
//  pointLocator->AutomaticOn();
//  pointLocator->SetNumberOfPointsPerBucket(2);
  pointLocator->SetDataSet( testPoints );
  pointLocator->BuildLocator();

/*
  double bounds[6] = { -0.03, 0.03, 0.02, 0.06, 0.01, 0.05 };
  pointLocator->InitPointInsertion( testPoints, bounds, 250 );

  for (vtkIdType id = 0; id < testPoints->GetNumberOfPoints(); id++ ) {
    double p[3];
    testPoints->GetPoint(id, p);
    pointLocator->InsertPoint(id, p);
  }
*/
  pointLocator->Print(std::cout);

  vtkSmartPointer<vtkPointLocator> ordered =
    vtkSmartPointer<vtkPointLocator>::New();

  // start with the first points in non-ordered list
  for (vtkIdType id = 0; id < testPoints->GetNumberOfPoints(); id++ ) {
    double x[3];
    x[0] = testPoints->GetPoint(id)[0];
    x[1] = testPoints->GetPoint(id)[1];
    x[2] = testPoints->GetPoint(id)[2];
    vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();
    pointLocator->FindClosestNPoints(4, x, list);
    int unique = 0;
    vtkIdType nclose = 0;
    std::cout << "point : " << x[0] <<" " << x[1] <<" " << x[2] << std::endl;
    while ( unique != 1 ) {
      double p[3];
      testPoints->GetPoint( list->GetId(nclose), p);
      std::cout << "closest point : " << p[0] << " " << p[1] << " " << p[2] << std::endl;
      vtkIdType uid;
      unique = ordered->InsertUniquePoint( p, uid);
      nclose ++;
    }
  }

  vtkSmartPointer<vtkPoints> leftPoints2 =
    vtkSmartPointer<vtkPoints>::New();

  leftPoints2->DeepCopy(ordered->GetPoints());

  vtkSmartPointer<vtkCellArray> leftVertices =
    vtkSmartPointer<vtkCellArray>::New();

  vtkSmartPointer<vtkCellArray> leftLines =
    vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkCellArray> aCellArray =
            vtkSmartPointer<vtkCellArray>::New();

  aCellArray->InsertNextCell(leftIntersection->GetOutput()->GetPoints()->GetNumberOfPoints());

  for ( vtkIdType i = 0; i < leftIntersection->GetOutput()->GetPoints()->GetNumberOfPoints(); i++ )
  {
    vtkIdType pointId =
      leftPoints2->InsertNextPoint(leftIntersection->GetOutput()->GetPoints()->GetPoint(i));
    aCellArray->InsertCellPoint(pointId);
  }

  for ( vtkIdType i = 0; i < leftIntersection->GetOutput()->GetPoints()->GetNumberOfPoints(); i++ )
  {
    vtkSmartPointer<vtkLine> line =
      vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId ( 0, i );
    line->GetPointIds()->SetId ( 1, i+1 );
    leftLines->InsertNextCell ( line );
  }

  // Create a polydata to store the polys.
  vtkSmartPointer<vtkPolyData> davidPolys =
          vtkSmartPointer<vtkPolyData>::New();
  davidPolys->SetPoints(leftPoints2);
  davidPolys->SetPolys(aCellArray);

  // Create a polydata to store the contour.
        vtkSmartPointer<vtkPolyData> contour =
                vtkSmartPointer<vtkPolyData>::New();
        contour->SetPoints(leftPoints2);
        contour->SetLines(aCellArray);
        contour->BuildCells();
        contour->BuildLinks();

        // Triangulate the grid points
        vtkSmartPointer<vtkContourTriangulator> davidTriangulator =
                vtkSmartPointer<vtkContourTriangulator>::New();
#if VTK_MAJOR_VERSION <= 5
        davidTriangulator->SetInput(contour);
#else
        davidTriangulator->SetInputData(contour);
#endif
        davidTriangulator->Update();

  // Visualize
  vtkSmartPointer<vtkPolyDataMapper> davidmeshMapper =
          vtkSmartPointer<vtkPolyDataMapper>::New();
  davidmeshMapper->SetInputConnection(davidTriangulator->GetOutputPort());

  vtkSmartPointer<vtkActor> davidmeshActor =
          vtkSmartPointer<vtkActor>::New();
  davidmeshActor->SetMapper(davidmeshMapper);
  davidmeshActor->GetProperty()->SetRepresentationToWireframe();

  vtkSmartPointer<vtkFeatureEdges> davidedgeFilter =
          vtkSmartPointer<vtkFeatureEdges>::New();
#if VTK_MAJOR_VERSION <= 5
  davidedgeFilter->SetInputConnection(polys->GetProducerPort());
#else
  davidedgeFilter->SetInputData(davidPolys);
#endif
  davidedgeFilter->Update();

  vtkSmartPointer<vtkPolyDataMapper> davidcontourMapper =
          vtkSmartPointer<vtkPolyDataMapper>::New();
  davidcontourMapper->SetInputConnection(davidedgeFilter->GetOutputPort());

  vtkSmartPointer<vtkActor> davidcontourActor =
          vtkSmartPointer<vtkActor>::New();
  davidcontourActor->SetMapper(davidcontourMapper);
  davidcontourActor->GetProperty()->SetColor(1,0,0);

  // Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> davidrenderer =
          vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> davidrenderWindow =
          vtkSmartPointer<vtkRenderWindow>::New();
  davidrenderWindow->AddRenderer(davidrenderer);
  vtkSmartPointer<vtkRenderWindowInteractor> davidrenderWindowInteractor =
          vtkSmartPointer<vtkRenderWindowInteractor>::New();
  davidrenderWindowInteractor->SetRenderWindow(davidrenderWindow);

  // Add the actor to the scene
  davidrenderer->AddActor(davidmeshActor);
  davidrenderer->AddActor(davidcontourActor);
  davidrenderer->SetBackground(.3, .6, .3); // Background color green

  // Render and interact
  davidrenderWindow->Render();
  davidrenderWindowInteractor->Start();


 // ***********************************************************************************

  vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
  line->GetPointIds()->SetId ( 0, leftIntersection->GetOutput()->GetPoints()->GetNumberOfPoints());
  line->GetPointIds()->SetId ( 1, 0 );
  leftLines->InsertNextCell ( line );

  vtkSmartPointer<vtkPolyData> left =
    vtkSmartPointer<vtkPolyData>::New();
//  left->DeepCopy(leftIntersection->GetOutput());
  left->SetPoints( leftPoints2 );
  left->SetVerts( leftVertices );
//  left->SetLines( leftLines );

  WritePolyData(left, "left.vtk");

  // *********************************************************************************
/*
  vtkSmartPointer<vtkConvexHull2D> conv =
    vtkSmartPointer<vtkConvexHull2D>::New();
  conv->SetInputData(left);
  conv->SetMinHullSizeInWorld(0.0001);
  conv->SetMinHullSizeInDisplay(0.0001);
  conv->SetHullShape(0);
  conv->OutlineOn();
  vtkSmartPointer<vtkPoints> outPoints =
    vtkSmartPointer<vtkPoints>::New();
  double minimumHullSize=0.1;
  conv->CalculateConvexHull(leftPoints2, outPoints, minimumHullSize);
  conv->Update();
  WritePolyData(conv->GetOutput(), "conv.vtk");
  writeBoundaryPoints( outPoints, "convPts.vtk" );
*/
  // *********************************************************************************

  // *********************************************************************************
  vtkSmartPointer<vtkSurfaceReconstructionFilter> sf =
    vtkSmartPointer<vtkSurfaceReconstructionFilter>::New();
  sf->SetInputData(left);
  sf->SetNeighborhoodSize(5);
  sf->SetSampleSpacing 	(0.0001);
  sf->Update();
  vtkSmartPointer<vtkContourFilter> cf =
    vtkSmartPointer<vtkContourFilter>::New();
  cf->SetInputConnection(sf->GetOutputPort());
  cf->SetValue(0, 0.0);
  vtkSmartPointer<vtkReverseSense> reverse =
    vtkSmartPointer<vtkReverseSense>::New();
  reverse->SetInputConnection(cf->GetOutputPort());
  reverse->ReverseCellsOn();
  reverse->ReverseNormalsOn();
  vtkSmartPointer<vtkPolyData> tmp =
    vtkSmartPointer<vtkPolyData>::New();
  tmp->DeepCopy(cf->GetOutput());
  WritePolyData(tmp, "sf.ply");
  // *********************************************************************************

  // *********************************************************************************
  vtkSmartPointer<vtkDelaunay2D> delaunay =
   vtkSmartPointer<vtkDelaunay2D>::New();

  delaunay->SetInputData( left );
  delaunay->SetAlpha(0.002);
  delaunay->SetTolerance(1e-09);
  delaunay->SetProjectionPlaneMode( VTK_BEST_FITTING_PLANE );
  delaunay->SetOffset(1.0);
  //delaunay->SetSourceData( mainPoly );
  delaunay->Update();

  vtkSmartPointer<vtkContourTriangulator> contourTri =
   vtkSmartPointer<vtkContourTriangulator>::New();
  contourTri->SetInputData(left);
  contourTri->Update();
  WritePolyData(contourTri->GetOutput(), "contTri.ply");

// **************************************************************************************

// Visualize
vtkSmartPointer<vtkPolyDataMapper> meshMapper =
  vtkSmartPointer<vtkPolyDataMapper>::New();
meshMapper->SetInputConnection( contourTri->GetOutputPort());

vtkSmartPointer<vtkActor> meshActor =
  vtkSmartPointer<vtkActor>::New();
meshActor->SetMapper(meshMapper);

vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
  vtkSmartPointer<vtkVertexGlyphFilter>::New();
glyphFilter->SetInputData(left);
glyphFilter->Update();

vtkSmartPointer<vtkPolyDataMapper> pointMapper =
  vtkSmartPointer<vtkPolyDataMapper>::New();
pointMapper->SetInputConnection(glyphFilter->GetOutputPort());

vtkSmartPointer<vtkActor> pointActor =
  vtkSmartPointer<vtkActor>::New();
pointActor->GetProperty()->SetColor(1,0,0);
pointActor->GetProperty()->SetPointSize(3);
pointActor->SetMapper(pointMapper);

vtkSmartPointer<vtkRenderer> renderer =
  vtkSmartPointer<vtkRenderer>::New();
vtkSmartPointer<vtkRenderWindow> renderWindow =
  vtkSmartPointer<vtkRenderWindow>::New();
renderWindow->AddRenderer(renderer);
vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
  vtkSmartPointer<vtkRenderWindowInteractor>::New();
renderWindowInteractor->SetRenderWindow(renderWindow);

renderer->AddActor(meshActor);
renderer->AddActor(pointActor);
renderer->SetBackground(.3, .6, .3); // Background color green

renderWindow->Render();
renderWindowInteractor->Start();

// **************************************************************************************

  // Add ids to the points and cells of the sphere
  vtkSmartPointer<vtkIdFilter> cellIdsFilter = vtkSmartPointer<vtkIdFilter>::New();
  cellIdsFilter->SetInputConnection(cpd->GetOutputPort());
  cellIdsFilter->SetCellIds(true);
  cellIdsFilter->SetPointIds(false);
  cellIdsFilter->SetIdsArrayName("CellIds");
  cellIdsFilter->Update();

  WriteDataSet( cellIdsFilter->GetOutput(), "CellIdFilter.vtp" );

  vtkSmartPointer<vtkIdFilter> pointIdsFilter = vtkSmartPointer<vtkIdFilter>::New();
  pointIdsFilter->SetInputConnection(cellIdsFilter->GetOutputPort());
  pointIdsFilter->SetCellIds(false);
  pointIdsFilter->SetPointIds(true);
  pointIdsFilter->SetIdsArrayName("PointIds");
  pointIdsFilter->Update();

  vtkDataSet* cpdWithIds = pointIdsFilter->GetOutput();

  WriteDataSet( cpdWithIds, "BothIdsData.vtp" );

  // clip box
  vtkSmartPointer<vtkBox> box = vtkSmartPointer<vtkBox>::New();
  box->SetBounds(-0.02, 0.02, 0.03, 0.05, 0.02, 0.04);

  // vtk clip polydata
  vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
  clipper->SetInputData( cpdWithIds ); // ***** adds the data with id instead of cpd
  clipper->SetClipFunction(box);
  clipper->SetValue(0);
  clipper->Update();

  WriteDataSet( clipper->GetOutput(), "clipper.vtp" );

  vtkPolyData* clipped = clipper->GetOutput();

  std::cout << "clipped has " << clipped->GetNumberOfPoints() << " points." << std::endl;
  std::cout << "Clipped has " << clipped->GetNumberOfCells() << " cells." << std::endl;

  vtkIdTypeArray* clippedCellIds =
    vtkIdTypeArray::SafeDownCast(clipped->GetCellData()->GetArray("CellIds"));

//  std::cout << "Clipped cell (id, value): " << std::endl;
/*
  for(vtkIdType i = 0; i < clippedCellIds->GetNumberOfTuples(); i++)
  {
    std::cout << "(" << i << ", " << clippedCellIds->GetValue(i) << "), ";
  }
  std::cout << std::endl;
*/

  // ********************************************* //
  // extract feature edges
  vtkSmartPointer<vtkFeatureEdges> featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
  featureEdges->SetInputData( clipped ); // input is already labled
  featureEdges->BoundaryEdgesOn();
  featureEdges->FeatureEdgesOff();
  featureEdges->NonManifoldEdgesOff();
  featureEdges->ManifoldEdgesOff();

  // strips (polyline + points) hierarchy
  vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
  stripper->SetInputConnection( featureEdges->GetOutputPort() );
  stripper->Update();

  // change the polylines into boundary polygons polydata bpd
  vtkSmartPointer<vtkPolyData> bpd = vtkSmartPointer<vtkPolyData>::New();
  bpd->SetPoints(stripper->GetOutput()->GetPoints());
  bpd->SetPolys(stripper->GetOutput()->GetLines());

  writeBoundaryPoints( featureEdges->GetOutput()->GetPoints(), "pts.vtk" );
  WritePolyData( featureEdges->GetOutput(), "fedge.vtk");
//  WritePolyData( strips->GetOutput(), "strips.vtk");
//  WritePolyData( featureEdges->GetOutput(), "featureEdges->GetOutput().vtk");
//  writeClosedClipSTL( featureEdges->GetOutput(), "boundary.stl" );

  std::cout << "boundary polydata has " << bpd->GetNumberOfPoints() << " points." << std::endl;
  std::cout << "boundary polydata has " << bpd->GetNumberOfCells() << " cells." << std::endl;

  vtkSmartPointer<vtkIdFilter> bCellIdsFilter = vtkSmartPointer<vtkIdFilter>::New();
  bCellIdsFilter->SetInputData( bpd ); // we are trying to create a labeled bpd
  bCellIdsFilter->SetCellIds(true);
  bCellIdsFilter->SetPointIds(false);
  bCellIdsFilter->SetIdsArrayName("bCellIds");
  bCellIdsFilter->Update();

  WriteDataSet( bCellIdsFilter->GetOutput(), "bCellIdsFilter.vtp" );

  vtkSmartPointer<vtkIdFilter> bPointIdsFilter = vtkSmartPointer<vtkIdFilter>::New();
  bPointIdsFilter->SetInputConnection(bCellIdsFilter->GetOutputPort());
  bPointIdsFilter->SetCellIds(false);
  bPointIdsFilter->SetPointIds(true);
  bPointIdsFilter->SetIdsArrayName("bPointIds");
  bPointIdsFilter->Update();

  vtkDataSet* bpdWithIds =
    bPointIdsFilter->GetOutput();
  WriteDataSet( bpdWithIds, "bBothIdsFilter.vtp" );

  vtkIdTypeArray* bCellIds = vtkIdTypeArray::SafeDownCast( bpdWithIds->GetCellData()->GetArray("bCellIds") );
  /*
  std::cout << "boundary cell ids: " << std::endl;
  for(vtkIdType i = 0; i < bCellIds->GetNumberOfTuples(); i++)
  {
    std::cout << bCellIds->GetValue(i) << ", ";
  }
  std::cout << std::endl;
  */
  vtkIdTypeArray* bPointIds = vtkIdTypeArray::SafeDownCast( bpdWithIds->GetPointData()->GetArray("bPointIds") );
/*
  std::cout << "boundary points ids: " << std::endl;
  for(vtkIdType i = 0; i < bPointIds->GetNumberOfTuples(); i++)
  {
    std::cout << bPointIds->GetValue(i) << ", ";
  }
  std::cout << std::endl;
*/
  vtkSmartPointer<vtkPolyData> bpdFinal = vtkSmartPointer<vtkPolyData>::New();
  bpdFinal->DeepCopy(bpdWithIds); // ***** adds the data with id instead of bpd

  // plays the role of topology
  vtkSmartPointer<vtkCellArray> boxPolygons =
    vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPoints> boxPoints =
      vtkSmartPointer<vtkPoints>::New();
  boxPoints->InsertNextPoint(0.02, 0.03, 0.02);
  boxPoints->InsertNextPoint(0.02, 0.05, 0.02);
  boxPoints->InsertNextPoint(-0.02, 0.05, 0.02);
  boxPoints->InsertNextPoint(-0.02, 0.03, 0.02);
  boxPoints->InsertNextPoint(0.02, 0.03, 0.04);
  boxPoints->InsertNextPoint(0.02, 0.05, 0.04);
  boxPoints->InsertNextPoint(-0.02, 0.05, 0.04);
  boxPoints->InsertNextPoint(-0.02, 0.03, 0.04);

  vtkSmartPointer<vtkPolygon> newleftPolygon =
    vtkSmartPointer<vtkPolygon>::New();
  newleftPolygon->GetPointIds()->SetNumberOfIds(4);
  newleftPolygon->GetPointIds()->SetId(0, 3);
  newleftPolygon->GetPointIds()->SetId(1, 0);
  newleftPolygon->GetPointIds()->SetId(2, 4);
  newleftPolygon->GetPointIds()->SetId(3, 7);

  vtkSmartPointer<vtkPolygon> topPolygon =
    vtkSmartPointer<vtkPolygon>::New();
  topPolygon->GetPointIds()->SetNumberOfIds(4);
  topPolygon->GetPointIds()->SetId(0, 4);
  topPolygon->GetPointIds()->SetId(1, 5);
  topPolygon->GetPointIds()->SetId(2, 6);
  topPolygon->GetPointIds()->SetId(3, 7);

  vtkSmartPointer<vtkPolygon> rightPolygon =
    vtkSmartPointer<vtkPolygon>::New();
  rightPolygon->GetPointIds()->SetNumberOfIds(4);
  rightPolygon->GetPointIds()->SetId(0, 2);
  rightPolygon->GetPointIds()->SetId(1, 1);
  rightPolygon->GetPointIds()->SetId(2, 5);
  rightPolygon->GetPointIds()->SetId(3, 6);

  vtkSmartPointer<vtkPolygon> bottomPolygon =
    vtkSmartPointer<vtkPolygon>::New();
  bottomPolygon->GetPointIds()->SetNumberOfIds(4);
  bottomPolygon->GetPointIds()->SetId(0, 0);
  bottomPolygon->GetPointIds()->SetId(1, 1);
  bottomPolygon->GetPointIds()->SetId(2, 2);
  bottomPolygon->GetPointIds()->SetId(3, 3);

  boxPolygons->InsertNextCell(newleftPolygon);
  boxPolygons->InsertNextCell(topPolygon);
  boxPolygons->InsertNextCell(bottomPolygon);
  boxPolygons->InsertNextCell(rightPolygon);

  vtkSmartPointer<vtkPolyData> boxPoly =
    vtkSmartPointer<vtkPolyData>::New();
  // Add the geometry and topology to the polydata
  boxPoly->SetPoints ( boxPoints );
  boxPoly->SetPolys ( boxPolygons );

/*
  vtkSmartPointer<vtkTriangleFilter> boxTriangulated =
    vtkSmartPointer<vtkTriangleFilter>::New();
  boxTriangulated->SetInputData( boxPoly);
  boxTriangulated->Update();

  WritePolyData(boxTriangulated->GetOutput(), "boxtri.ply");
*/

  vtkSmartPointer<vtkDelaunay2D> bdelaunay =
   vtkSmartPointer<vtkDelaunay2D>::New();

  bdelaunay->SetInputData( bpdFinal );
  bdelaunay->SetAlpha(0.003);
  bdelaunay->SetTolerance(1e-09);
  bdelaunay->SetProjectionPlaneMode( VTK_BEST_FITTING_PLANE );
  bdelaunay->SetOffset(1.0);
  delaunay->SetSourceData( bpdFinal );

  bdelaunay->Update();

  WritePolyData(bdelaunay->GetOutput(), "bdel.ply");
 // *************************************************************** //

  vtkSmartPointer<vtkContourTriangulator> bcontourTri =
    vtkSmartPointer<vtkContourTriangulator>::New();
  bcontourTri->SetInputData(bpdWithIds);
  bcontourTri->Update();
  WritePolyData(bcontourTri->GetOutput(), "bcontTri.ply");

  // ************************************************************** //

  vtkSmartPointer<vtkSurfaceReconstructionFilter> bsf =
    vtkSmartPointer<vtkSurfaceReconstructionFilter>::New();
  bsf->SetInputData(bpdFinal);
  bsf->SetNeighborhoodSize(20);
  bsf->SetSampleSpacing 	(0.0001);
  bsf->Update();
  vtkSmartPointer<vtkContourFilter> bcf =
    vtkSmartPointer<vtkContourFilter>::New();
  bcf->SetInputConnection(sf->GetOutputPort());
  bcf->SetValue(0, 0.0);
  vtkSmartPointer<vtkReverseSense> breverse =
    vtkSmartPointer<vtkReverseSense>::New();
  breverse->SetInputConnection(cf->GetOutputPort());
  breverse->ReverseCellsOn();
  breverse->ReverseNormalsOn();
  vtkSmartPointer<vtkPolyData> btmp =
    vtkSmartPointer<vtkPolyData>::New();
  btmp->DeepCopy(bcf->GetOutput());
  WritePolyData(btmp, "bsf.ply");

  // ************************************************************** //

  // Visualize
  vtkSmartPointer<vtkPolyDataMapper> bmeshMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  bmeshMapper->SetInputConnection(bdelaunay->GetOutputPort());

  vtkSmartPointer<vtkActor> bmeshActor =
    vtkSmartPointer<vtkActor>::New();
  bmeshActor->SetMapper(bmeshMapper);

  vtkSmartPointer<vtkVertexGlyphFilter> bglyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  bglyphFilter->SetInputData(bpdFinal);
  bglyphFilter->Update();

  vtkSmartPointer<vtkPolyDataMapper> bpointMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  bpointMapper->SetInputConnection(bglyphFilter->GetOutputPort());

  vtkSmartPointer<vtkActor> bpointActor =
    vtkSmartPointer<vtkActor>::New();
  bpointActor->GetProperty()->SetColor(1,0,0);
  bpointActor->GetProperty()->SetPointSize(3);
  bpointActor->SetMapper(bpointMapper);

  vtkSmartPointer<vtkRenderer> brenderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> brenderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  brenderWindow->AddRenderer(brenderer);
  vtkSmartPointer<vtkRenderWindowInteractor> brenderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  brenderWindowInteractor->SetRenderWindow(brenderWindow);

  brenderer->AddActor(bmeshActor);
  brenderer->AddActor(bpointActor);
  brenderer->SetBackground(.3, .6, .3); // Background color green

  brenderWindow->Render();
  brenderWindowInteractor->Start();

  // ************************************************************** //
  // data is labeled with ids there are two cells with many points.
  // we want to make surfaces between the points of two cells

  vtkSmartPointer<vtkIdList> cell0_PointIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cell1_PointIds = vtkSmartPointer<vtkIdList>::New();

  vtkIdType cell0_Id = 0;
  vtkIdType cell1_Id = 1;
  bpdFinal->GetCellPoints( cell0_Id, cell0_PointIds );
  bpdFinal->GetCellPoints( cell1_Id, cell1_PointIds );

/*
  std::cout << "cell 0 point ids : " << std::endl;
  for(vtkIdType i = 0; i < cell0_PointIds->GetNumberOfIds(); i++)
  {
    std::cout << cell0_PointIds->GetId(i) << ",";
  }
  std::cout << std::endl << std::endl;

  std::cout << "cell 1 point ids : " << std::endl;
  for(vtkIdType i = 0; i < cell0_PointIds->GetNumberOfIds(); i++)
  {
    std::cout << cell1_PointIds->GetId(i) << ",";
  }
  std::cout << std::endl << std::endl;
*/

  // Create line pairs on two cells and add to lines.
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

  for (vtkIdType i = 0; i < 400; i++) {

    vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
    line0->GetPointIds()->SetId(0, cell0_PointIds->GetId(i));
    line0->GetPointIds()->SetId(1, cell0_PointIds->GetId(i+1));
    lines->InsertNextCell(line0);

    vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
    line1->GetPointIds()->SetId(0, cell1_PointIds->GetId(i));
    line1->GetPointIds()->SetId(1, cell1_PointIds->GetId(i+1));
    lines->InsertNextCell(line1);
  }

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  //polydata->SetPoints(bpdFinal->GetPoints());
  polydata->DeepCopy(bpdFinal);
  polydata->SetLines(lines);

  WritePolyData(polydata, "polyline.vtk");

  vtkSmartPointer<vtkRuledSurfaceFilter> ruledSurfaceFilter = vtkSmartPointer<vtkRuledSurfaceFilter>::New();
  ruledSurfaceFilter->SetInputData(polydata);
  ruledSurfaceFilter->SetResolution(10, 10);
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
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(p);
  vtkSmartPointer<vtkXMLPolyDataWriter> pWriter =  vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  pWriter->SetFileName(name.c_str());
  pWriter->SetInputData(polydata);
  pWriter->SetDataModeToAscii();
  std::cout << "writing stuff .." << std::endl;
  pWriter->Write();
}

void writeBoundaryPolyData(vtkSmartPointer<vtkPolyData> data, std::string name){
  //sf
  vtkSmartPointer<vtkDataSetSurfaceFilter> sf = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  sf->SetInputData(data);
  sf->Update();
  // stl writer
  vtkSmartPointer<vtkSTLWriter> sw = vtkSmartPointer<vtkSTLWriter>::New();
  sw->SetFileName(name.c_str());
  sw->SetInputConnection(sf->GetOutputPort());
  sw->SetFileTypeToBinary();
  std::cout << "writing .. " << std::endl;
  sw->Write();
}

void writeClosedClipSTL(vtkSmartPointer<vtkPolyData> data, std::string name){
  // write the detected boundary edges
  vtkSmartPointer<vtkSTLWriter> sw2 = vtkSmartPointer<vtkSTLWriter>::New();
  sw2->SetFileName(name.c_str());
  std::cout << "writing boundary poly data .. " << std::endl;
  sw2->SetInputData(data);
  sw2->Write();
}

void WritePolyData(vtkPolyData* const polyData, const std::string& filename)
{
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
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
    vtkSmartPointer<vtkDataSetWriter> writer = vtkSmartPointer<vtkDataSetWriter>::New();
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(dataSet);
#else
    writer->SetInputData(dataSet);
#endif
    writer->SetFileName(filename.c_str());
    writer->Write();
}
