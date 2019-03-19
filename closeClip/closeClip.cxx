/**
***********************************************
author: Milad Kiaee, kiaeedar@ualberta

discription: vtk code to close-clip a surface mesh (STL)
the resulting features are being connected to provide
a closed surtace.
usuage: exec ./closeClip input.stl
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

#include <vtkPointLocator.h>

#include <vtkMath.h>

#include <iostream>
#include <vector>

void WritePolyData(vtkPolyData* const polyData, const std::string& filename);
void WriteDataSet(vtkDataSet* const dataSet, const std::string& filename);
void writeBoundaryPoints( vtkSmartPointer<vtkPoints> p, std::string name );
void writeBoundaryPolyData( vtkSmartPointer<vtkPolyData> d , std::string name );
void writeClosedClipSTL( vtkSmartPointer<vtkPolyData> d, std::string name );
void WriteConnectingSurf(double[3], double[3], double[3], double[3],
  vtkSmartPointer<vtkCleanPolyData>, std::string);

int main (int argc, char *argv[])
{
  // PolyData to process
  std::string input_name(argv[1]);
  std::cout << "Reading stl file : " << input_name << std::endl;
  vtkSmartPointer<vtkSTLReader> stlReader =
    vtkSmartPointer<vtkSTLReader>::New();
  stlReader->SetFileName(input_name.c_str());
  stlReader->Update();;
  vtkSmartPointer<vtkPolyData> mainPoly;
  mainPoly = stlReader->GetOutput(); //or try shallowcopy
  //std::cout << "main poly" << std::endl;
  //mainPoly->Print(std::cout);

  // clean polyData
  vtkSmartPointer<vtkCleanPolyData> cpd =
    vtkSmartPointer<vtkCleanPolyData>::New();
  cpd->SetInputData(mainPoly);
  cpd->Update();

  // *******************************************
  // this can be also be implemented in a function in an automatic manner
  double left_p1[3] = {0.02, 0.03, 0.02};
  double left_p2[3] = {0.02, 0.03, 0.04};
  double left_p3[3] = {-0.02, 0.03, 0.04};
  double left_p4[3] = {-0.02, 0.03, 0.02};

  double top_p1[3] = {0.02, 0.03, 0.04};
  double top_p2[3] = {0.02, 0.05, 0.04};
  double top_p3[3] = {-0.02, 0.05, 0.04};
  double top_p4[3] = {-0.02, 0.03, 0.04};

  double right_p1[3] = {0.02, 0.05, 0.04};
  double right_p2[3] = {0.02, 0.05, 0.02};
  double right_p3[3] = {-0.02, 0.05, 0.02};
  double right_p4[3] = {-0.02, 0.05, 0.04};

  double bottom_p1[3] = {-0.02, 0.03, 0.02};
  double bottom_p2[3] = {-0.02, 0.05, 0.02};
  double bottom_p3[3] = {0.02, 0.05, 0.02};
  double bottom_p4[3] = {0.02, 0.03, 0.02};

  // *******************************************

  WriteConnectingSurf( left_p1, left_p2, left_p3, left_p4, cpd, "left");
  WriteConnectingSurf( top_p1, top_p2, top_p3, top_p4, cpd, "top");
  WriteConnectingSurf( right_p1, right_p2, right_p3, right_p4, cpd, "right");
  WriteConnectingSurf( bottom_p1, bottom_p2, bottom_p3, bottom_p4, cpd, "bottom");

  // **************************************************
  // ********************************************************

  // Add ids to the points and cells of the sphere
  vtkSmartPointer<vtkIdFilter> cellIdsFilter = vtkSmartPointer<vtkIdFilter>::New();
  cellIdsFilter->SetInputConnection(cpd->GetOutputPort());
  cellIdsFilter->SetCellIds(true);
  cellIdsFilter->SetPointIds(false);
  cellIdsFilter->SetIdsArrayName("CellIds");
  cellIdsFilter->Update();

  //WriteDataSet( cellIdsFilter->GetOutput(), "CellIdFilter.vtp" );

  vtkSmartPointer<vtkIdFilter> pointIdsFilter = vtkSmartPointer<vtkIdFilter>::New();
  pointIdsFilter->SetInputConnection(cellIdsFilter->GetOutputPort());
  pointIdsFilter->SetCellIds(false);
  pointIdsFilter->SetPointIds(true);
  pointIdsFilter->SetIdsArrayName("PointIds");
  pointIdsFilter->Update();

  vtkDataSet* cpdWithIds = pointIdsFilter->GetOutput();

  //WriteDataSet( cpdWithIds, "BothIdsData.vtp" );

  // clip box
  vtkSmartPointer<vtkBox> box = vtkSmartPointer<vtkBox>::New();
  box->SetBounds(-0.02, 0.02, 0.03, 0.05, 0.02, 0.04);

  // vtk clip polydata
  vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
  clipper->SetInputData( cpdWithIds ); // ***** adds the data with id instead of cpd
  clipper->SetClipFunction(box);
  clipper->SetValue(0);
  clipper->Update();

//  WriteDataSet( clipper->GetOutput(), "clipper.vtp" );

  writeClosedClipSTL(clipper->GetOutput(), "clipper.stl");

  vtkPolyData* clipped = clipper->GetOutput();

  std::cout << "clipped has " << clipped->GetNumberOfPoints() << " points." << std::endl;
  std::cout << "clipped has " << clipped->GetNumberOfCells() << " cells." << std::endl;

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

void WriteConnectingSurf(double p1[3], double p2[3], double p3[3], double p4[3],
  vtkSmartPointer<vtkCleanPolyData> cpd, std::string input){

  std::cout << "finding " << input << " .." << std::endl;

  vtkSmartPointer<vtkPolygon> polygon =
    vtkSmartPointer<vtkPolygon>::New();

  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  points->InsertNextPoint( p1[0], p1[1], p1[2] );
  points->InsertNextPoint( p2[0], p2[1], p2[2] );
  points->InsertNextPoint( p3[0], p3[1], p3[2] );
  points->InsertNextPoint( p4[0], p4[1], p4[2] );
  polygon->GetPoints()->DeepCopy(points);
  polygon->GetPointIds()->SetNumberOfIds(4);
  polygon->GetPointIds()->SetId(0, 0);
  polygon->GetPointIds()->SetId(1, 1);
  polygon->GetPointIds()->SetId(2, 2);
  polygon->GetPointIds()->SetId(3, 3);

  // plays the role of topology
  vtkSmartPointer<vtkCellArray> polygons =
    vtkSmartPointer<vtkCellArray>::New();
  polygons->InsertNextCell(polygon);

  vtkSmartPointer<vtkPolyData> poly =
    vtkSmartPointer<vtkPolyData>::New();
  // Add the geometry and topology to the polydata
  poly->SetPoints ( points );
  poly->SetPolys ( polygons );

  vtkSmartPointer<vtkTriangleFilter> polyTriangulated =
    vtkSmartPointer<vtkTriangleFilter>::New();
  polyTriangulated->SetInputData(poly);
  polyTriangulated->Update();

  // *******************************************

  // calculating intersection
  vtkSmartPointer<vtkIntersectionPolyDataFilter> intersection =
   vtkSmartPointer<vtkIntersectionPolyDataFilter>::New();
  intersection->SetInputData( 0, polyTriangulated->GetOutput() );
  intersection->SetInputData( 1, cpd->GetOutput() );
  intersection->Update();

  std::string st = input + ".ply";
  WritePolyData(intersection->GetOutput(), st.c_str());

  vtkSmartPointer<vtkPoints> testPoints =
    vtkSmartPointer<vtkPoints>::New();
  testPoints->DeepCopy( intersection->GetOutput()->GetPoints() );

  vtkSmartPointer<vtkPolyData> testPolyData =
    vtkSmartPointer<vtkPolyData>::New();
  testPolyData->SetPoints(testPoints);

  vtkSmartPointer<vtkPointLocator> pointLocator =
    vtkSmartPointer<vtkPointLocator>::New();
  pointLocator->SetDataSet( testPolyData );

  double bounds[6] = { -0.02, 0.02, 0.03, 0.05, 0.02, 0.04 };
  pointLocator->InitPointInsertion( testPoints, bounds);
  pointLocator->BuildLocator();

  vtkSmartPointer<vtkIdList> order = vtkSmartPointer<vtkIdList>::New();
  vtkIdType first = 0;
  order->InsertId( first, first ); // both ordered and unordered start from 0

  vtkIdType N = 50;

  std::cout << "number of points: " << testPoints->GetNumberOfPoints() << std::endl;
  for (vtkIdType i = first; i < testPoints->GetNumberOfPoints()-1; i++ ) {
    std::cout << "(" << i << ")";
    // current point in order list
    double x[3];
    x[0] = testPoints->GetPoint(order->GetId(i))[0];
    x[1] = testPoints->GetPoint(order->GetId(i))[1];
    x[2] = testPoints->GetPoint(order->GetId(i))[2];

    vtkSmartPointer<vtkIdList> closestList = vtkSmartPointer<vtkIdList>::New();
    pointLocator->FindClosestNPoints( N, x, closestList );
    vtkIdType nclose = 0;
    vtkIdType loc = 10000; //some large number
    while ( loc != i+1 ) {

      if ( i == 0 ) {
        nclose ++;
        loc = order->InsertUniqueId( closestList->GetId(nclose) );
        std::cout << closestList->GetId(nclose) << "--" << loc << "||";

      } else if ( i > 0 ) {
        double prev[3];
        prev[0] = testPoints->GetPoint(order->GetId( i - 1 ))[0];
        prev[1] = testPoints->GetPoint(order->GetId( i - 1 ))[1];
        prev[2] = testPoints->GetPoint(order->GetId( i - 1 ))[2];
        double v1[3];
        v1[0] = x[0] - prev[0];
        v1[1] = x[1] - prev[1];
        v1[2] = x[2] - prev[2];

        double refId = 10000;
        double refTheta = 1000;

        for (vtkIdType j=0; j<N; j++){

          vtkIdType candidId = closestList->GetId( j );
          double candid[3];
          candid[0] = testPoints->GetPoint(candidId)[0];
          candid[1] = testPoints->GetPoint(candidId)[1];
          candid[2] = testPoints->GetPoint(candidId)[2];

          if ( order->IsId( candidId ) == -1 ){
            //std::cout << "-- -- -- candid " << j <<" : " << candid[0] << ", "
            //  << candid[1] << ", " << candid[2] << std::endl;
            double v2[3];
            v2[0] = candid[0] - x[0];
            v2[1] = candid[1] - x[1];
            v2[2] = candid[2] - x[2];

            double v1n = vtkMath::Norm(v1);
            //std::cout << "norm v1 = " << v1n << std::endl;
            double v2n = vtkMath::Norm(v2);
            //std::cout << "norm v2 = " << v2n << std::endl;
            double theta = vtkMath::AngleBetweenVectors( v1, v2 )/vtkMath::Pi();

            if ( v2n/v1n < 4 || theta < 0.1 || v2n < 0.0005) {
              refId = candidId;
              break;
            } else if (theta < refTheta) {
              refId = candidId;
              refTheta = theta;
            }
          }
        }

        loc = order->InsertUniqueId( refId );
        std::cout << refId << "--";

      }
    }
  }
  std::cout << std::endl;

  vtkSmartPointer<vtkCellArray> lines =
    vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkCellArray> aCellArray =
    vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPoints> points2 =
    vtkSmartPointer<vtkPoints>::New();

  aCellArray->InsertNextCell(order->GetNumberOfIds());

  for ( vtkIdType i = 0; i < order->GetNumberOfIds(); i++ )
  {
    vtkIdType id = points2->InsertNextPoint(testPoints->GetPoint(order->GetId(i)));
    aCellArray->InsertCellPoint(id);
  }

  for ( vtkIdType i = 0; i < order->GetNumberOfIds() - 1 ; i++ )
  {
    vtkSmartPointer<vtkLine> line =
      vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId ( 0, i );
    line->GetPointIds()->SetId ( 1, i+1 );
    lines->InsertNextCell ( line );
  }

  vtkSmartPointer<vtkLine> lastline =
    vtkSmartPointer<vtkLine>::New();
  lastline->GetPointIds()->SetId ( 0, order->GetNumberOfIds() - 1 );
  lastline->GetPointIds()->SetId ( 1, 0 );
  lines->InsertNextCell ( lastline );

  // Create a polydata to store the polys.
  vtkSmartPointer<vtkPolyData> poly2 =
    vtkSmartPointer<vtkPolyData>::New();
  poly2->SetPoints(points2);
  poly2->SetPolys(aCellArray);
  poly2->SetLines(lines);

  vtkSmartPointer<vtkTriangleFilter> poly2Tri =
    vtkSmartPointer<vtkTriangleFilter>::New();
  poly2Tri->SetInputData( poly2 );
  poly2Tri->Update();

  //WritePolyData(poly2, "leftPolys2.ply");
  std::string out_name = input + ".stl";
  writeClosedClipSTL( poly2Tri->GetOutput(), out_name.c_str() );

  std::cout << " - - - - - " << std::endl;
}
