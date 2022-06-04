/* This makes the surface plot of the solution to Laplace's equation.
 * This code was cloned from the VTK example at
 * https://kitware.github.io/vtk-examples/site/Cxx/Plotting/SurfacePlot
 * and then heavily modified for my purposes.
 */
#include "make_vtk_plot.h"

#include <vtkCamera.h>
#include <vtkChartXYZ.h>
#include <vtkContextMouseEvent.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkFloatArray.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPen.h>
#include <vtkPlotSurface.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTable.h>
#include <vtkVector.h>
#include <vtkAutoInit.h> 

VTK_MODULE_INIT(vtkRenderingContextOpenGL2);
VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle)
  
int make_vtk_plot(int Nx, int Ny, float *x, float *y, float *z)
{
  vtkNew<vtkNamedColors> colors;

  vtkNew<vtkChartXYZ> chart;
  chart->SetGeometry(vtkRectf(70.0, 70.0, 730, 530));

  vtkNew<vtkPlotSurface> plot;

  vtkNew<vtkContextView> view;
  view->GetRenderWindow()->SetSize(900, 700);
  view->GetScene()->AddItem(chart);
  view->GetRenderer()->SetBackground(colors->GetColor3d("Grey").GetData());
  
  //view->GetRenderer()->GetActiveCamera()->SetRoll(30);
  view->GetRenderer()->GetActiveCamera()->Elevation(60);  
  //view->GetRenderer()->GetActiveCamera()->Azimuth(0);
  //view->GetRenderer()->GetActiveCamera()->SetViewUp (1.0, 0.0, 0.0);
  //view->GetRenderer()->ResetCamera();  

  // Create a VTK table which holds the data used for plotting.
  vtkNew<vtkTable> table;
  for (vtkIdType i = 0; i < Nx; i++)
  {
    vtkNew<vtkFloatArray> arr;
    table->AddColumn(arr);
  }

  // Copy input solution vector into VTK table.
  table->SetNumberOfRows(Ny);
  for (vtkIdType i = 0; i < Nx; i++) {
    for (vtkIdType j = 0; j < Ny; j++) {
      table->SetValue(i, j, z[LINDEX(Nx, Ny, i, j)]  );
    }
  }

  // Set up the surface plot we wish to visualize and add it to the chart.
  plot->SetInputData(table);
  plot->GetPen()->SetColorF(colors->GetColor3d("Red").GetData());
  chart->AddPlot(plot);

  view->GetRenderWindow()->SetMultiSamples(0);
  view->GetInteractor()->Initialize();
  view->GetRenderWindow()->SetWindowName("Temperature plot");
  view->GetRenderWindow()->Render();

  view->GetInteractor()->Start();

  return EXIT_SUCCESS;
}
